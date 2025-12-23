import time
import psutil
import os
import multiprocessing
import argparse
from collections import Counter, defaultdict

import numpy as np
from scipy.stats import norm
from Bio import SeqIO as BioSeqIO
from Bio.Seq import Seq
from Bio import SeqRecord
import csv
import traceback
import cProfile
import pstats

# Import custom modules
from outpost_detector import OutpostDetectorOptimized
from suffix_tree import SuffixTree
from compute_irreducible_words import IrreducibleWordComputer
from build_assembly_graph import build_assembly_graph_for_eulerian

from linear_unitig_extractor import extract_unitigs_linear     
from prefix_kmp_scaffolder import prefix_kmp_scaffolder          
from graph_guided_scaffolder import graph_guided_scaffolder  

# --- Imports for 3-Step Assembly Process ---
from graph_analyzer import GraphAnalyzer
from contig_filter import filter_non_maximal_contigs
# The new 3-step modules
from unitig_extractor import UnitigExtractor
from unitig_graph_builder import UnitigGraphBuilder
from unitig_graph_resolver import UnitigGraphResolver
from contig_merger import merge_contigs

# --- import for Graph Simplification ---
from graph_collapsing import simplify_graph_by_collapsing

# --- Other (currently unused) assemblers from original file ---
from extract_contigs_eulerian import extract_contigs_eulerian_robust
from all_paths_extractor import AllPathsExtractor
from graph_simplifier import GraphSimplifier

# Required imports (add to your existing imports)
import mmap
import struct
import tempfile
import shutil
import atexit
from concurrent.futures import ProcessPoolExecutor, as_completed

# Constants for index struct
_IDX_STRUCT = "QQ"   # offset (unsigned long long), length (unsigned long long)
_IDX_SIZE = struct.calcsize(_IDX_STRUCT)


# ----------------------------
# 1) Build mmap-backed read store
# ----------------------------
def build_mmap_reads(reads, tmp_prefix=None):
    """
    Pack reads into two files:
      - <tmp_prefix>.bin : concatenated bytes of all reads (ASCII)
      - <tmp_prefix>.idx : index table of (offset, length) entries (unsigned long long)
    Returns: (idx_path, bin_path, num_reads, tmp_dir)
    tmp_dir should be removed by caller when finished (cleanup helper provided).
    """
    if tmp_prefix is None:
        tmp_dir = tempfile.mkdtemp(prefix="assembler_reads_")
        base = tmp_dir + os.sep + "reads"
    else:
        tmp_dir = tempfile.mkdtemp(prefix="assembler_reads_")
        base = tmp_dir + os.sep + tmp_prefix

    bin_path = base + ".bin"
    idx_path = base + ".idx"
    offset = 0

    # Write files in binary mode
    with open(bin_path, "wb") as bf, open(idx_path, "wb") as ix:
        for r in reads:
            b = r.encode("ascii")
            bf.write(b)
            ix.write(struct.pack(_IDX_STRUCT, offset, len(b)))
            offset += len(b)

    return idx_path, bin_path, len(reads), tmp_dir


# ----------------------------
# 2) Worker globals & helper accessor
# ----------------------------
# Worker-process global variables (initialized by init_worker)
_g_idx_mmap = None
_g_bin_mmap = None
_g_num_reads = None
_g_forward_detector = None
_g_backward_detector = None
_g_irreducible_computer = None

def _get_read_from_mmap(i):
    """Return the i-th read as a python string using the worker-local mmap(s)."""
    global _g_idx_mmap, _g_bin_mmap
    off, ln = struct.unpack_from(_IDX_STRUCT, _g_idx_mmap, i * _IDX_SIZE)
    return _g_bin_mmap[off: off + ln].decode("ascii")


# ----------------------------
# 3) Worker initializer (runs once per worker process)
# ----------------------------
def init_worker(idx_path, bin_path, num_reads, debug=False):
    """
    Per-process initializer for ProcessPoolExecutor.
    Memory-maps idx and bin files and constructs per-worker SuffixTree and detectors.
    """
    global _g_idx_mmap, _g_bin_mmap, _g_num_reads
    global _g_forward_detector, _g_backward_detector, _g_irreducible_computer

    _g_num_reads = num_reads

    # Open and mmap index and bin files read-only
    f_idx = open(idx_path, "rb")
    _g_idx_mmap = mmap.mmap(f_idx.fileno(), 0, access=mmap.ACCESS_READ)

    f_bin = open(bin_path, "rb")
    _g_bin_mmap = mmap.mmap(f_bin.fileno(), 0, access=mmap.ACCESS_READ)

    # Build local reads list for suffix-tree construction (worker-local memory)
    reads_local = []
    for i in range(num_reads):
        off, ln = struct.unpack_from(_IDX_STRUCT, _g_idx_mmap, i * _IDX_SIZE)
        reads_local.append(_g_bin_mmap[off: off + ln].decode("ascii"))

    # Build and store the forward and reversed suffix trees and detectors
    reversed_reads_local = [r[::-1] for r in reads_local]

    # Build suffix trees (may be memory heavy per worker; unavoidable unless SuffixTree is serializable)
    st = SuffixTree(reads_local)
    rst = SuffixTree(reversed_reads_local)

    _g_forward_detector = OutpostDetectorOptimized(st)
    _g_backward_detector = OutpostDetectorOptimized(rst)
    _g_irreducible_computer = IrreducibleWordComputer()

    if debug:
        proc = psutil.Process(os.getpid())
        # report worker memory usage to stdout (no queue)
        print(f"[worker {os.getpid()}] initialized. RSS: {proc.memory_info().rss / (1024**2):.2f} MB")


# ----------------------------
# 4) Per-read worker function (top-level, picklable)
# ----------------------------
def _process_read_mmap(i, debug=False):
    """
    Worker task: uses worker-local globals initialized by init_worker.
    Returns: (right_x, left_x, mapped_irreducible_words, irreducible_word_lengths, summary)
    Same return shape as your previous implementation.
    """
    global _g_forward_detector, _g_backward_detector, _g_irreducible_computer

    try:
        read = _get_read_from_mmap(i)
        n = len(read)

        # compute right outposts
        right_x = []
        for j in range(n):
            rx, _ = _g_forward_detector.find_outpost_and_right_optimized(i, j, debug)
            right_x.append(rx)

        # compute left outposts using reversed detectors
        left_x = [0] * n
        for j_original in range(1, n + 1):
            j_reversed = n - j_original
            lx, _ = _g_backward_detector.find_outpost_and_left_optimized(i, j_reversed, debug)
            left_x[j_original - 1] = lx

        # compute irreducible intervals and map to words
        intervals = _g_irreducible_computer.compute_irreducible_words(right_x, left_x)
        mapped = []
        lengths = []
        for start, end in intervals:
            if isinstance(start, int) and isinstance(end, int) and 1 <= start <= n and 1 <= end <= n and start <= end:
                w = read[start - 1: end]
                mapped.append((w, i))
                lengths.append(len(w))
            else:
                mapped.append(("Invalid Indices", i))

        summary = {"right_x": right_x, "left_x": left_x}
        return right_x, left_x, mapped, lengths, summary

    except Exception as e:
        # Return the canonical error result structure
        tb = traceback.format_exc()
        return [], [], [("Error_in_Read_Processing", i)], [], {"error_message": str(e), "traceback": tb}


# ----------------------------
# 5) Orchestration using ProcessPoolExecutor + mmap
# ----------------------------
def compute_irreducible_words_parallel_mmap(reads, args, tmp_prefix=None):
    """
    Build mmap-backed files, spin up a ProcessPoolExecutor with an initializer
    that mmaps files and builds worker-local detectors, and execute _process_read_mmap
    across read indices. Returns the same tuple you used previously.
    """
    # 1) create mmap files (temporary directory returned for cleanup)
    idx_path, bin_path, num_reads, tmp_dir = build_mmap_reads(reads, tmp_prefix)

    # Ensure temp_dir gets cleaned on process exit of main
    def _cleanup_tmp_dir():
        try:
            shutil.rmtree(tmp_dir)
        except Exception:
            pass
    atexit.register(_cleanup_tmp_dir)

    # 2) launch ProcessPoolExecutor
    n_workers = max(1, int(args.threads))
    results = [None] * num_reads

    initargs = (idx_path, bin_path, num_reads, args.debug)

    # Submit tasks — pass only the read index to avoid pickling reads
    with ProcessPoolExecutor(max_workers=n_workers, initializer=init_worker, initargs=initargs) as exe:
        future_to_idx = {exe.submit(_process_read_mmap, i, args.debug): i for i in range(num_reads)}

        for fut in as_completed(future_to_idx):
            idx = future_to_idx[fut]
            try:
                results[idx] = fut.result()
            except Exception as exc:
                # ensure consistent fallback row
                results[idx] = ([], [], [("Error_in_Read_Processing", idx)], [], {"error_message": str(exc), "traceback": traceback.format_exc()})

    # 3) Aggregate results in the same shape as before
    all_right_x = [r[0] for r in results]
    all_left_x = [r[1] for r in results]
    all_mapped_irreducible_words_with_indices = [r[2] for r in results]
    all_irreducible_word_lengths = [length for r in results for length in r[3]]
    irreducible_word_counts_per_read = [len([w for w in r[2] if w[0] not in ("Invalid Indices", "Error_in_Read_Processing")]) for r in results]
    all_outposts_summary = [r[4] for r in results]

    # remove temporary files now that the data has been read into workers (workers have their own mmaps)
    try:
        shutil.rmtree(tmp_dir)
    except Exception:
        pass

    return all_right_x, all_left_x, all_mapped_irreducible_words_with_indices, all_irreducible_word_lengths, \
           irreducible_word_counts_per_read, all_outposts_summary


def main():
    start_time_total = time.time()

    parser = argparse.ArgumentParser(description="Find outposts in reads from a FASTQ file and analyze irreducible word lengths.")
    parser.add_argument("--input", type=str, help="Path to the FASTQ file.")
    parser.add_argument("--debug", action="store_true", help="Enable debug output.")
    parser.add_argument("--plot", action="store_true", help="Generate a histogram of irreducible word lengths.")
    parser.add_argument("--output_gfa", type=str, help="Path to save the GFA output.")
    parser.add_argument("--output_fasta", type=str, help="Path to save the FASTA output.")
    parser.add_argument("--output_irreducible_words_csv", type=str, help="Path to save all mapped irreducible words as a CSV file.")
    parser.add_argument("--output_common_irreducible_words_csv", type=str, help="Path to save common irreducible words as a CSV file with detailed info.")
    parser.add_argument("--output_right_x_csv", type=str, help="Path to save Right_x outposts and words as a CSV file.")
    parser.add_argument("--threads", type=int, default=1, help="Number of CPUs to use (default: 1).")
    parser.add_argument("--memory", type=str, default="4G", help="Maximum memory to use (e.g., 4G, 8G).")
    args = parser.parse_args()

    print(f"CPUs: {args.threads}")

    # Initialize process for memory tracking in main
    process = psutil.Process(os.getpid())
    def get_memory_usage_mb():
        return process.memory_info().rss / (1024 * 1024)

    print("Start reading file...")
    start_time_read = time.time()
    try:
        reads = []
        for record in BioSeqIO.parse(args.input, "fastq"):
            reads.append(str(record.seq))
    except FileNotFoundError:
        print(f"Error: FASTQ file '{args.input}' not found.")
        return
    except Exception as e:
        print(f"Error parsing FASTQ file: {e}")
        return

    if not reads:
        print("No reads were found in the fastq file")
        return
    end_time_read = time.time()
    print("Reading done.")
    print(f"Number reads: {len(reads)}")

    # --- Use multiprocessing.Manager to create a shared list of reads ---
    # This avoids copying the 'reads' list for each worker's `initargs`
    manager = multiprocessing.Manager()
    shared_reads = manager.list(reads) # Create a shared list from your reads

    print("Start construction Suffix Tree...")
    start_time_st = time.time()
    end_time_st = time.time()
    print("Suffix Tree construction will happen in worker processes.")
    st_construction_time = end_time_st - start_time_st
    print(f"Time for initial setup of shared reads: {st_construction_time:.2f} seconds")


    # Memory info for the main process (before parallel work)
    process = psutil.Process()
    memory_info_main = process.memory_info()
    virtual_memory_main = memory_info_main.vms
    resident_set_size_main = memory_info_main.rss
    print(f"Amount virtual memory occupied by main process (after reading, before parallel): {virtual_memory_main / (1024 ** 2):.2f} MB")
    print(f"Resident set size occupied by main process (after reading, before parallel): {resident_set_size_main / (1024 ** 2):.2f} MB")

        # after reading fasta/fastq into `reads`
    print("\nComputing Right_x, Left_x, and MUS Words (mmap + ProcessPoolExecutor)...")
    start_time_irreducible = time.time()

    # Use a short custom prefix (optional) to help identify temp files
    tmp_prefix = "reads_mmap"

    all_right_x, all_left_x, all_mapped_irreducible_words_with_indices, all_irreducible_word_lengths, \
        irreducible_word_counts_per_read, all_outposts_summary = compute_irreducible_words_parallel_mmap(
            reads, args, tmp_prefix=tmp_prefix
        )

    end_time_irreducible = time.time()
    time_irreducible = end_time_irreducible - start_time_irreducible
    print(f"Time taken to compute Right_x, Left_x, and Mapped MUS Words (mmap parallel): {time_irreducible:.2f} seconds")

    # Outposts (Right Outposts from Right_x and Left Outposts from Left_x)
    
    # --- Print all_right_x and their associated words ---
    print("\n--- All Right Outposts (Right_x) and Words ---")
    for read_idx, rx_list in enumerate(all_right_x):
       current_read = reads[read_idx]
       print(f"Read {read_idx}:")
       for pos_idx, rx_val in enumerate(rx_list):
           if rx_val is not None:
               # Check if rx_val is positive AND less or equal to the length of current_read 
               if rx_val > 0 and rx_val <= len(current_read):
                   rx_word = current_read[pos_idx: rx_val]
               else:
                   rx_word = "" # Assign empty string if conditions are not met
               print(f"Position {pos_idx} (0-based): Length {len(rx_word)}, Word: '{rx_word}")
           else:
               print(f"  Position {pos_idx} (0-based): Error/No data for Right_x")          	
    print("---------------------------------------------")
    
     # --- Save All Right Outposts (Right_x) and Words to CSV ---
    if args.output_right_x_csv:
        print(f"\nSaving all Right Outposts (Right_x) and Words to CSV: {args.output_right_x_csv}")
        try:
            with open(args.output_right_x_csv, 'w', newline='') as csvfile:
                csv_writer = csv.writer(csvfile)
                # Write header row for Right_x data
                csv_writer.writerow(['Read_Index', 'Position_0_Based', 'Right_x_Length', 'Word'])

                for read_idx, rx_list in enumerate(all_right_x):
                    current_read = reads[read_idx]
                    for pos_idx, rx_val in enumerate(rx_list):
                        rx_word_to_write = ""
                        if rx_val is not None:
                            if rx_val > 0 and rx_val <= len(current_read):
                                rx_word_to_write = current_read[pos_idx: rx_val]
                            # If rx_val is not positive or not within bounds, rx_word_to_write remains ""
                        
                        # Only write the row if the word is not an empty string
                        if len(rx_word_to_write) > 0:
                            csv_writer.writerow([read_idx, pos_idx, len(rx_word_to_write), rx_word_to_write])
                        
                        # Write the row. If rx_val is None, write "N/A" for length.
                        #csv_writer.writerow([read_idx, pos_idx, len(rx_word_to_write), rx_word_to_write])
            print("Right Outposts (Right_x) and Words saved successfully.")
        except IOError as e:
            print(f"Error saving Right Outposts (Right_x) to CSV: {e}")
    else:
        print("\nNo output CSV file specified for Right Outposts (Right_x).")
    # --- End Save All Right Outposts (Right_x) and Words to CSV ---


    # --- Print all_left_x and their associated words ---
    print("\n--- All Left Outposts (Left_x) and Words ---")
    for read_idx, lx_list in enumerate(all_left_x):
       current_read = reads[read_idx]
       print(f"Read {read_idx}:")
       for pos_idx, lx_val in enumerate(lx_list):
           if lx_val is not None:
               lx_word = current_read[pos_idx - lx_val + 1: pos_idx + 1] if lx_val > 0 else ""
               print(f"  Position {pos_idx} (0-based): Length {lx_val}, Word: '{lx_word}'")
           else:
               print(f"  Position {pos_idx} (0-based): Error/No data for Left_x")
    print("---------------------------------------------")

    # --- All Outposts Summary (as previously) ---
    print("\n--- All Outposts Summary ---")
    for i, outpost_data in enumerate(all_outposts_summary):
       print(f"Read {i}: {outpost_data}")
    print("---------------------------------------------")

    # Extract words from (word, read_index) tuples for direct display/count
    all_irreducible_words_only = [word_info[0] for sublist in all_mapped_irreducible_words_with_indices for word_info in sublist]

    print("\n--- All Mapped Irreducible Words per Read ---")
    # For printing, we still want to show words per read without the extra index info
    for i, words_with_indices_in_read in enumerate(all_mapped_irreducible_words_with_indices):
       words_only_in_read = [word_info[0] for word_info in words_with_indices_in_read]
       print(f"Read {i}: {words_only_in_read}")
    print("---------------------------------------------")

    # --- Save All Irreducible Words to CSV ---
    if args.output_irreducible_words_csv:
        print(f"\nSaving all irreducible words to CSV: {args.output_irreducible_words_csv}")
        try:
            with open(args.output_irreducible_words_csv, 'w', newline='') as csvfile:
                csv_writer = csv.writer(csvfile)
                csv_writer.writerow(['Read_Index', 'Irreducible_Word_Index', 'Irreducible_Word', 'Length'])
                for i, words_with_indices_in_read in enumerate(all_mapped_irreducible_words_with_indices):
                    for j, word_info in enumerate(words_with_indices_in_read):
                        word = word_info[0] # Extract the word from the tuple
                        word_length = len(word) if isinstance(word, str) and word not in ["Invalid Indices", "Error_in_Read_Processing"] else 0
                        csv_writer.writerow([i, j, word, word_length])
            print("Irreducible words saved successfully.")
        except IOError as e:
            print(f"Error saving irreducible words to CSV: {e}")
    else:
        print("\nNo output CSV file specified for all irreducible words.")
    # --- End Save All Irreducible Words to CSV ---

    # --- Calculate common irreducible words with occurrences ---
    # Use defaultdict to store counts and a set of read indices for each word
    common_word_details = defaultdict(lambda: {"count": 0, "reads": set()})

    for sublist in all_mapped_irreducible_words_with_indices:
        for word, read_idx in sublist:
            if word not in ["Invalid Indices", "Error_in_Read_Processing"]:
                common_word_details[word]["count"] += 1
                common_word_details[word]["reads"].add(read_idx)

    # Filter for words that appear more than once
    final_common_irreducible_words = {
        word: details for word, details in common_word_details.items() if details["count"] > 1
    }

    print("\n--- Common Irreducible Words Across Reads (Count > 1) ---")
    if final_common_irreducible_words:
       # Sort by count descending for better readability in console
       sorted_common_words_for_display = sorted(final_common_irreducible_words.items(), key=lambda item: item[1]["count"], reverse=True)
       for word, details in sorted_common_words_for_display:
           # Convert set of read indices to a sorted list for consistent display
           read_indices_str = ", ".join(map(str, sorted(list(details["reads"]))))
           print(f"Word: '{word}', Count: {details['count']}, Length: {len(word)}, Occurs in Reads: [{read_indices_str}]")
    else:
       print("No irreducible words are common across more than one read.")

    # --- Save Common Irreducible Words to CSV ---
    if args.output_common_irreducible_words_csv:
        print(f"\nSaving common irreducible words to CSV: {args.output_common_irreducible_words_csv}")
        try:
            with open(args.output_common_irreducible_words_csv, 'w', newline='') as csvfile:
                csv_writer = csv.writer(csvfile)
                # Updated header for common words CSV
                csv_writer.writerow(['Word', 'Length', 'Count', 'Reads_Occurred_In'])

                # Sort by count descending for CSV output
                sorted_common_words_for_csv = sorted(final_common_irreducible_words.items(), key=lambda item: item[1]["count"], reverse=True)
                for word, details in sorted_common_words_for_csv:
                    word_length = len(word)
                    count = details["count"]
                    # Convert the set of read indices to a comma-separated string
                    reads_occurred_in = ", ".join(map(str, sorted(list(details["reads"]))))
                    csv_writer.writerow([word, word_length, count, reads_occurred_in])
            print("Common irreducible words saved successfully.")
        except IOError as e:
            print(f"Error saving common irreducible words to CSV: {e}")
    else:
        print("\nNo output CSV file specified for common irreducible words.")
    # --- End Save Common Irreducible Words to CSV ---


    # Ensure build_assembly_graph_for_eulerian receives a list of lists of words (not tuples)
    # This assumes build_assembly_graph_for_eulerian expects [[word1, word2], [word3], ...]
    # We need to flatten the (word, read_index) tuples back to just words.
    print("\nBuilding the assembly graph...")
    start_time_assembly_graph = time.time()
    
     # --- Get memory usage BEFORE building the graph ---
    process = psutil.Process(os.getpid())
    mem_before_graph = process.memory_info().rss
    
    all_mapped_irreducible_words_for_graph = [[word_info[0] for word_info in sublist if word_info[0] not in ["Invalid Indices", "Error_in_Read_Processing"]] for sublist in all_mapped_irreducible_words_with_indices]
    assembly_graph_euler, edge_data, node_prefix_suffix = build_assembly_graph_for_eulerian(all_mapped_irreducible_words_for_graph, reads)
    
    end_time_assembly_graph = time.time()
    
    # --- Get memory usage AFTER building the graph ---
    mem_after_graph = process.memory_info().rss
    # --- Calculate and print the time and memory delta ---
    time_assembly_graph = end_time_assembly_graph - start_time_assembly_graph
    mem_graph_delta = mem_after_graph - mem_before_graph
    
    print(f"Time taken to build the assembly graph: {time_assembly_graph:.2f} seconds")
    print(f"Memory (RSS delta) used to build the assembly graph: {mem_graph_delta / (1024 ** 2):.2f} MB")

    print("\n--- Assembly Graph (Eulerian) ---")
    for node, edges in assembly_graph_euler.items():
        print(f"Node: '{node}' -> {edges}")

    print("\n--- Edge Data (Eulerian) ---")
    for edge, data_list in edge_data.items():
        print(f"Edge: {edge} -> {data_list}")

    
    
    print("\n--- Running Contig Extraction Strategies ---")

    # --- STEP 0.5 + STEP 1: Linear-time Unitig Extraction (collapse + extract) ---
    print("\n--- Step 0.5+1: Extracting Maximal Non-Branching Paths (Unitigs) ---")
    start_time_u = time.time()
    start_mem_u = get_memory_usage_mb()

    # original graph and edge_data are used as base_graph/base_edge_data
    # assembly_graph_euler and edge_data come from your earlier build_assembly_graph_for_eulerian call
    unitigs = extract_unitigs_linear(assembly_graph_euler, edge_data)

    end_time_u = time.time()
    end_mem_u = get_memory_usage_mb()

    print(f"Found {len(unitigs)} unitigs.")
    for u in unitigs[:10]:   # print first 10 for brevity
        seq_snip = u.sequence if hasattr(u, "sequence") else ""
        print(f"  {u.id}: {u.path} (seq: '{seq_snip[:40]}')")

    print(f"Step 0.5+1 Time: {end_time_u - start_time_u:.2f}s | Memory Delta: {end_mem_u - start_mem_u:.2f} MB")
    print("---------------------------------")

    # --- STEP 2: Build the Unitig Graph (use the original base graph/edge_data) ---
    print("\n--- Step 2: Building Unitig Graph ---")
    start_time_step2 = time.time()
    start_mem_step2 = get_memory_usage_mb()

    builder = UnitigGraphBuilder(assembly_graph_euler, edge_data, unitigs)
    unitig_graph, unitig_edge_data = builder.build_graph()

    end_time_step2 = time.time()
    end_mem_step2 = get_memory_usage_mb()

    print("Unitig Graph Edges:")
    for u_id, neighbors in unitig_graph.items():
        for v_id, flow in neighbors.items():
            print(f"  {u_id} -> {v_id} (flow: {flow})")

    print(f"Step 2 Time: {end_time_step2 - start_time_step2:.2f}s | Memory Delta: {end_mem_step2 - start_mem_step2:.2f} MB")
    print("---------------------------------")

    # --- STEP 3: Resolve the Unitig Graph (unchanged) ---
    print("\n--- Step 3: Resolving Unitig Graph ---")
    start_time_step3 = time.time()
    start_mem_step3 = get_memory_usage_mb()

    resolver = UnitigGraphResolver(unitig_graph, unitig_edge_data, unitigs, edge_data)
    resolved_contigs = resolver.resolve_and_assemble()   # this returns list[str]

    end_time_step3 = time.time()
    end_mem_step3 = get_memory_usage_mb()

    print(f"Resolved {len(resolved_contigs)} contigs:")
    for i, c in enumerate(resolved_contigs[:10]):
        print(f"  Resolved Contig {i+1}: {c[:60]}")

    print(f"Step 3 Time: {end_time_step3 - start_time_step3:.2f}s | Memory Delta: {end_mem_step3 - start_mem_step3:.2f} MB")
    print("----------------------------------")

    # --- STEP 4: Overlap Merging (Scaffolding) using optimized prefix-KMP ---
    print("\n--- Step 4: Overlap Merging (Scaffolding) ---")
    start_time_step4 = time.time()
    start_mem_step4 = get_memory_usage_mb()

    # Pre-filter to reduce candidate count
    pre_filtered = filter_non_maximal_contigs(resolved_contigs)

    # Merge — tuned parameters:
    MIN_OVERLAP_THRESHOLD = 2     # use 3..k-1; increase for more conservative merging
    SEED_LEN = 3                 # good default for moderate contig lengths

    merged_contigs = prefix_kmp_scaffolder(pre_filtered, min_overlap=MIN_OVERLAP_THRESHOLD, seed_len=SEED_LEN)

    end_time_step4 = time.time()
    end_mem_step4 = get_memory_usage_mb()

    print(f"Merged into {len(merged_contigs)} contigs:")
    for i, c in enumerate(merged_contigs[:10]):
        print(f"  Merged Contig {i+1}: {c[:60]}")

    print(f"Step 4 Time: {end_time_step4 - start_time_step4:.2f}s | Memory Delta: {end_mem_step4 - start_mem_step4:.2f} MB")
    print("----------------------------------")

    # Final printout
    print(f"\n--- Final Assembled Maximal Contigs ({len(merged_contigs)}):")
    for i, c in enumerate(merged_contigs):
        print(f"  {i+1}: {c[:120]}")



    if args.output_fasta:
    
        output_fasta = args.output_fasta
        if not output_fasta.lower().endswith(".fasta"):
            output_fasta += ".fasta"

        fasta_records = []

        for i, contig in enumerate(merged_contigs):
            record = SeqRecord.SeqRecord(
                Seq(contig),
                id=f"contig_{i+1}",
                description=f"length={len(contig)}"
            )
            fasta_records.append(record)

    # Write to FASTA
        with open(output_fasta, "w") as f:
            BioSeqIO.write(fasta_records, f, "fasta")

        print(f"\nContigs saved to FASTA: {output_fasta}")

    else:
        print("\nNo output FASTA file specified.")

    end_time_total = time.time()
    total_time = end_time_total - start_time_total
    print(f"\nTotal time taken: {total_time:.2f} seconds")

if __name__ == "__main__":

    multiprocessing.freeze_support()
    main()

    
    
    

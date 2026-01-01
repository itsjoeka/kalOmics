# Genome Assembly Pipeline

A high-performance de novo genome assembly pipeline using irreducible word decomposition and graph-based algorithms.

## Features

- **Memory-efficient parallel processing** using mmap-backed storage
- **Suffix tree-based outpost detection** for identifying assembly units
- **Irreducible word computation** (MUS - Maximal Unique Segments)
- **Multi-strategy assembly**:
  - Linear unitig extraction
  - Eulerian graph-based assembly
  - KMP-optimized scaffolding
- **Graph simplification** and resolution algorithms
- **Comprehensive output formats**: FASTA, GFA, CSV reports

## Installation

### Prerequisites

- Python 3.8 or higher
- 4GB RAM minimum (8GB+ recommended for larger datasets)

### From Source

```bash
# Clone the repository
git clone https://github.com/yourusername/genome-assembler.git
cd genome-assembler

# Install dependencies
pip install -r requirements.txt

# Install the package
pip install -e .
```

### Using pip (after publishing)

```bash
pip install genome-assembler
```

## Quick Start

### Basic Assembly

```bash
genome-assembler --input reads.fastq --output_fasta assembly.fasta --threads 4
```

### Complete Analysis with Reports

```bash
genome-assembler \
  --input reads.fastq \
  --output_fasta assembly.fasta \
  --output_gfa assembly.gfa \
  --output_irreducible_words_csv mus_words.csv \
  --output_common_irreducible_words_csv common_mus.csv \
  --output_right_x_csv right_outposts.csv \
  --threads 8 \
  --memory 16G \
  --debug
```

## Usage

### Command-line Options

```
Required Arguments:
  --input PATH              Input FASTQ file containing reads

Output Options:
  --output_fasta PATH       Output assembled contigs in FASTA format
  --output_gfa PATH         Output assembly graph in GFA format
  --output_irreducible_words_csv PATH
                           All irreducible words (MUS) per read
  --output_common_irreducible_words_csv PATH
                           Common irreducible words across reads
  --output_right_x_csv PATH Right outpost positions and sequences

Performance Options:
  --threads INT            Number of CPU threads (default: 1)
  --memory SIZE            Maximum memory usage (e.g., 4G, 8G)

Analysis Options:
  --debug                  Enable detailed debug output
  --plot                   Generate histogram of word lengths
```

### Example Workflows

#### 1. Quick Assembly (Single Thread)

```bash
genome-assembler --input small_dataset.fastq --output_fasta contigs.fasta
```

#### 2. High-Performance Assembly (Multi-threaded)

```bash
genome-assembler \
  --input large_dataset.fastq \
  --output_fasta assembly.fasta \
  --threads 16 \
  --memory 32G
```

#### 3. Research Mode (Full Diagnostics)

```bash
genome-assembler \
  --input reads.fastq \
  --output_fasta assembly.fasta \
  --output_irreducible_words_csv all_mus.csv \
  --output_common_irreducible_words_csv common_mus.csv \
  --output_right_x_csv outposts.csv \
  --debug \
  --plot
```

## Algorithm Overview

This assembler implements a novel approach based on irreducible word decomposition:

1. **Outpost Detection**: Identifies unique sequence positions using suffix trees
2. **Irreducible Word Extraction**: Computes maximal unique substrings (MUS)
3. **Graph Construction**: Builds de Bruijn-style graphs from MUS words
4. **Unitig Extraction**: Identifies maximal non-branching paths
5. **Graph Resolution**: Resolves complex branching structures
6. **Scaffolding**: Merges contigs using prefix-KMP algorithm

For detailed algorithm description, see [docs/algorithm.md](docs/algorithm.md).

## Output Files

### FASTA Output
Standard multi-FASTA format with assembled contigs:
```
>contig_1 length=1234
ATCGATCGATCG...
>contig_2 length=567
GCTAGCTAGCTA...
```

### CSV Reports

**Irreducible Words CSV**: All MUS words identified per read
**Common Words CSV**: MUS words appearing in multiple reads
**Right Outposts CSV**: Right outpost positions and sequences

### GFA Output (if specified)
Assembly graph in GFA format for visualization with Bandage.

## Performance Considerations

- **Memory Usage**: Approximately 4-8× the size of input FASTQ
- **CPU Utilization**: Scales linearly up to ~16 threads
- **Disk I/O**: Temporary files created in system temp directory
- **Recommended**: SSD storage for large datasets (>1GB)

## Testing

Run the test suite:

```bash
# Install test dependencies
pip install pytest pytest-cov

# Run all tests
pytest tests/

# Run with coverage
pytest --cov=genome_assembler tests/
```

## Citation

If you use this software in your research, please cite:

```bibtex
@phdthesis{yourname2024,
  title={Genome Assembly Using Irreducible Word Decomposition},
  author={Your Name},
  year={2024},
  school={Your University}
}
```

## Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

See [CONTRIBUTING.md](CONTRIBUTING.md) for detailed guidelines.

## License

This project is licensed under the MIT License - see [LICENSE](LICENSE) file for details.

## Acknowledgments

- BioPython community for sequence I/O utilities
- Contributors to suffix tree algorithms
- Your PhD advisor and committee members

## Support

- **Issues**: [GitHub Issues](https://github.com/yourusername/genome-assembler/issues)
- **Discussions**: [GitHub Discussions](https://github.com/yourusername/genome-assembler/discussions)
- **Email**: your.email@university.edu

## Roadmap

- [ ] Support for paired-end reads
- [ ] Integration with quality score filtering
- [ ] GPU acceleration for suffix tree construction
- [ ] Long-read assembly support (PacBio, Nanopore)
- [ ] Web interface for small datasets

# References

# contributors

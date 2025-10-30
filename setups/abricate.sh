conda config --add channels bioconda
conda create -n abricate -y -c bioconda abricate
conda activate abricate
abricate -v
abricate-get_db --db vfdb
abricate-get_db --db card
abricate-get_db --db resfinder
abricate-get_db --db ncbi
abricate --setupdb
abricate --csv --db vfdb *.fasta > vfdb_abricate.csv

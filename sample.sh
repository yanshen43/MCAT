cd "${0%/*}"

# runs with all sequences
./orange_pipeline.py data/positive_genes.fasta data/n_genes.fasta

# runs with only S. pombe sequences
./orange_pipeline.py data/pombe+.fasta data/pombe-.fasta

# runs with only sequences associated with the mad2 gene
./orange_pipeline.py data/mad2.txt data/negative_mad2.txt
#!/bin/bash
#put files and this script in the same directory as dynamit
for filename in *.fasta; do
	python -m dynamit generated_dna_config.txt "$filename" -o dynamit_results
done

#!/bin/bash

for filename in /home/ccogs/gasp/test_data/syn/*
do
SEQIN="$filename"

INPUT_TOPDIR=/home/ccogs/gasp/test_data/syn
OUTPUT_TOPDIR=/home/ccogs/gasp/test_results/xmotif_synthetic_results
BASENAME=${SEQIN%.*}
OUTDIR=${OUTPUT_TOPDIR}/$BASENAME

rm -rf "$OUTDIR"
mkdir -p "$OUTDIR"


	# starts XXmotif with the following options: 
	#    --zoops : zero or one occurrences per sequence (DEFAULT)
	#    --format FASTA: input sequences are in fasta format (DEFAULT)
	#    --type ALL: use FIVEMERS, PALINDROMIC and TANDEMIC seeds (DEFAULT)
	#    --localization : incorporate P-values for positioned motifs
	#    --downstream: defines the number of nucleotides downstream of the anchor point (sequences are cut from -300 to +100 w.r.t TSS)
	#    --merge-motif-threshold LOW: leads to a nonredundant small results list, but has the risk to loose some motifs
	#    --XXmasker: mask homology, repeats and low complexity regions in input sequences
	./XXmotif $OUTDIR "$SEQIN" \
		--zoops --format FASTA --type ALL --localization \
		--no-graphics --downstream 100 --merge-motif-threshold MEDIUM --XXmasker \


done

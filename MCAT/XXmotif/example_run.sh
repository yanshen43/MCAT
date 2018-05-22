#!/bin/sh

TYPE=$1

case "$TYPE" in
	1)
		SEQIN=dmel_segmentation_modules.mfasta
		;;
	2)
		SEQIN=hsap_core_promoters_ribosomal.fasta
		;;
	3)
		SEQIN=hsap_core_promoters_all.fasta
		;;
	*)
		echo "Usage: $0 NUMBER, specify either
	(1) D.mel segmentation modules (see paper Figure 5),
	(2) H.sap core promoters (ribosomal, see paper Figure 6B), or
	(3) H.sap core promoters (all, see paper Figure 6A)."
		exit 1
esac

INPUT_TOPDIR=example_input
OUTPUT_TOPDIR=example_output
BASENAME=${SEQIN%.*}
OUTDIR=${OUTPUT_TOPDIR}/$BASENAME

rm -rf "$OUTDIR"
mkdir -p "$OUTDIR"

if [ "$TYPE" -eq 1 ]; then

	# starts XXmotif with the following options:
	#    --mops : multiple occurrences per sequence
	#    --revcomp: reverse compliment
	#    --format MFASTA: input sequences are in multiple fasta format
	#    --maxMultiplesequences 13: use only the alignment of 13 species (12 Drosophila species + 1 outgroup)
	#    --type ALL: use FIVEMERS, PALINDROMIC and TANDEMIC seeds (DEFAULT)
	#    --gaps 1: scan all fivemers including at most 1 gap position in the seeds phase
	#    --merge-motif-threshold LOW: leads to a nonredundant small results list, but has the risk to loose some motifs
	#    --XXmasker: mask homology, repeats and low complexity regions in input sequences
	./XXmotif $OUTDIR "$INPUT_TOPDIR/$SEQIN" \
		--mops --revcomp --format MFASTA --maxMultipleSequences 13 \
		--type ALL --gaps 1 --merge-motif-threshold MEDIUM --XXmasker \
		2>&1 | tee $OUTDIR/$BASENAME.xxlog

elif [ "$TYPE" -eq 2 ]; then

	# starts XXmotif with the following options: 
	#    --zoops : zero or one occurrences per sequence (DEFAULT)
	#    --format FASTA: input sequences are in fasta format (DEFAULT)
	#    --type ALL: use FIVEMERS, PALINDROMIC and TANDEMIC seeds (DEFAULT)
	#    --localization : incorporate P-values for positioned motifs
	#    --downstream: defines the number of nucleotides downstream of the anchor point (sequences are cut from -300 to +100 w.r.t TSS)
	#    --merge-motif-threshold LOW: leads to a nonredundant small results list, but has the risk to loose some motifs
	#    --XXmasker: mask homology, repeats and low complexity regions in input sequences
	./XXmotif $OUTDIR "$INPUT_TOPDIR/$SEQIN" \
		--zoops --format FASTA --type ALL --localization \
		--downstream 100 --merge-motif-threshold MEDIUM --XXmasker \
		2>&1 | tee $OUTDIR/$BASENAME.xxlog

elif [ "$TYPE" -eq 3 ]; then

	# starts XXmotif with the following options: 
	#    --zoops : zero or one occurrences per sequence (DEFAULT)
	#    --format FASTA: input sequences are in fasta format (DEFAULT)
	#    --type ALL: use FIVEMERS, PALINDROMIC and TANDEMIC seeds (DEFAULT)
	#    --localization : incorporate P-values for positioned motifs
	#    --downstream: defines the number of nucleotides downstream of the anchor point (sequences are cut from -300 to +100 w.r.t TSS)
	#    --merge-motif-threshold LOW: leads to a nonredundant small results list, but has the risk to loose some motifs
	#    --XXmasker: mask homology, repeats and low complexity regions in input sequences
	./XXmotif $OUTDIR "$INPUT_TOPDIR/$SEQIN" \
		--zoops --format FASTA --type ALL --localization \
		--downstream 100 --merge-motif-threshold MEDIUM --XXmasker \
		2>&1 | tee $OUTDIR/$BASENAME.xxlog

fi

#!/bin/bash
#put files and this script in the same directory as emd
#make sure this is running on baobab
for filename in *.fasta
do
	#cleanup
	rm -r step5regulon
	rm -r RMDS*.o* runmotif-MD*.sh RMDS*.log step6*

	#set up environment
	mkdir step5regulon
	cp "$filename" step5regulon/

	#run tests
	emdrunPX.pl -f emd.cfg -w 15 -n 5
	sleep 20s
	emdMotif.pl -f step5regulon/"$filename" -c emd.cfg -n 5 >> emd_results/"$filename"_results.out
	
done

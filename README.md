This is a pipeline for finding motifs in fasta files.
It can be run from the command line as follows:

usage: orange_pipeline_refine.py [-h] [-w W] [--nmotifs NMOTIFS] [--iter ITER] [-c C]
                                 [-s S] [-d] [-ff] [-v V]
                                 positive_seq negative_seq

positional arguments:
  positive_seq       the fasta file for the positive sequences
  negative_seq       the fasta file for the negative sequences

optional arguments:
  -h, --help         show this help message and exit
  -w W               motif width
  --nmotifs NMOTIFS  max number of distinct motifs to look for
  --iter ITER        max number of iterations (for DECOD tool)
  -c C               config filename
  -s S               statfile directory
  -d                 turn debug mode on
  -ff                force (re)filtering
  -v V               number of results to include in visual output


This pipeline makes use of 6 publicly available 3rd party tools:
DECOD
MEME
CMF
BioProspector
Weeder
XXmotif

This distribution includes the source for all of these tools,
as well as 64-bit and 32-bit executables compiled on linux.
If you are using Windows or MacOS, you might need to recompile
some of the tools on your platform and replace the relevant
executables.
All of these sources were downloaded from their respective websites
and are unaltered with the exception of 2 small file path changes in
Weeder to account for our pipeline not providing files the way weeded
expected.

Requirements:
Java 1.8
Python 2.7
    SciPy 0.19.0
    NumPy 1.12.1
    MatPlotLib 1.2.0
perl 5.16.3
Chostscript 9.07

![alt text](https://github.com/yanshen43/MCAT/blob/master/vis.png)


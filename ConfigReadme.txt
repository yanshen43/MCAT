orange_pipeline can now read options from a config file!

On the command line, you can use the -c option to specify the config file name

The format for the config file is as follows:

To select which tool you want to edit the command line options for type 'tool:'
followed immediately by the tool name in all uppercase on a new line. For 
example, 'tool:MEME' would specify that the following options should be passed
to MEME.

To add a command line option to the currently selected tool type '-opt=val' on
a new line where opt is the option you want to set and val is the value that you
want to set it to. For example, '-oc=results/meme/' would add '-oc results/meme/'
to the command line for the currently selected tool. If there is no value
associated with a given option, type '-opt=', ending the line after the '='.

To delete a command line option from the currently selected tool type '!-opt' on
a new line where opt is the option that you want to delete. For example, '!-oc'
would remove the oc option from the command line of the currently selected tool.


Any options that you specify in a config file will overwrite the default for
that option for that tool. For example, the default value for -maxm for weeder
is 10. If you put '-maxm=7' under 'tool:WEEDER' in the config file, then the
pipeline will have '-maxm 7' on weeder's command line instead of '-maxm 10'.


Here is an example of a config file which sets meme's output directory to
'memeOut/', prevents the -dna option from being passed to meme, adds the -rna
option to meme, sets the max motif length for cmf to 13, and sets
BioProspecotor's number of motifs to report to 7.

tool:MEME
-oc=memeOut
!-dna
-rna=
tool:CMF
-u=13
tool:BIOPROSPECTOR
-r=7



Here are all the command line options as reported by each tool.


Command line options for Weeder:

-f input_file 
input_file must be in Multi-FASTA format.

-O frequency_file_organism_code

When no organism code for oligo frequency files is provided it is assumed to be HS (Homo sapiens).

Oligo frequency files for the following organisms are available in the standard Weeder 2.0 package: 

Homo sapiens - Code: HS
Mus musculus - Code: MM
Drosophila melanogaster - Code: DM
Saccharomyces cerevisiae - Code: SC
Arabidopsis thaliana - Code: AT

Other frequency files  may be added to the FreqFiles directory by using the "Frequency maker" program
available at http://www.beaconlab.it/modtools

-chipseq
This flag activates the ChIP-Seq heuristic speeding-up the computation.

-top <num> (DEFAULT: 100)
If the -chipseq parameter is used Weeder 2.0 scans all the input sequences for occurrences of the oligos contained in the top <num> input sequences.
Increase this value when your input has many more than <num> sequences to improve the chance of finding motifs enriched only 
in a subset of your input sequences.

-b <num> (DEFAULT: 50)
Weeder 2.0 builds occurrences matrix profiles and outputs (if other conditions are met) only the top <num> scoring motifs
for each motif length. Increase this value to have more (lowest scoring) motifs in the output (see also -maxm).

-maxm <num> (DEFAULT: 25)
To limit the output length, Weeder 2.0 reports only the top <num> scoring motifs with their associated occurrences
matrix and occurrences list. Increase <num> to have longer outputs with more lowest scoring motifs.

-ss
Single strand mode.

ADVANCED OPTIONS

-sim <num> (DEFAULT: 0.95 MIN: 0 MAX: 1)
Similarity threshold for the redundancy filter. This filter removes from the output those motifs that are too similar to other motifs
already reported. Values close to 0 mean a stricter filter and vice versa values close to 1 impose a looser filter.
Set <num> to 1 to disable the filter altogether. Set it to 0 to have in the output only the top scoring oligo for each one of
the possible oligo lengths (6, 8 and 10).

-em <num> (DEFAULT: 1 MIN: 0 MAX: 100)
Weeder 2.0 has a built-in expectation maximization (EM) matrix profiles refinement step.
<num> defines the number of EM cycles to be performed by Weeder 2.0.
One (default) or few EM cycles should be sufficient to "clean" matrix profiles without overfitting them.


Command line options for Meme:
(note that there isn't an option for specifying the input file, so this necessarily comes from
the pipline's positive_seq parameter)

[-h]                    print this message
[-objfun classic|nc|smhg|cv|nz|ll]      objective function (default: classic)
[-o <output dir>]       name of directory for output files
                        will not replace existing directory
[-oc <output dir>]      name of directory for output files
                        will replace existing directory
[-text]                 output in text format (default is HTML)
[-dna]                  sequences use DNA alphabet
[-rna]                  sequences use RNA alphabet
[-protein]              sequences use protein alphabet
[-alph <alph file>]     sequences use custom alphabet
[-mod oops|zoops|anr]   distribution of motifs
[-nmotifs <nmotifs>]    maximum number of motifs to find
[-evt <ev>]             stop if motif E-value greater than <evt>
[-nsites <sites>]       number of sites for each motif
[-minsites <minsites>]  minimum number of sites for each motif
[-maxsites <maxsites>]  maximum number of sites for each motif
[-wnsites <wnsites>]    weight on expected number of sites
[-w <w>]                motif width
[-minw <minw>]          minimum motif width
[-maxw <maxw>]          maximum motif width
[-nomatrim]             do not adjust motif width using multiple
                        alignment
[-wg <wg>]              gap opening cost for multiple alignments
[-ws <ws>]              gap extension cost for multiple alignments
[-noendgaps]            do not count end gaps in multiple alignments
[-bfile <bfile>]        name of background Markov model file
[-revcomp]              allow sites on + or - DNA strands
[-pal]                  force palindromes (requires -dna)
[-maxiter <maxiter>]    maximum EM iterations to run
[-distance <distance>]  EM convergence criterion
[-psp <pspfile>]        name of positional priors file
[-prior dirichlet|dmix|mega|megap|addone]
                        type of prior to use
[-b <b>]                strength of the prior
[-plib <plib>]          name of Dirichlet prior file
[-maxwords <maxwords>]  maximum number of words to test as EM starts
[-spfuzz <spfuzz>]      fuzziness of sequence to theta mapping
[-spmap uni|pam]        starting point seq to theta mapping type
[-cons <cons>]          consensus sequence to start EM from
[-heapsize <hs>]        size of heaps for widths where substring 
                        search occurs
[-x_branch]             perform x-branching
[-w_branch]             perform width branching
[-allw]                 include all motif widths from min to max
[-bfactor <bf>]         branching factor for branching search
[-maxsize <maxsize>]    maximum dataset size in characters
[-nostatus]             do not print progress reports to terminal
[-p <np>]               use parallel version with <np> processors
[-time <t>]             quit before <t> CPU seconds consumed
[-sf <sf>]              print <sf> as name of sequence file
[-V]                    verbose mode
[-version]              display the version number and exit
[-seed <seed>]          seed for random numbers for shuffling and 
                        sampling
[-shuf <kmer>]          preserve frequencies of k-mers of size <kmer> 
                        when shuffling
[-ctfrac <ctfrac>]      fraction of control sequences to use


Command line options for CMF:

-w,  the length of the motif seed
-m,  numer of mismatches in seed
-F,  fdr cutoff
-t,  the number of top seeds to test
-i1, first set of sequences
-i2, second set of sequences
-d,  1: enrichment only in i1, 2: enrichment in both datasets
-l,  lower bound on length of motifs
-u,  upper bound on length of motifs
-o,  output of seed statistics
-f,  folder for all other output (must exist before running)
-c,  if > 0 filter out seeds based on cg content (ex -c 4, filter out seeds with more than 4 C or G's)

Command line options for DECOD:

-nogui
-pos <positive_sequence_file>
-neg <negative_sequence_file>
-w <width>
-nmotif <#motifs to find>
-niter <#iterations>
-o <output_file>
[-exact]
[-withrepeat]
[-strand forward|both]
[-pmotif <#Motif per positive sequence>]
[-countupperlimit <Upper limit of kmer count in a sequence>]

Command line options for BioProspector:

-i    seqfile
-W    <motif width (default 10)>
-o    <output file (default stdout)>
-f    <background distribution file (default seqfile)>
-b    <background sequence file (default input sequences)>
-n    <number of times trying to find motif (default 40)>
-r    <number of top motifs to report (default 5)>
-w    <second motif block width for two-block motif (default 0)>
-p 1  [if two-block motif is palindrome (default 0)]
-G    <max gap between two motif blocks (default 0)>
-g    <min gap between two motif blocks (default 0)>
-d 1  [if only need to examine forward (default 2)]
-a 1  [if every sequence contains the motif (default 0)]
-h 1  [if want more degenerate sites (default fewer sites)]
-e    <expected bases per motif site in the sequences 
          (will use Bayes motif scoring, don't specify if unknown)>

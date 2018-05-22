DaECOD: Fast and accurate discriminative DNA motif finding
=========================================================
Version 1.01 2011.10.24

* You need to have Java JRE installed.

* To use GUI version, simply double click the jar file (DECOD-20110613.jar) 
  after unzip or call:

	java -jar DECOD-20111024.jar

* To use command line version, call

	java -jar DECOD-20110613.jar -nogui -pos <positive filename> -neg 
    <negative filename> [-w <motif width, default 8>] [-nmotif <#motifs 
    to find, default 10>] [-niter <#iterations, default 50>] [-c <motif 
	carinality, default 20>] -o <output filename> [-strand <forward|both>]
	[-exact] [-withrepeat] [-pmotif <#motif occurrence per sequence, 
    default 1.0>] [-countupperlimit <Upper limit of kmer count in
    a sequence, default 2>] [-help]

	Note:

	1. "-nogui" must immediately follow DECOD-20110613.jar (i.e. be the 
       first parameter) in order to use the command-line version.

	2. Input sequences must be in fasta format. Can have "N"s to indicate 
       masked or unknown bases

	3. "-c" can be used to specify motif cardinality if desired, default 20. 

	3. If "-exact" is added, the program will do the exact calculation 
       (recommended only for w<=8). Otherwise the speedup calculations will
       be used by default

	4. If you do NOT want to further avoid reporting motifs with simple 
       compositions (one- or di- nucleotide repeats), you can add the 
       "-withrepeat" option. Otherwise by default the program will ignore 
       such simple kmers.

	5. Use "-strand" to specify whether to search for motif only on the 
	   forward strand ("forward") or both strands ("both") of the input
	   sequences. By default the option is set to "both"

	6. User "-pmotif" to set the assumed number of times that the motif is 
       present in each positive sequence. By default this is set to 1 (
       once per positive sequence). This can only be changed in the 
       command-line version of the program.

	7. Use "-countupperlimit" to set an upper limit to the count of each
       kmer in a sequence. The purpose of this is to attenuate the effect
       of simple repetitive kmers. By default this is set to 2. This can 
       only be changed in the command-line version of the program.

	8. If "-help" is added, a simple summary will be printed.

* This program makes use of a few methods from the Apache Commons Mathematics 
  Library which is published under Apache License Version 2.0. The license 
  files are provided together under the ./apache-commons-math-2.1-license-info
  directory.

* Please cite the following publication if you used this software in your 
  work:

  Peter Huggins, Shan Zhong, Idit Shiff, et al., DECOD: Fast and Accurate 
  Discriminative Motif Finding.

* Please contact Shan Zhong at szhong@andrew.cmu.edu for bug report.


Shan Zhong
Jun 14, 2011


============
Version log:

* 1.01: Released 2011.10.24
   - Fixed wrong formatting of decimal points in international regional settings
   - Allowed parameters to be set when calling the program in command line but asking to run in GUI mode

* V1.0: 2011.06.13 Initial release.


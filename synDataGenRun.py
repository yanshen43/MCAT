# Written by Jeff Robertson (thejar@vt.edu)
# Uses seq_gen.py written by Jake Martinez (jrm98@vt.edu)
# and orange_pipeline.py

from seq_gen import generate
from subprocess import call
import time
import random
import shutil
import os
import sys
import pdb

def main():
    timestamp = time.strftime("%y-%m-%d_%H-%M-%S")
    random.seed(time.time())
    
    resdir = "datasets_" + timestamp
    os.mkdir(resdir)
    with open("dataGenRun_" + timestamp + '.log', 'w') as logFile:
        # print and log s
        def out(s):
            print(str(s))
            print >>logFile, time.asctime(), str(s)
            logFile.flush()
        try:
            # number of fasta files/datasets to generate and run
            dataSetCount = 50
            
            # motif width range (inclusive)
            wRange = (7, 13)
            
            # number of sequences
            N = 35
            
            # sequence length
            n = 1000

            # maximum starting location of motif in sequence
            maxloc = -1
            # minimum starting location of motif in sequence
            minloc = 0
            
            # output directory, used when calling pipeline and for saving results
            directory = "./"
            
            # max number of motifs in a single sequence
            nmotifs = 1
            
            # test with dyad motifs?
            dyad = False

            out("Starting synDataGenRun with:")
            out(timestamp)
            out(dataSetCount)
            out(wRange)
            out(N)
            out(n)
            out(maxloc)
            out(minloc)
            out(directory)
            out(nmotifs)
            
            for dataSetNum in xrange(dataSetCount):
                # output filename, used when calling pipeline and for saving results
                output = time.strftime("%y-%m-%d_%H-%M-%S.fasta")

                # motif length
                w = random.randint(wRange[0], wRange[1])
                # max SNPs/motif
                e = random.randint(1, w/4)
                # prob of a SNP in the motif
                ep = random.randint(25, 80) / 100.0
                # prob motif is in sequence
                P = random.randint(65, 100) / 100.0
                out("Generating data set #" + str(dataSetNum) + " with:")
                out(output)
                out(w)
                out(e)
                out(ep)
                out(P)

                # gen positive fasta
                generate(w, N, n, e, ep, P, minloc, maxloc, output, directory, False, nmotifs, dyad)
                # gen negative fasta
                generate(w, N, n, e, ep, P, minloc, maxloc, "neg_"+output, directory, True, nmotifs, dyad)
                
                #get motif, is there a better way?
                with open(directory+output, 'r') as fin:
                    motif = fin.readline().split("'")[1]
                out("datasets generated with motif " + motif + ", running pipeline:")
                cmdLine = ('./orange_pipeline.py -c synData.cfg -w ' + str(len(motif)) + ' ' + directory + output +
                    ' ' + directory + 'neg_' + output)
                out(cmdLine)
                with open(resdir+'/'+motif+'.log', 'w') as lout:
                    call([cmdLine], shell=True, stdout=lout)
                out("pipeline finished, copying data and results to "+resdir+'/'+motif)
                # file not necessary, 2/3 the space
                try:
                    os.unlink('results/cmfSeeds')
                except:
                    # if the file didn't exist, that means the pipeline didn't run (or at least cmf) and something is wrong
                    pdb.set_trace()
                shutil.copytree('results', resdir+'/'+motif)
                shutil.copyfile(directory+output, resdir+'/'+motif+'.fasta')
                out("done copying")
                # delete fasta files
                for f in filter(lambda a: output.split('.')[0] in a and 'fasta' in a,
                        os.listdir(directory)):
                    out("deleting " + directory+f)
                    os.unlink(directory+f)
                
        except:
            out(sys.exc_info())
        

if __name__ == "__main__":
    main()

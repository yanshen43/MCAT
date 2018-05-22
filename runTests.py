import os
import traceback
import shutil
import time
import subprocess
from re import compile
from sys import argv
# Written by Jeff Robertson (thejar@vt.edu)
def main(argv):
    if len(argv) != 3:
        print "Invocation: %s input_directory output_directory" % argv[0]
        return 0
    count = 0
    pat = compile("^[ACGTacgt]+\.fasta$")
    outDir = argv[2]#"midTestResults71/"
    inDir = argv[1]#"midTestFiles7/"
    #outDir = "testTestResults/"
    #inDir = "testTestFiles/"

    with open("logs/test_"+time.strftime("%y-%m-%d_%H:%M:%S")+".log", 'w') as logFile:
        def out(s=""):
            print s
            print >>logFile, time.asctime(), s
            logFile.flush()
        motifs = []
        for testFile in os.listdir(inDir):
            # only run for positive files
            if not pat.match(testFile):
                continue
            motif = testFile.split('.')[0]
            # don't run with long motifs
            if len(motif) < 7 or len(motif) > 15:
                continue
            motifs.append(motif)
        out("running " + str(len(motifs)) + " tests")
        
        for motif in motifs:
            cmdline = ["./orange_pipeline.py","-w",str(len(motif)),'-v','4',
                    inDir+motif+".fasta",inDir+motif+"_neg.fasta"]
            out()
            out("TEST: " + str(count))
            out(" ".join(cmdline))
            
            with open(outDir+motif+".out", 'w') as fout:
                fout.write(time.asctime() + "\n")
                fout.write(" ".join(cmdline) + "\n\n")
                fout.flush()
                subprocess.call(cmdline, stdout=fout, stderr=fout)
                try:
                    os.unlink('results/cmfSeeds')
                except Exception as e:
                    print e
                    print traceback.print_exc()

            out("DONE: " + str(count))
            out("copying results to "+outDir+motif)
            shutil.copytree("results", outDir+motif)
            count += 1

        out(str(count) + " tests completed")    

if __name__ == '__main__':
    main(argv)

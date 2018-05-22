#! /usr/bin/python
import pdb
from sys import argv
from os import listdir
if len(argv) != 2:
    print "Invocation: %s [directory]"
    quit()
resDir = argv[1]
good = 0
bad = 0
def match(s1, s2):
    return max([sum([1 if c1 == c2 else 0 for c1, c2 in zip(s1, s2[p:]+s2[:p])])
        for p in xrange(len(s2))])
for outfile in listdir(resDir):
    try:
        if outfile[-4:] == '.out':
            motif = outfile[:-4]
            print 'parsing results from', motif+'.fasta',
            with open(resDir+outfile) as oin:
                s = 0
                matched = False
                for line in oin:
                    if line.strip() == 'Finished discovery':
                        s = 1
                    elif s == 1 and len(line.strip()) > 1 and line[0] == '[':
                        #pdb.set_trace()
                        results = map(lambda a:a[1], eval(line))
                        print results
                        matched = False
                        for res in results:
                            if match(res, motif.upper()) >= len(motif)*.75:
                                matched = True
                                print res, 'matches the target motif'
                                break
                        break
                if matched:
                    good += 1
                elif s:
                    bad += 1
                    print 'not found'#motif
            print
    except Exception as e:
        print e
        #print line
        pdb.set_trace()
print good, 'target motifs found out of', good+bad

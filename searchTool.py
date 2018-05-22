#! /usr/bin/python

# Written by Jeff Robertson (thejar@vt.edu)
# originally for CS4884 Spring 2016 Class project

from compareTool import compare, diff
from sys import stdout
import math
import pdb

def searchFasta(fName, resultsDir, motif, motifnum):
    print "~",motifnum,"#"
    sInstances = {0:[], -1:[], -2:[]}
    threshold = 0
    sequenceName = ""
    mLen = len(motif)

    with open(fName, "r") as f:
        fileText = f.read().split("\n")
    fLen = len(fileText)
    lineNum = -1

    while lineNum < fLen - 1:
        # retrieve line
        lineNum += 1
        line = fileText[lineNum]

        if len(line) == 0:
            # ignore blank lines
            continue

        if line[0] == '>':
            # beginning of new sequence, save sequence name
            sequenceName = line[1:].strip()
            continue
        
        # retrieve sequence
        sequence = ""
        while lineNum < fLen and len(fileText[lineNum]) != 0 and fileText[lineNum][0] != '>':
            sequence += fileText[lineNum].strip()
            lineNum += 1
#        print "|" + fileText[lineNum] + "|" + str(lineNum)
        if lineNum < fLen and len(fileText[lineNum]) > 0 and fileText[lineNum][0] == '>':
            lineNum -= 1
        sLen = len(sequence)
        
        # search sequence
        sPos = 0
        while sPos < sLen - mLen:
            match = sequence[sPos:sPos+mLen]
            score = compare(motif, match)
            # if match is better
            if score > min(sInstances.keys()):
                if diff(motif, match) < diff(motif, sequence[sPos+1:sPos+mLen+1]):
                    sPos += 1
                    continue
                if not score in sInstances:
                    # add score to found instances and remove previous min
                    sInstances[score] = []
                    del sInstances[min(sInstances.keys())]
                sInstances[score].append([match, score, sequenceName, str(sLen-sPos), sLen])
                sPos += mLen
            #end if
            sPos += 1
        #end while
    #end while
    print motif + "\n Instances:"
    instances = []
    for key in sInstances:
        for instance in sInstances[key]:
            instances.append(instance)
    
    # write to file for visualisation
    with open(resultsDir + "motif" + str(motifnum) + ".instances", "w") as fout:
        #pdb.set_trace()
        for instance in instances:
            fout.write(">"+instance[2]+"\t"+str(instance[3])+"|"+str(instance[4])\
                    +"\n"+instance[0]+"\n\n")
    print makePWM(mLen, [x[0] for x in instances])[0]
    
# constructs a PWM from the given motif instances
def makePWM(mLen, instances, backgroundRatios=[.25,.25,.25,.25, 1]):
    #pdb.set_trace()
    PWMstr = ""
    logLikelihood = 0
    print "~" + str(instances) + "!"
    letters = {'A':0, 'C':1, 'G':2, 'T':3, 'N':4}
    pwm = [[0.0 for x in xrange(mLen)] for y in xrange(len(letters))]
    inc = 1.0/len(instances)
    
    for instance in instances:
        #fout.write(str(instance)+"\n")
        for pos, char in enumerate(instance):
            if char in letters:
                pwm[letters[char]][pos] += inc
    
    for let in "ACGTN":
        row = let + ": "
        for posVal in pwm[letters[let]]:
            row += "{:.2f} | ".format(posVal)
        PWMstr += row+"\n"
    # currently discounting N in computing loglikelihood
    for col in [[pwm[x][y] for x in range(4)] for y in range(len(pwm[0]))]:
        logLikelihood += math.log(max(col)/backgroundRatios[col.index(max(col))], 2)
#        pdb.set_trace()
    return PWMstr, logLikelihood
#        print pwm



def main():
    fastaFile = "data/positive_genes.fasta"
    searchFasta(fastaFile, "./", "ACTCTATG", 1)

if __name__ == "__main__":
    main()

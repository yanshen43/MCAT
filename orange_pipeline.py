#! /usr/bin/python

# Written (mostly) by Jeff Robertson (thejar@vt.edu)
# with notable help from Yanshen Yang (yanshen@vt.edu)
# originally for CS4884 Spring 2016 Class project

import threading
import subprocess
import time
import argparse
import fasta2bp
import pdb
import struct
import os
import sys
import pprint
import shutil
from compareTool import compare
from searchTool import searchFasta
from visTool import visualizeOutput
from statistics import pval
from filterTool import filterTool
from collections import defaultdict as DD, Counter
from json import dump
import statistics
import statisticspwm
from statfile import StatFile


from math import log
from scipy import stats
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

# define constants and paths, could be read from command line
NUM_OF_TOOLS = 6
MOTIF_LEN = 12
NUM_MOTIFS = 3 #number of motifs expected per sequence
DECOD_ITER = 50 #number of iterations to perform
SEQ_CARD = 20 #****number of total motif occurences (?)****
opmode = str(8*struct.calcsize("P")) + "bit"
config = {'MEME':{}, 'DECOD':{}, 'CMF':{}, 'WEEDER':{}, 'BIOPROSPECTOR':{},\
        'dustmasker':{}, 'filterTool':{}, 'XXmotif':{}}

DATA_DIR = "data/"
POS_SEQ = "mad2.txt"
FILTERED_SEQ = ''
NEG_SEQ = "negative_mad2.txt"
RES_DIR = "results/"

times = {}
foundMotifs = {}
# {"toolname":[[toolResult1, pos1, pos2, ...], ...]}
foundMotifsSeqs = {}
# {"toolname":{toolResults1: {"seq1": [pos1, pos2, ...], "seq2": [pos1, pos2, ...], ...}, toolResult2...}}
mlist = []
probabilities = {"A":0, "C":0, "G":0, "T":0, "N": 0}
FORCE_FILTER = False

DEBUG = False

# functions for running each tool
########################
def runDustmasker():
    times["dustmasker"] = {}
    times["dustmasker"]["Start"] = time.time()
    if not DEBUG:
        with open("results/dustmasker.out", "w") as fout:
            # only filter if file hasn't been filtered in a prev run or FORCE_FILTER is set
            if os.path.exists(config['filterTool']['foname']) and not FORCE_FILTER:
                fout.write('not running dustmasker because filtered fasta file already exists\n')
                fout.write('run pipeline with -ff to force run dustmasker\n')
            else:
                cmdLine = './dustmasker/dustmasker '
                cmdLine += " ".join([opt + ' ' + config["dustmasker"][opt] for opt in config["dustmasker"]])
                fout.write(cmdLine + '\n')
                fout.flush()
                subprocess.call([cmdLine], stdout=fout, shell=True)
                print 'Filtered fasta file with dustmasker'
                with open("results/filterTool.out", "w") as fout:
                    filterTool(config['dustmasker']['-out'], config['filterTool']['foname'], fout)
        FILTERED_SEQ = config['filterTool']['foname']
    times["dustmasker"]["End"] = time.time()

def runCMF():
    times["CMF"] = {}
    times["CMF"]["Start"] = time.time()
    if not DEBUG:
        with open("results/cmf.out", "w") as fout:
            cmdLine = './cmf/cmf '
            cmdLine += " ".join([opt + ' ' + config['CMF'][opt] for opt in config['CMF']])
            fout.write(cmdLine + '\n')
            fout.flush()
            subprocess.call([cmdLine], stdout=fout, shell=True)
    times["CMF"]["End"] = time.time()
    print "CMF finished in <= {:.3f} seconds".format(times["CMF"]["End"] - times["CMF"]["Start"])

def runWeeder():
    times["Weeder"] = {}
    times["Weeder"]["Start"] = time.time()
    if not DEBUG:
        with open(RES_DIR+"weeder.out", "w") as fout:
            cmdLine = './weeder/weeder2 '
            cmdLine += " ".join([opt + ' ' + config['WEEDER'][opt] for opt in config['WEEDER']])
            fout.write(cmdLine + '\n')
            fout.flush()
            subprocess.call([cmdLine], stdout=fout, stderr=fout, shell=True)
    times["Weeder"]["End"] = time.time()
    print "Weeder finished in <= {:.3f} seconds".format(times["Weeder"]["End"] - times["Weeder"]["Start"])

def runDECOD():
    times["DECOD"] = {}
    times["DECOD"]["Start"] = time.time()
    if not DEBUG:
        with open("results/decod.out", "w") as fout:
            cmdLine = 'java -jar ./DECOD/DECOD-20111024.jar '
            cmdLine += " ".join([opt + ' ' + config['DECOD'][opt] for opt in config['DECOD']])
            fout.write(cmdLine + '\n')
            fout.flush()
            subprocess.call([cmdLine], stdout=fout, shell=True)
    times["DECOD"]["End"] = time.time()
    print "DECOD finished in <= {:.3f} seconds".format(times["DECOD"]["End"] - times["DECOD"]["Start"])

def runMEME():
    times["MEME"] = {}
    times["MEME"]["Start"] = time.time()
    if not DEBUG:
        with open("results/meme.out", "w") as fout, open("results/meme/meme.txt", 'w') as mout:
            cmdLine = './meme/' + opmode + '/bin/meme ' + FILTERED_SEQ + ' '
            cmdLine += " ".join([opt + ' ' + config['MEME'][opt] for opt in config['MEME']])
            fout.write(cmdLine + '\n')
            fout.flush()
            subprocess.call([cmdLine], stdout=mout, stderr=fout, shell=True)
    times["MEME"]["End"] = time.time()
    print "MEME finished in <= {:.3f} seconds".format(times["MEME"]["End"] - times["MEME"]["Start"])

def runBioProspector():
    times["BioProspector"] = {}
    times["BioProspector"]["Start"] = time.time()
    if not DEBUG:
        with open("results/bioprospector.out", "w") as fout:
            fasta2bp.convert(FILTERED_SEQ, "bp_files/bp_"+get_filename(FILTERED_SEQ))
            fasta2bp.convert(NEG_SEQ, "bp_files/bp_"+get_filename(NEG_SEQ))
            cmdLine = './BioProspector/' + opmode + '/BioProspector '
            cmdLine += " ".join([opt + ' ' + config['BIOPROSPECTOR'][opt] for opt in config['BIOPROSPECTOR']])
            fout.write(cmdLine + '\n')
            fout.flush()
            subprocess.call([cmdLine], stdout=fout, stderr=fout, shell=True)
    times["BioProspector"]["End"] = time.time()
    print "BioProspector finished in <= {:.3f} seconds".format(times["BioProspector"]["End"] - times["BioProspector"]["Start"])

def runXXmotif():
    times['XXmotif'] = {}
    times['XXmotif']['Start'] = time.time()
    if not DEBUG:
        with open("results/XXmotif.out", "w") as fout:
            cmdLine = './XXmotif/XXmotif'
            cmdLine += " ".join([opt + ' ' + config['XXmotif'][opt] for opt in config['XXmotif']])
            fout.write(cmdLine + '\n')
            fout.flush()
            subprocess.call([cmdLine], stdout=fout, stderr=fout, shell=True)
    times['XXmotif']['End'] = time.time()
    print "XXmotif finished in <= {:.3f} seconds".format(times["XXmotif"]["End"] - times["XXmotif"]["Start"])

# runs all tools in parallel
########################
def runTools():
    t = [threading.Thread() for x in range(NUM_OF_TOOLS)]
    t[0].run = runCMF
    t[1].run = runWeeder
    t[2].run = runDECOD
    t[3].run = runMEME
    t[4].run = runBioProspector
    t[5].run = runXXmotif

    for th in t:
        th.start()
    print "Started pipeline"
    for th in t:
        th.join()
    print "Finished discovery"



# parses the output from each tool and combines it into consistantly formatted file
########################

def parseCMF(fout):
    try:
        with open(RES_DIR + "output.txt", "r") as CMFResults:
            foundMotifs["CMF"] = []
            foundMotifsSeqs["CMF"] = DD(lambda: DD(list))
            alreadyFound = {}
            readmode = 0 # 0=looking for motif,1=reading pos,2=looking for pos
            for line in CMFResults:
                if readmode == 2:
                    readmode = 1
                    continue
                if readmode == 1:
                    if len(line) < 2:
                        readmode = 0
                        continue
                    seqName, mpos = line.split("\t")[:2]
                    seqName = seqName.strip()
                    #pdb.set_trace()
                    foundMotifs["CMF"][-1][1].append(SEQ_LENS[seqName[1:]] - int(mpos))
                    foundMotifsSeqs["CMF"][motif][seqName[1:]].append(int(mpos))
                    fout.write("\t"+str(SEQ_LENS[seqName[1:]] - int(mpos)))
                    readmode = 2
                    continue
                if line[0:7] == "MOTIF:\t" and not line[7:] in alreadyFound:
                    motif = line[7:-1]
                    fout.write("\nCMF\t"+line[7:-1])
                    foundMotifs["CMF"] += [[line[7:-1], []]]
                    #foundMotifsSeqs["CMF"][motif] = DD(list)
                    #alreadyFound[line[7:]] = 0
                elif "Positive Sites:" in line:
                    readmode = 1
    except IOError:
        print "Error opening CMF's results file"

def parseWeeder(fout):
    try:
        with open(RES_DIR + get_filename(FILTERED_SEQ) + ".w2", "r") as _weederResults:
            foundMotifs["Weeder"] = []
            foundMotifsSeqs["Weeder"] = {}
            weederResults = _weederResults.readlines()[6:]
            for line in weederResults:
                if "Matrix" in line:
                    result = line.split()
                    motif = result[2]
                    fout.write("\nweeder\t"+result[2])
                    foundMotifs["Weeder"] += [[result[2], []]]
                    foundMotifsSeqs["Weeder"][motif] = DD(list)
                elif line[0] == ">":
                    seqName, mpos = line.split('\t')[0::4]
                    seqName = seqName.strip()
                    foundMotifs["Weeder"][-1][1].append(SEQ_LENS[seqName[1:]] - int(mpos))
                    foundMotifsSeqs["Weeder"][motif][seqName[1:]].append(int(mpos))
                    fout.write("\t"+str(SEQ_LENS[seqName[1:]] - int(mpos)))
    except IOError:
        print "Error opening Weeder's results file"

def parseDECOD(fout):
    try:
        with open(RES_DIR + "decod_found_motifs.txt", "r") as _decodResults:
            foundMotifs["DECOD"] = []
            foundMotifsSeqs["DECOD"] = {}
            alreadyFound = {}
            decodResults = _decodResults.readlines()
            lineNum = 0
            chars = ["A", "C", "G", "T"]
            for motif_num in xrange(int(config["DECOD"]["-nmotif"])):
                PWM = {}
                # parse PWM
                for line in decodResults[lineNum:]:
                    lineNum+=1
                    if len(line) > 1 and line[0] in chars:
                        PWM[line[0]] = line.strip("ACGT []\n").split()
                        if line[0] == "T":
                            break
                # construct motif from PWM
                motif = ""
                for i in xrange(len(PWM['A'])):
                    col = [PWM[x][i] for x in chars]
                    motif += chars[col.index(max(col))]
                
                foundMotifsSeqs["DECOD"][motif] = DD(list)
                # parse instances
                if not motif in alreadyFound:
                    fout.write("\nDECOD\t"+motif)
                    # find beginning of instances
                    while decodResults[lineNum][0] != ">":
                        lineNum+=1

                    # save instances
                    positions = []
                    while decodResults[lineNum][0] == ">":
                        line = decodResults[lineNum]
                        lineNum+=1
                        if "|revcom" in line:
                            continue
                        seqName = line[1:line.find("\t")]
                        seqName = seqName.strip()
                        pos = int(line.split("\t")[1])
                        foundMotifsSeqs["DECOD"][motif][seqName].append(pos)
                        pos = SEQ_LENS[seqName] - pos
                        positions.append(pos)
                        fout.write("\t" + str(pos))
                    foundMotifs["DECOD"] += [[motif, positions]]
                    alreadyFound[motif] = 0
    except IOError:
        print "Error opening DECOD's results file"

def parseMEME(fout):
    try:
        with open(RES_DIR + "/meme/meme.txt", "r") as memeResults:
            foundMotifs["MEME"] = []
            foundMotifsSeqs["MEME"] = {}
            readMode = 0 # 0 is looking for next motif, 1 is looking for pos's, 2 is reading pos's
            for line in memeResults:
                if readMode == 1:
                    readMode = 2
                    continue
                if readMode == 2:
                    if "----------------------" in line:
                        readMode = 0
                        continue
                    seqName, sPos = line.split()[0:2]
                    seqName = seqName[0:19].strip()
                    foundMotifs["MEME"][-1][1].append(SEQ_LENS[seqName] - int(sPos))
                    foundMotifsSeqs["MEME"][motif][seqName].append(int(sPos))
                    fout.write("\t" + str(SEQ_LENS[seqName] - int(sPos)))
                if "Multilevel" in line:
                    motif = line.strip().replace("Multilevel","").replace(" ","")
                    fout.write("\nMEME\t"+motif)
                    foundMotifs["MEME"] += [[motif, []]]
                    foundMotifsSeqs["MEME"][motif] = DD(list)
                elif "Sequence name" in line and "Start" in line:
                    readMode = 1
    except IOError:
        print "Error opening MEME's results file"


def parseBioProspector(fout):
    # read sequence mappings
    mappings = {}
    mfname = FILTERED_SEQ.split('/')[-1]
    with open("bp_files/bp_" + mfname + ".mappings", "r") as mapIn:
        for line in mapIn:
            line = line.split(">>")
            mappings[line[0]] = line[1].strip()
    with open(RES_DIR + "bp_output.txt", "r") as bpResults:
        foundMotifs["BioProspector"] = []
        foundMotifsSeqs["BioProspector"] = {}
        for line in bpResults:
            if "Motif #" in line:
                motif = line[line.index("(")+1:line.index("/")]
                fout.write("\nBioProspector\t"+motif)
                foundMotifs["BioProspector"] += [[motif, []]]
                foundMotifsSeqs["BioProspector"][motif] = DD(list)
            elif ">seq" in line:
                slen = int(line.split()[2])
                mpos = int(line.split()[-1])
                if line.split()[-2] == 'r':
                    mpos -= len(motif)
                foundMotifs["BioProspector"][-1][1].append(slen - mpos)
                foundMotifsSeqs["BioProspector"][motif][mappings[line.split()[0][1:]]].append(mpos)
                fout.write("\t" + str(slen - mpos))

def parseXXmotif(fout):
    try:
        foundMotifs['XXmotif'] = []
        foundMotifsSeqs['XXmotif'] = {}
        reportedMotifs = [None]
        lets = 'ACGT'
        baseFname = RES_DIR + '.'.join(FILTERED_SEQ.split('/')[-1].split('.')[:-1])
        #baseFname = 'XXmotif/results/hsap_core_promoters_all'
        # parse motifs
        #pdb.set_trace()
        with open(baseFname+'.pwm') as XXmotifResults:
            cm = []
            for i, line in enumerate(XXmotifResults):
                if i % 6 == 0 or i % 6 == 5:
                    if len(cm) == 4:
                        motif = ''.join([lets[p.index(max(p))] for p in zip(cm[0], cm[1], cm[2], cm[3])])
                        foundMotifs['XXmotif'].append([motif, []])
                        foundMotifsSeqs['XXmotif'][motif] = DD(list)
                    cm = []
                else:
                    cm.append(map(float, line.split()))
        seqs = []
        # parse sequence mappings
        with open(baseFname+'_sequence.txt') as XXmotifSequences:
            for i in xrange(4):
                XXmotifSequences.next()
            for line in XXmotifSequences:
                seqs.append(line.split('\t')[-1].strip())
        # parse instance locations
        with open(baseFname+'_Pvals.txt') as XXmotifPoss:
            for i in xrange(4):
                XXmotifPoss.next()
            mnum = 0
            for line in XXmotifPoss:
                if len(line) < 5:
                    continue
                if line[:6] == 'Motif ':
                    mnum = int(line.split()[1][:-1])-1
                else:
                    pos = int(line.split('\t')[4])
                    seq = seqs[int(line.split('\t')[3])-1]
                    foundMotifs['XXmotif'][mnum][1].append(pos)
                    foundMotifsSeqs['XXmotif'][foundMotifs['XXmotif'][mnum][0]][seq].append(pos)
        for motif, poss in foundMotifs['XXmotif']:
            fout.write('\nXXmotif\t'+motif)
            for pos in poss:
                fout.write('\t'+str(pos))
    except IOError as e:
        print e, e.filename
        print "Error opening XXmotif's result file"

# calls the parse functions for each tool
# consolidates results from tool outputs to results.out
def parseResults():
    with open("results/results.out", "w") as fout:
        parseCMF(fout)
        parseWeeder(fout)
        parseDECOD(fout)
        parseMEME(fout)
        parseBioProspector(fout)
        parseXXmotif(fout)
    # dump raw results to json
    with open("results/results.json", "w") as jout:
        dump(foundMotifsSeqs, jout, indent=2)

def match(s1, s2):
    l1, l2 = len(s1), len(s2)
    for i in xrange(min(l1, l2)):
        if s1[i] != s2[i]:
            return 0
    if l1 == l2:
        return l1
    return min(l1, l2) - 0.5

def best(lst, s):
    """best(lst, s) Returns the item in lst which has the largest prefix equal to s"""
    besti = ''
    bestv = -1
    for item in lst:
        v = match(item, s)
        if v > bestv:
            besti = item
            bestv = v
    return besti

# poll the results from all tools and sort them by their score
def voteRank(sequences, motifs):
    poll = {}
    for seq in sequences:
        poll[seq] = [0.0] * len(sequences[seq])
    
    # perform poll
    for tool in motifs:
        for motif in motifs[tool]:
            for seq in motifs[tool][motif]:
                sequence = best(sequences, seq)
                for pos in motifs[tool][motif][seq]:
                    for i in xrange(pos, pos + len(motif)):
                        try:
                            # instead of weighting all results the same (1), we
                            # could bias based on tool or number of results or something like that
                            poll[sequence][i - 1] += 1
                        except Exception as e:
                            print e
                            print 'It appears a tool has reported finding a motif',\
                                'outside the bounds of a sequence'
                            print 'such as finding a motif of length 10 at position',\
                                '195 in a sequence with length 200'
                            pdb.set_trace()
    # add up votes for each motif
    ress = DD(int)
    for tool in motifs:
        for motif in motifs[tool]:
            for seq in motifs[tool][motif]:
                for pos in motifs[tool][motif][seq]:
                    for p in xrange(pos, pos + len(motif)):
                        ress[motif] += poll[best(sequences, seq)][p-1]
    # sort motifs by number of votes
    return sorted(map(lambda a: list(a[::-1]), ress.iteritems()))

# framework for ranking motifs
def voteRankMotifs():
    seqs = {}
    with open(POS_SEQ, 'r') as fin:
        seqName = ""
        for line in fin:
            if line[0] == '>':
                seqName = line.strip()[1:]
                seqs[seqName] = ""
            else:
                seqs[seqName] = seqs[seqName] + line.strip()
    # actually do the ranking
    results = voteRank(seqs, foundMotifsSeqs)

    # take the top result
    mlist.append(results[-1])
    for score, newMotif in results[::-1][1:]:
        # try the next highest ranked motif
        maxSimilarity = 0
        for oldMotif in mlist:
            maxSimilarity = max(maxSimilarity, compare(newMotif, oldMotif[1]))
        if maxSimilarity < 4:
            # include it if it isn't too similar to previously found motifs
            mlist.append([score, newMotif])
        # stop if requested number of results have been collected
        if len(mlist) > NUM_VIS - 1:
            break
    print mlist
    for mnum, m in enumerate(mlist):
        searchFasta(FILTERED_SEQ, RES_DIR, m[-1], mnum+1)

# not used, old clustering/ranking method
def combineMotifs():
    print "Combining motifs"
    results = []
    for tool in foundMotifs:
        others = list(set(foundMotifs) - set([tool]))
        for motif in foundMotifs[tool]:
            if 'TATATA' in motif[0]:
                continue
            mscore = 0
            compares = 0
            for oTool in others:
                for oMotif in foundMotifs[oTool]:
                    compares += 1
                    mscore += compare(motif[0], oMotif[0], pos1=motif[1], pos2=oMotif[1])
            results += [[mscore/compares, motif[0]]]
    results.sort()
    with open("results/compareScores.out", 'w') as cout:
        for score, motif in results:
            cout.write(motif + '\t%.3f\n'%score)
    mlist.append(results[-1])
    for score, newMotif in results[::-1][1:]:
        maxSimilarity = 0
        for oldMotif in mlist:
            maxSimilarity = max(maxSimilarity, compare(newMotif, oldMotif[1]))
        if maxSimilarity < 4:
            mlist.append([score, newMotif])
        # stop if requested number of results have been collected
        if len(mlist) > NUM_VIS - 1:
            break
    print mlist
    for mnum, m in enumerate(mlist):
        searchFasta(FILTERED_SEQ, RES_DIR, m[-1], mnum+1)
    return results

def get_filename(s):
    index = s.rfind('/')
    if index < 0:
        return s
    else:
        return s[index+1:]

# goes through positive fastafile and gets the length of each sequence
# stores in global SEQ_LENS by seq name
def getSeqLens():
    with open(FILTERED_SEQ, "r") as fastafile:
        global SEQ_LENS
        SEQ_LENS = {}
        currSeqLen = 0
        for line in fastafile:
            if line[0] == '>':
                if currSeqLen != 0 and not seqName in SEQ_LENS:
                    SEQ_LENS[seqName] = currSeqLen
                    SEQ_LENS[seqName[0:19]] = currSeqLen
                seqName = line[1:].strip()
                currSeqLen = 0
            elif len(line) > 1:
                currSeqLen += len(line)
            else:
                SEQ_LENS[seqName] = currSeqLen
                #beacuse weeder only likes the first 24 chars of the seq name
                SEQ_LENS[seqName[0:19]] = currSeqLen 
                # and meme only uses the first token of the seq name
                SEQ_LENS[seqName.split()[0].strip()] = currSeqLen
        SEQ_LENS[seqName] = currSeqLen
        SEQ_LENS[seqName[0:19].strip()] = currSeqLen
        SEQ_LENS[seqName.split()[0].strip()] = currSeqLen
           
# goes through positive fastafile and gets background probabilities
def getProbabilities():
    with open(FILTERED_SEQ, "r") as fastafile:
        lets = "ACGTN"
        global probabilities
        changeProbs = {x:{y:0 for y in lets} for x in lets}
        changeProbs['>'] = {x:0 for x in lets}
        last = ''
        seqNames = Counter()
        for line in fastafile:
            if line[0] == '>':
                seqName = line[1:].strip()
                seqNames.update([seqName])
            if line[0] in lets:
                for let in line.strip():
                    probabilities[let] += 1
                    changeProbs[last][let] += 1
                    last = let
            else:
                last = line[0]
        s = sum(probabilities.values())
        for let in lets:
            if float(s) > 0:
                probabilities[let] = probabilities[let]/float(s)
            else:
                probabilities[let] = 0
        for l1 in changeProbs:
            s = 0.0
            for l2 in changeProbs[l1]:
                s += changeProbs[l1][l2]
            if s > 0.1:
                for l2 in changeProbs[l1]:
                    changeProbs[l1][l2] /= s
        print probabilities
        with open(RES_DIR + "matrix.txt", 'w') as mout:
            for let in lets:
                mout.write(str(changeProbs['>'][let]) + '\n')
            for l1 in lets:
                for l2 in lets:
                    mout.write(str(changeProbs[l1][l2]) + '\n')
        print changeProbs
        # check for sequence names which are duplicates
        if seqNames.most_common(1)[0][1] > 1:
            print "Error: multiple sequences in the positive fasta file have the same name"
            for seqName, count in seqNames.most_common():
                if count < 2:
                    break
                print count, seqName
            quit(1)

# generates p-Values for the list of motifs in mlist
def getPValues(motifs, sequences,pvals=None,method=None):
    print "Generating p-Values"
    pvalues = []
    for motif in motifs:
        if method is None:
            if pvals is not None:
                p_val = pval(motif,sequences,pvals=pvals)
            else:
                p_val = pval(motif,sequences)
        else:
            if pvals is not None:
                p_val = pval(motif,sequences,pvals=pvals,method=method)
            else:
                p_val = pval(motif,sequences,method=method)
        print motif, p_val
        pvalues.append( (motif, p_val) )
    return pvalues

# generates p-Values using QFAST
def getPValuesQFAST(motifs, sequences):
    print "Generating p-Values using QFAST"
    pvalues = []
    for motif in motifs:
        pwm = statisticspwm.searchFasta(FILTERED_SEQ, RES_DIR, motif)
        p_val = statisticspwm.pval(motif,sequences, pwm)
        print motif, p_val
        pvalues.append( (motif, p_val) )
    return pvalues

def computeZScores(sf):
    sys.stdout.write('Computing Z-scores')
    sys.stdout.flush()
    begin = time.time()
    with open('results/zscores.txt','w') as f:
        zscores = statistics.zscore_all(MOTIF_LEN,FILTERED_SEQ)
        pprint.pprint(zscores, f)
    end = time.time()
    sys.stdout.write('...Done. (%f sec)\n'%(end-begin))

    sf.add_vals('zscore',zscores)
    return zscores

def computePVals(sf,generateAll=False):
    # filtered data
    sys.stdout.write('Computing p-Values')
    sys.stdout.flush()
    begin = time.time()
    with open('results/pvals.txt','w') as f:
        pvs_out = statistics.pval_all(MOTIF_LEN,FILTERED_SEQ)
        pvs = [(k,v) for k,v in pvs_out.iteritems()]
        pvs = sorted(pvs, key=lambda x: float(x[1])) # sort by pvalue
        pprint.pprint(pvs,f)
    end = time.time()
    sys.stdout.write('...Done. (%f sec)\n'%(end-begin))

    sf.add_vals('pval',pvs_out)

    if generateAll:
        sys.stdout.write('Computing p-Values (Second-order Markov)')
        sys.stdout.flush()
        begin = time.time()
        with open('results/pvals_order2.txt','w') as f:
            pvs_out2 = statistics.pval_all(MOTIF_LEN,FILTERED_SEQ,2)
            pvs2 = [(k,v) for k,v in pvs_out2.iteritems()]
            pvs2 = sorted(pvs2, key=lambda x: float(x[1])) # sort by pvalue
            pprint.pprint(pvs2,f)
        end = time.time()
        sys.stdout.write('...Done. (%f sec)\n'%(end-begin))

        sf.add_vals('pval_2',pvs_out2)

        # unfiltered data
        sys.stdout.write('Computing p-Values [UNFILTERED]')
        sys.stdout.flush()
        begin = time.time()
        with open('results/pvals_unfiltered.txt','w') as f:
            pvs_out_u = statistics.pval_all(MOTIF_LEN,POS_SEQ)
            pvs = [(k,v) for k,v in pvs_out_u.iteritems()]
            pvs = sorted(pvs, key=lambda x: float(x[1])) # sort by pvalue
            pprint.pprint(pvs,f)
        end = time.time()
        sys.stdout.write('...Done. (%f sec)\n'%(end-begin))

        sys.stdout.write('Computing p-Values (Second-order Markov) [UNFILTERED]')
        sys.stdout.flush()
        begin = time.time()
        with open('results/pvals_order2_unfiltered.txt','w') as f:
            pvs_out2_u = statistics.pval_all(MOTIF_LEN,POS_SEQ,2)
            pvs = [(k,v) for k,v in pvs_out2_u.iteritems()]
            pvs = sorted(pvs, key=lambda x: float(x[1])) # sort by pvalue
            pprint.pprint(pvs,f)
        end = time.time()
        sys.stdout.write('...Done. (%f sec)\n'%(end-begin))

    return pvs_out

def computeQVals(sf):
    sys.stdout.write('Computing q-Values')
    sys.stdout.flush()
    begin = time.time()
    with open('results/qvals.txt','w') as f:
        qvals = statistics.qval_all(MOTIF_LEN,FILTERED_SEQ)
        pprint.pprint(qvals,f)
    end = time.time()
    sys.stdout.write('...Done. (%f sec)\n'%(end-begin))

    sf.add_vals('qval',qvals)

def computePValEsts(sf, pvs_out):
    sys.stdout.write('Computing p-Value estimates (QFAST,Fisher,Stouffer)')
    sys.stdout.flush()
    begin = time.time()
    with open('results/pvals_est.txt','w') as f:
        pvals_est = statistics.pval_est(pvs_out,top=50)
        pprint.pprint(pvals_est,f)
    end = time.time()
    sys.stdout.write('...Done. (%f sec)\n'%(end-begin))

    pv_q = {x:pvals_est[x][0] for x in pvals_est}
    sf.add_vals('p_qfast',pv_q)

    pv_f = {x:pvals_est[x][1] for x in pvals_est}
    sf.add_vals('p_fisher',pv_f)

    pv_s = {x:pvals_est[x][2] for x in pvals_est}
    sf.add_vals('p_stouffer',pv_s)


def scoreMotifs():
    sf_filename = 'stat_'+str(int(time.time()*100))
    sf = StatFile(STAT_DIR+'/'+sf_filename)

    computeZScores(sf)
    pvs_out = computePVals(sf)
    computeQVals(sf)
    computePValEsts(sf, pvs_out)

    sf.save()

    motifs_comb = [x[1] for x in mlist]
    pvals_comb = getPValues(motifs_comb,FILTERED_SEQ)
    #pvals_comb_qfast = getPValuesQFAST(motifs_comb, FILTERED_SEQ)

    for mnum, (motif, pval) in enumerate(pvals_comb):
        mlist[mnum] += [pval, statistics.zscore(motif, FILTERED_SEQ)]
    print mlist


    # analyze correlations
    zscores = [0 for x in mlist]
    motifs = [0 for x in mlist]
    pvals = [0 for x in mlist]
    ranks = [0 for x in mlist]
    for x in range(len(mlist)):
        zscores[x] = mlist[x][3]
        motifs[x] = mlist[x][1]
        pvals[x] = 1#log(mlist[x][2])
        ranks[x] = x+1


    # compare rankings to p-values
    fig = plt.figure(figsize=(7,7))
    ax = fig.add_subplot(111)

    ax.scatter(ranks,pvals)
    ax.set_xlabel('Rankings')
    ax.set_ylabel('log of p-Values')
    ax.autoscale(tight=True)

    # compute spearman rank correlation
    print('Computing spearman rank correlation for Rankings vs p-Values...')
    corr = stats.spearmanr(ranks,pvals)
    print('   '+str(corr))
    plt.title('spearman = %f'%corr[0])

    filename = 'stat/rankings_vs_pvalues.png'
    print('Generating graph and writing to "'+str(filename)+'"')
    fig.savefig(filename)


    # compare rankings to z scores
    fig = plt.figure(figsize=(7,7))
    ax = fig.add_subplot(111)

    ax.scatter(ranks,zscores)
    ax.set_xlabel('Rankings')
    ax.set_ylabel('Z Scores')
    ax.autoscale(tight=True)

    # compute spearman rank correlation
    print('Computing spearman rank correlation for Rankings vs Z scores...')
    corr = stats.spearmanr(ranks,zscores)
    print('   '+str(corr))
    plt.title('spearman = %f'%corr[0])

    filename = 'stat/rankings_vs_zscores.png'
    print('Generating graph and writing to "'+str(filename)+'"')
    fig.savefig(filename)


    # compare z scores to pvals
    fig = plt.figure(figsize=(7,7))
    ax = fig.add_subplot(111)

    ax.scatter(zscores,pvals)
    ax.set_xlabel('Z scores')
    ax.set_ylabel('log of p-Values')
    ax.autoscale(tight=True)

    # compute spearman rank correlation
    print('Computing spearman rank correlation for Rankings vs Z scores...')
    corr = stats.spearmanr(zscores,pvals)
    print('   '+str(corr))
    plt.title('spearman = %f'%corr[0])

    filename = 'stat/zscores_vs_pvals.png'
    print('Generating graph and writing to "'+str(filename)+'"')
    fig.savefig(filename)

    return sf

# outputs p-Values to a file, sorted by p-Value ascending
# params: pvalues=array of tuples (motif,pval)
#         filename=string for output file name
#         threshold=cutoff for the maximum allowed p-Value
def outputPValues(pvalues, filename, threshold):
    fout = open(filename, 'w')

    pvalues = sorted(pvalues, key=lambda x: float(x[1])) # sort by pvalue
    for x in pvalues:
        motif = x[0]
        pvalue = x[1]
        if pvalue < threshold:
            fout.write(str(motif)+' '+str(pvalue)+'\n')

# create folders for final and intermediate results
# if they don't already exist
def makeFolders():
    if os.path.exists("results/") and not DEBUG:
        shutil.rmtree(os.getcwd()+"/results")
        os.mkdir("results")
        os.mkdir("results/meme")
    elif not os.path.exists("results/"):
        os.mkdir("results")
        os.mkdir("results/meme")
    if not os.path.exists("bp_files/"):
        os.mkdir("bp_files")
    if not os.path.exists("stat/"):
        os.mkdir("stat")

def defaultOptions():
    # set default CMF command line options
    config['CMF']['-i1'] = FILTERED_SEQ
    config['CMF']['-i2'] = NEG_SEQ
    config['CMF']['-l'] = str(MOTIF_LEN - 2)
    config['CMF']['-u'] = str(MOTIF_LEN + 2)
    config['CMF']['-w'] = str(MOTIF_LEN)
    config['CMF']['-t'] = str(5)
    config['CMF']['-f'] = RES_DIR.strip('/')
    config['CMF']['-o'] = RES_DIR + 'cmfSeeds'
    
    # set default weeder command line options
    config['WEEDER']['-f'] = FILTERED_SEQ
    config['WEEDER']['-O'] = 'SC'
    #config['WEEDER']['-b'] = '3'
    config['WEEDER']['-maxm'] = '10'

    # set default DECOD command line options
    config['DECOD']['-pos'] = FILTERED_SEQ
    config['DECOD']['-neg'] = NEG_SEQ
    config['DECOD']['-w'] = str(MOTIF_LEN)
    config['DECOD']['-nmotif'] = str(NUM_MOTIFS)
    config['DECOD']['-niter'] = str(DECOD_ITER)
    config['DECOD']['-c'] = str(SEQ_CARD)
    config['DECOD']['-o'] = RES_DIR + "decod_found_motifs.txt"
    config['DECOD']['-nogui'] = ''
    
    # set default MEME command line options
    config['MEME']['-oc'] = RES_DIR + 'meme'
    config['MEME']['-nmotifs'] = str(NUM_MOTIFS)
    config['MEME']['-w'] = str(MOTIF_LEN)
    config['MEME']['-dna'] = ''
    config['MEME']['-text'] = ''
    
    # set default bioprospector command line options
    config['BIOPROSPECTOR']['-W'] = str(MOTIF_LEN)
    config['BIOPROSPECTOR']['-i'] = 'bp_files/bp_' + get_filename(FILTERED_SEQ)
    config['BIOPROSPECTOR']['-b'] = 'bp_files/bp_' + get_filename(NEG_SEQ)
    config['BIOPROSPECTOR']['-o'] = RES_DIR + 'bp_output.txt'
    config['BIOPROSPECTOR']['-n'] = str(DECOD_ITER)

    # set default dustmasker command line options
    p1 = POS_SEQ.rfind('.')
    config['dustmasker']['-in'] = POS_SEQ
    config['dustmasker']['-out'] = POS_SEQ[:p1] + '.marked' + POS_SEQ[p1:]
    config['dustmasker']['-outfmt'] = 'fasta'
    config['filterTool']['foname'] = FILTERED_SEQ

    # set default XXmotif command line options
    config['XXmotif'][''] = RES_DIR + ' ' + FILTERED_SEQ
    config['XXmotif']['--no-graphics'] = ''
    config['XXmotif']['--zoops'] = ''
    config['XXmotif']['--format'] = 'FASTA'
    config['XXmotif']['--type'] = 'ALL'
    #config['XXmotif']['--localization'] = ''
    config['XXmotif']['--merge-motif-threshold'] = 'LOW'
    config['XXmotif']['--negSet'] = NEG_SEQ

# parse a config file and save the options listed therein
def parseConfig(filename):
    try:
        with open(filename, 'r') as cin:
            tool = ''
            for line in cin:
                # change which tool the options will be applied to
                if line[:5].lower() == 'tool:':
                    tool = line.strip().upper()[5:]
                # delete a default option
                elif line[0] == '!':
                    del config[tool][line.strip()[1:]]
                # add (or replace) an option
                else:
                    option, value = line.strip().split('=')
                    config[tool][option] = value
    except IOError:
        print "Error reading from config file: '" + filename + "'!"
        quit()

# parse and save command line args, parse config file if one is specified
def parseArgs():
    opt = argparse.ArgumentParser()
    opt.add_argument('positive_seq', default='', help='the fasta file for the positive sequences')
    opt.add_argument('negative_seq', default='', help='the fasta file for the negative sequences')
    opt.add_argument('-w', default=12, type=int, help='motif width')
    #opt.add_argument('-o', default='', type=str, help='output directory') # output dir is hardcoded in Weeder...
    opt.add_argument('--nmotifs', default=3, type=int, help='max number of distinct motifs to look for')
    opt.add_argument('--iter', default=5, type=int, help='max number of iterations (for DECOD tool)')
    opt.add_argument('-c', default='', help='config filename')
    opt.add_argument('-s', default='stat', help='statfile directory')
    opt.add_argument('-d', action='store_true', help='turn debug mode on')
    opt.add_argument('-ff', action='store_true', help='force (re)filtering')
    opt.add_argument('-v', default=3, type=int, help='number of results to include in visual output')
    args=opt.parse_args()

    global MOTIF_LEN
    global NUM_MOTIFS
    global DECOD_ITER 
    global SEQ_CARD
    global POS_SEQ
    global NEG_SEQ
    global DEBUG
    global FORCE_FILTER
    global FILTERED_SEQ
    global STAT_DIR
    global NUM_VIS

    MOTIF_LEN = args.w
    NUM_MOTIFS = args.nmotifs #number of motifs expected per sequence
    DECOD_ITER = args.iter #number of iterations to perform
    SEQ_CARD = 5 #****number of total motif occurences (?)****
    POS_SEQ = args.positive_seq # "mad2.txt"
    NEG_SEQ = args.negative_seq # "negative_mad2.txt"
    DEBUG = args.d
    FORCE_FILTER = args.ff
    STAT_DIR = args.s if args.s[-1] != '/' else args.s[:-1] # trim slash
    NUM_VIS = args.v
    
    p1 = POS_SEQ.rfind('.')
    FILTERED_SEQ = POS_SEQ[:p1] + '.filtered'+ POS_SEQ[p1:]
    
    # populate tool's options lists with defaults
    defaultOptions()

    # override (or augment) options from config file
    if args.c != '':
        parseConfig(args.c)

def main():
    parseArgs()
    if DEBUG:
        print 'WARNING: DEBUG MODE IS ON, NO TOOLS WILL RUN!'
        print >>sys.stderr, 'WARNING: DEBUG MODE IS ON, NO TOOLS WILL RUN!'

    # setup
    makeFolders()
    runDustmasker()
    getProbabilities()
    getSeqLens()

    # motif finding
    runTools()

    # post-processing
    parseResults()
    #rankings = combineMotifs()
    voteRankMotifs()

    sf = scoreMotifs()

    # reporting/visualization
    visualizeOutput(RES_DIR, mlist, probabilities)

if __name__ == "__main__":
    main()

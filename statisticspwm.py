# statisticspwm.py
# author: Zhen Guo
# description: p value calculation using QFAST

from math import *
from compareTool import *
from scipy import stats,interpolate
#from multiprocessing import Pool
import sys
sys.setrecursionlimit(1000000000)

# p-Value
#
# parameters:
#     motif=a string for the motif to be evaluated
#     sequences=a string for the filename of a .fasta file with sequences
#     pwm=the position weight matrix of motif
def pval(motif,sequences, pwm):
	mlen = len(motif)	
	with open(sequences, "r") as f:
		text = f.read()
		text = text.replace(' ','')
		text = text.split("\n")
	file_length = len(text)
	ln_num = -1	
	num_seq = 0

	#compute the parameters for qfast calculation:
	letters = ['A', 'C', 'G', 'T', 'N']
	p_fast = 1
	x = {}
	for i in range(mlen+1):
		x[i] = []
	x_compute(0, 0.0, pwm, x)
	M = m_compute(mlen, x, pwm)

	while ln_num < file_length - 1:
		# retrieve line
		ln_num += 1
		line = text[ln_num]

		if len(line) == 0:
			# ignore blank lines
			continue

		if line[0] == '>':
			# beginning of new sequence, save sequence name
			seq_name = line[1:]
			continue

		# retrieve sequence
		sequence = ""
		while ln_num < file_length and len(text[ln_num]) != 0:
			sequence += text[ln_num]
			ln_num += 1
		num_seq += 1
		seq_len = len(sequence)

		# best score position
		best_score = 0
		for i in range(seq_len - mlen + 1):
			score = 0
			for j in range(mlen):
				score += pwm[letters.index(sequence[i+j])][j]
			if score > best_score:
				best_score = score

		# compute the position p value
		pval = M[mlen][x[mlen].index(best_score)]
		seq_pval = sequence_pval(pval, seq_len)
		p_fast = p_fast * seq_pval	#p_fast is the product of p values for num_seq sequences

	return qfast(num_seq, p_fast) 

def sequence_pval(pval, seq_len):
	tmp = pval * seq_len
	if tmp < 0.001:		#good approximation if value is small
		return tmp
	else:	
		return 1.0 - pow(1.0 - pval, seq_len)
	
def qfast(num_seq, p_fast):
	i = 1
	if (num_seq == 0):
		return 1.0
	if (p_fast == 0):
		return 0.0
	mlnp = -log(p_fast)
	q_fast = t = p_fast
	for i in range(num_seq - 1):
		t = t * mlnp/(i + 1)
		q_fast += t
	return q_fast

def x_compute(i, xvalue, pwm, x):
	if xvalue not in x[i]:
		x[i].append(xvalue)	

		if (i == len(pwm[0])):
			return
		
		for j in range(len(pwm)):
			x_compute(i+1, xvalue+pwm[j][i], pwm, x)

def m_compute(mlen, x, pwm, p=[.25,.25,.25,.25, 0]):
	M = {0:[1.0]}
	for k in range(mlen):
		xvalues = []
		mlist = []
		for item in x[k+1]:
			mvalue = 0
			for j in range(len(pwm)):
				if item-pwm[j][k] in x[k]:
					mvalue += M[k][x[k].index(item-pwm[j][k])] * p[j]
			mlist.append(mvalue)
		M[k+1] = mlist
	return M

def searchFasta(fName, resultsDir, motif):
    #print "~",motifnum,"#"
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

    instances = []
    for key in sInstances:
        for instance in sInstances[key]:
            instances.append(instance)
    
    return makePWM(mLen, [x[0] for x in instances])
    
# constructs a PWM from the given motif instances
def makePWM(mLen, instances, backgroundRatios=[.25,.25,.25,.25, 1]):
    #pdb.set_trace()
    logLikelihood = 0
    #print "~" + str(instances) + "!"
    letters = {'A':0, 'C':1, 'G':2, 'T':3, 'N':4}
    pwm = [[0.0 for x in xrange(mLen)] for y in xrange(len(letters))]
    inc = 1
    
    for instance in instances:
        #fout.write(str(instance)+"\n")
        for pos, char in enumerate(instance):
            if char in letters:
                pwm[letters[char]][pos] += inc
    return pwm



# statistics.py
# author: Jake Martinez
# description: contains functions that give various statistical information for
#    analysis of motifs found within a genetic sequence

from math import *
from compareTool import *
from scipy import stats,interpolate
import statistics
from multiprocessing import Pool

def mean(array):
	return float(sum(array)) / len(array)

def variance(array):
	m = mean(array)

	diffsq = 0.0
	for x in array:
		diffsq += (x - m) ** 2

	return diffsq / len(array)

def stdev(array):
	return variance(array) ** .5


# must be ungapped and same length
def edit_dist(seq1,seq2):
	dist = 0
	for i in range(len(seq1)):
		if seq1[i] != seq2[i]:
			dist += 1
	return dist

# p-Value
#
# parameters:
#     motif=a string for the motif to be evaluated
#     sequences=a string for the filename of a .fasta file with sequences
#	  order=markov model order to be used
#     dist=distribution to use for pvalue (binomial or hypergeometric)
#     pvals=use already generated observed pvalue set
#     threshold=maximum amount of error
def pval(motif,sequences,order=1,dist="binomial",method="qfast",pvals=None,threshold=None):
	mlen = len(motif)

	# default threshold to 20% of the motif length
	if threshold is None:
		threshold = int(len(motif) * .20)

	# generate pvalues if not given
	if pvals is None:
		pvals = pval_all(mlen, sequences, order=order, dist=dist)

	similar = [pvals[x] for x in pvals.keys() if edit_dist(x,motif) <= threshold]

	if len(similar) == 0:
		return 1

	if method == "fisher":
		return stats.combine_pvalues(similar, method="fisher")
	elif method == "stouffer":
		return stats.combine_pvalues(similar, method="stouffer")
	else:
		p_fast = 1
		for p in similar:
			p_fast = p_fast * p
		return qfast(len(similar), p_fast) #p_fast is the product of p values for similar sequences

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


def pval_est_helper(x,motif,threshold=2):
	return edit_dist(x,motif) <= threshold


	

# gets a better estimate of the p-Value by combining similar motifs into a single p-Value
def pval_est(pvals,top=None,threshold=3):
	p_est = dict()

	if top is not None:
		pvs = {x:pvals[x] for x in sorted(pvals,key=pvals.get,reverse=False)[:top]}
	else:
		pvs = pvals

	# pool = Pool(processes=4)
	for motif in pvs:
		# res = [(pool.apply_async(pval_est_helper, (x,motif)),pvals[x]) for x in pvals]
		# similar = [y for x,y in res if x.get()]

		similar = [pvals[x] for x in pvals if edit_dist(x,motif) <= threshold]

		p_fast = 1.0
		for p in similar:
			p_fast = p_fast * p
		qfst = qfast(len(similar), p_fast)
		fsher = stats.combine_pvalues(similar, method="fisher")[1]
		stffr = stats.combine_pvalues(similar, method="stouffer")[1]
		p_est[motif] = (qfst,fsher,stffr)

	return p_est

# p-Values (observed from sample)
#
# desc: 
# parameters:
#     motif=a string for the motif to be evaluated
#     sequences=a string for the filename of a .fasta file with sequences
def pval_all(mlen,sequences,order=1,dist="binomial"):
	threshold = 0
	seq_name = ""

	with open(sequences, "r") as f:
		text = f.read()
		text = text.replace(' ','')
		text = text.split("\n")
	file_length = len(text)
	ln_num = -1

	num_seq = 0
	num_spots = 0
	seq_matched = {}

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
		
		# search sequence
		pos = 0
		while pos < seq_len - mlen:
			num_spots += 1
			match = sequence[pos:pos+mlen]
			if match not in seq_matched:
				seq_matched[match] = 1
			elif match in seq_matched:
				seq_matched[match] += 1
			pos += 1

	matches = sorted(seq_matched, key=seq_matched.get, reverse=True)

	pvals = {}

	n = len(seq_matched)

	mkv_res = pval_markov('A'*mlen,sequences,order)
	mkv = mkv_res[1]
	bg_model = mkv_res[2]

	for match in matches:
		# print '   '+str(matches[i])+' ('+str(seq_matched[matches[i]])+')'
		k = seq_matched[match]
		# print('pval: k='+str(k)+'\tn='+str(n)+'\tP='+str(p))

		p_hat = markov2pval(match,mkv,bg_model,order)
		pvals[match] = stats.binom_test(k,n,p_hat,alternative='greater')


	return pvals

def permutations(perm):
	chars = ['A','G','C','T','N']
	res = []
	for c in chars:
		res += [x+c for x in perm]
	return res


def build_markov(order=1):
	chars = ['A','G','C','T','N']

	if order == 1:
		model = {x:{x:0 for x in chars} for x in chars}
		return model

	perm = chars
	for i in range(1,order):
		perm = permutations(perm)

	return {x:{x:0 for x in chars} for x in perm}


def reverse_compl(text):
	text_rcomp = ""
	for c in text:
		if c == 'A':
			text_rcomp = 'T' + text_rcomp
		elif c == 'T':
			text_rcomp = 'A' + text_rcomp
		elif c == 'C':
			text_rcomp = 'G' + text_rcomp
		elif c == 'G':
			text_rcomp = 'C' + text_rcomp
		elif c == 'N':
			text_rcomp = 'N' + text_rcomp
	return text_rcomp


def pval_markov(motif,sequences,order=1,two_stranded=False):
	model = build_markov(order)
	order_prob = {key:0 for key in model}

	threshold = 0
	seq_name = ""

	with open(sequences, "r") as f:
		text = f.read().split("\n")

	if two_stranded:
		rc = reverse_compl(text)
		file_length = len(rc)
		ln_num = -1
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
				sequence += line
				ln_num += 1
			
			seq_len = len(sequence)
			
			# search sequence and update markov model
			pos = 0
			while pos < seq_len - (order+1):
				match = sequence[pos:pos+order]
				order_prob[match] += 1
				model[match][sequence[pos+order+1]] += 1
				pos += 1

	file_length = len(text)
	ln_num = -1
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
			sequence += line
			ln_num += 1
		
		seq_len = len(sequence)
		
		# search sequence and update markov model
		pos = 0
		while pos < seq_len - (order+1):
			match = sequence[pos:pos+order]
			order_prob[match] += 1
			model[match][sequence[pos+order+1]] += 1
			pos += 1

	p_hat = 1

	# determine probability of initial chars in sequence with length equal to
	# the markov model order
	order_freq = order_prob[motif[0:order]]
	p_hat = order_freq / float(sum([order_prob[key] for key in order_prob]))

	final_markov = model
	for key in final_markov:
		subtotal = sum([final_markov[key][x] for x in final_markov[key]])
		for c in final_markov[key]:
			c_freq = final_markov[key][c]
			if subtotal <= 0:
				final_markov[key][c] = 0
				if c_freq > 0:
					print('statistics.py: a possible error may have occurred (line 274)')
			else:
				final_markov[key][c] = float(c_freq) / float(subtotal)

	# compute probabilities of remaining characters based on markov model
	for i in range(len(motif)-order):
		mkv_key = motif[i:i+order]
		next_char = motif[i+order]

		p_hat *= final_markov[mkv_key][next_char]

		# mkv_prob = model[motif[i:i+order]]
		# freq = mkv_prob[next_char]
		# prob = freq / float(sum([mkv_prob[key] for key in mkv_prob]))
		# p_hat *= prob

	return (p_hat,final_markov,order_prob) # p_hat, markov, bg model

# isn't completely accurate for order>1, doesn't account for first few characters
def markov2pval(motif,markov,bg_model,order=1):
	p_hat = 1

	order_freq = bg_model[motif[0:order]]
	p_hat = order_freq / float(sum([bg_model[key] for key in bg_model]))

	for i in range(len(motif)-order):
		mkv_key = motif[i:i+order]
		next_char = motif[i+order]

		p_hat *= markov[mkv_key][next_char]

	return p_hat

# z-score
#
# desc: statistic describing the overrepresentation of a motif; uses scipy
#		zscore function
# parameters: 
#     motif=motif string
#     sequences=sequences filename
def zscore(motif,sequences):
	threshold = 0
	seq_name = ""
	mlen = len(motif)

	with open(sequences, "r") as f:
		text = f.read().split("\n")
	file_length = len(text)
	ln_num = -1

	num_seq = 0
	num_spots = 0
	seq_matched = {}
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
		
		# search sequence
		pos = 0
		while pos < seq_len - mlen:
			num_spots += 1
			match = sequence[pos:pos+mlen]
			if match not in seq_matched:
				seq_matched[match] = 1
			elif match in seq_matched:
				seq_matched[match] += 1
			pos += 1

	matches = sorted(seq_matched, key=seq_matched.get, reverse=True)

	zscores = stats.zscore([seq_matched[x] for x in matches])

	if motif in matches:
		return zscores[matches.index(motif)]
	else:
		return 0


# z-score
#
# desc: statistic describing the overrepresentation of a motif; uses scipy
#		zscore function
# parameters: 
#     motif=motif string
#     sequences=sequences filename
def zscore_all(mlen,sequences):
	threshold = 0
	seq_name = ""

	with open(sequences, "r") as f:
		text = f.read().split("\n")
	file_length = len(text)
	ln_num = -1

	num_seq = 0
	num_spots = 0
	seq_matched = {}
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
		
		# search sequence
		pos = 0
		while pos < seq_len - mlen:
			num_spots += 1
			match = sequence[pos:pos+mlen]
			if match not in seq_matched:
				seq_matched[match] = 1
			elif match in seq_matched:
				seq_matched[match] += 1
			pos += 1

	matches = sorted(seq_matched, key=seq_matched.get, reverse=True)

	zscores = stats.zscore([seq_matched[x] for x in matches])

	results = {matches[i]:zscores[i] for i in range(len(zscores))}

	return results

# pi
#
# desc: used for q-value computation
# lambda is a threshold, usually 0.5 works
def pi(data, mlen, lmda = 0.5,pvals=None):
	if pvals is None:
		pvals = [v for k,v in pval_all(mlen,data).iteritems() if v <= lmda]
	return len(pvals) / (len(data)*(1-lmda))

# q-Value (estimated)
#
# desc: statistic describing the false discovery rate
# parameters: 
#     mlen=motif length
#     sequences=fasta filename
#     pvals=output from pval_all should be given here for some speedup
def qval_all(mlen,sequences,pvals=None):
	if pvals is None:
		pv = pval_all(mlen,sequences)
		pvals = sorted([pv[x] for x in pv])
	else:
		pv = pvals
		pvals = sorted([pv[x] for x in pv])

	lmbds = [x*.01 for x in range(0,95,1)]
	pis = [pi(sequences,mlen,lmda=lmbd,pvals=pvals) for lmbd in lmbds]

	tck = interpolate.splrep(lmbds,pis,k=3) # cubic spline with 3 deg freedom
	x = [1]

	pi_0 = interpolate.splev(x,tck)[0] # evaluate spline at f(1)
	m = len(pvals)

	qvals = [0 for x in pvals] # initialize to size of pval array
	qvals[m-1] = pi_0*m*pvals[m-1]/(m-1)

	i = m - 2
	while i >= 0:
		p_i = pvals[i]
		qvals[i] = min( (pi_0*m*p_i)/(i+1) , qvals[i+1] )
		i -= 1

	#build mapping again
	motifs = sorted(pv, key=pv.get)
	qv = {motifs[i]:qvals[i] for i in range(len(qvals))}

	return qv

# should not really be used because it has to redo a lot of computations
def qval(motif,sequences,qvals=None,pvals=None):
	if pvals is None:
		pvals = pval_all(len(motif),sequences)
		pvals = sorted([pvals[x] for x in pvals])
	if qvals is None:
		qvals = qval_all(len(motif),sequences,pvals=pvals)

	# pv = pval(motif,sequences)
	# ind = 0
	# while pvals[ind] < pv:
	# 	ind += 1
	return qvals[motif]



###### PLACEHOLDER FUNCTIONS ######
def loglik(motif,sequences):
	pass

def seqspec(motif,sequences):
	pass

def comb(n,r):
    f = factorial
    return f(n) / f(r) / f(n-r)

def normal_pdf(x, m, v):
	return 1.0/sqrt(2*pi*v) * exp(-(x-m)**2/(2*v))

def binomial_pdf(k, n, p):
	if n < 100:
		return comb(n, k) * p**k * p**(n-k)
	return normal_pdf(k, n*p, n*p*(1.0-p))



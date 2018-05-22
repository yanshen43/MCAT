#! /usr/bin/python
# randomized motif-embedding sequence generator
# by: jake martinez (jrm98)


import random
import argparse
import time
from pprint import pprint


def main():
	opt = argparse.ArgumentParser()
	opt.add_argument('-w', default=8, type=int, help='motif width')
	opt.add_argument('-N', default=20, type=int, help='number of sequences')
	opt.add_argument('-n', default=100, type=int, help='sequence length')
	opt.add_argument('-e', default=0, type=int, help='max # of SNPs in motif')
	opt.add_argument('--eprob', default=.65, type=float, help='probability of a SNP existing in the motif')
	opt.add_argument('--maxloc', default=-1, type=int, help='max starting location for motifs (-1 represents the end of the sequence)')
	opt.add_argument('--minloc', default=0, type=int, help='min starting location for motifs')
	opt.add_argument('-P', default=1, type=float, help='probability motif is present in sequence')
	opt.add_argument('-o', default=None, type=str, help='output filename')
	opt.add_argument('--dir', default='', type=str, help='output directory')
	opt.add_argument('--dyad', help='spaced dyad mode (overrides nmotifs parameter, -w option becomes the total width of the combined dyad)', action='store_true')
	opt.add_argument('--nmotifs', default=1, type=int, help='max number of distinct motifs in a single sequence')
	opt.add_argument('--nsites', default=1, type=int, help='max number of appearances of the motif in a single sequence (not implemented)')
	opt.add_argument('--negative', help='used to generate a sequence for a negative discriminative set', action='store_true')
	opt.add_argument('--template', default=None, type=str, help='template sequence file path')
	opt.add_argument('--bgprob', default=None, type=str, help='background probability file path')
	opt.add_argument('--markov', default=None, type=str, help='markov probability file')
	opt.add_argument('-v', help='turn on output messages', action='store_true')
	opt.add_argument('--verbose', help='turn on output messages', action='store_true')
	args=opt.parse_args()

	N = args.N
	n = args.n
	e = args.e
	ep = args.eprob
	P = args.P
	w = args.w
	minloc = args.minloc
	maxloc = args.maxloc
	output = args.o
	directory = args.dir
	negative = args.negative
	nmotifs = args.nmotifs
	dyad = args.dyad
	tmp = args.template
	bg = args.bgprob
	mkv = args.markov
	verbose = args.verbose or args.v

	# create a default filename if none was given
	if output is None:
		output = 'seq_w'+str(w)+'_N'+str(N)+'_n'+str(n)
		if e != 0:
			output += '_e'+str(e)
		if ep != .65:
			output += '_ep'+str(int(ep*100))
		if P != 1:
			output += '_P'+str(int(ep*100))
		output += '-'+str(int(time.time()))+'.fasta'

	if directory is None:
		directory = ''


	if tmp is not None:
		print('seq_gen: loading from template...')
		template = load_template(tmp,verbose=verbose)
		if verbose:
			print('seq_gen: loaded statistical model: ')
			pprint(template)
		generate(w,N,n,e,ep,P,minloc,maxloc,output,directory,negative,
			nmotifs,dyad,markov=template,verbose=verbose)

		#show the Markov model generated from the synthetic sequence
		if verbose:
			gen_template = load_template(directory+output,save=False,verbose=verbose)
			print('seq_gen: statistical model of output: ')
			pprint(gen_template)

			#get background probs from template?
			# gen_bgprob = load_bgprob(directory+output)
			# for key in gen_bgprob:
			# 	print('  '+key+': '+str(bgprob[key])+' => '+str(gen_bgprob[key]))



	elif mkv is not None:
		print('seq_gen: loading from Markov file...')
		markov = load_markov(mkv)
		if verbose:
			print('seq_gen: loaded statistical model: ')
			pprint(markov)
		generate(w,N,n,e,ep,P,minloc,maxloc,output,directory,negative,
			nmotifs,dyad,markov=markov,verbose=verbose)

		#show the Markov model generated from the synthetic sequence
		if verbose:
			gen_template = load_template(directory+output,save=False,verbose=verbose)
			print('seq_gen: statistical model of output: ')
			pprint(gen_template)

	elif bg is not None:
		print('seq_gen: loading from bgprob file...')
		bgprob = load_bgprob(bg)
		if verbose:
			print('seq_gen: loaded statistical model: ')
			pprint(bgprob)
		generate(w,N,n,e,ep,P,minloc,maxloc,output,directory,negative,
			nmotifs,dyad,bgprob=bgprob,verbose=verbose)

		#show the background probability model generated from the synthetic sequence
		if verbose:
			gen_bgprob = load_bgprob(directory+output)
			print('seq_gen: statistical model of output: ')
			for key in gen_bgprob:
				print('  '+key+': '+str(bgprob[key])+' => '+str(gen_bgprob[key]))

	else:
		generate(w,N,n,e,ep,P,minloc,maxloc,output,directory,negative,
		nmotifs,dyad,verbose=verbose)


def generate(w,N,n,e=0,ep=.65,P=1,minloc=0,maxloc=-1,
	output='seq.fasta',directory='',negative=False,
	nmotifs=1,dyad=False,bgprob={'A':.25,'T':.25,'C':.25,'G':.25},
	markov=None,verbose=False):

	# look for potential errors in input parameters
	psum = bgprob['A'] + bgprob['T'] + bgprob['C'] + bgprob['G']
	if psum > 1.0 or psum < .99:
		print('seq_gen: invalid background probabilities')
		return

	if len(directory) > 0 and directory[-1:] != '/':
		directory += '/'

	# reflect input parameters in output
	if verbose:
		if negative and not dyad:
			print('seq_gen: NEGATIVE=true')
		elif dyad:
			print('seq_gen: DYAD=true')
			print('seq_gen: w='+str(w)+'\tN='+str(N)+'\tn='+str(n)+'\te='+str(e)+'\tP='+str(P))
		else:
			print('seq_gen: w='+str(w)+'\tN='+str(N)+'\tn='+str(n)+'\te='+str(e)+'\tP='+str(P))
		print('\toutput => '+directory+output)

	sequences = []
	motif = []

	if dyad:
		w1 = w - random.randint(2,w-2) #requires total motif width of at least 4
		w2 = w - w1
		motif.append(gen_seq(w1))
		motif.append(gen_seq(w2))
		if verbose:
			print('seq_gen: motifs generated...')
		for i in range(N):
			sequence = gen_seq(n,bgprob=bgprob,template=markov)
			result = embed_motif(sequence, motif, maxerror=e, errorprob=ep, P=P, minstart=minloc, maxstart=maxloc, nmotifs=2)
			sequences.append((result,motif,i))
			# print result, motif, i
	else:
		for i in range(nmotifs):
			motif.append(gen_seq(w))

		if verbose:
			print('seq_gen: motifs generated...')

		for i in range(N):
			sequence = gen_seq(n,bgprob=bgprob,template=markov)
			if negative:
				result = embed_motif(sequence, motif, maxerror=0, errorprob=0, P=0, minstart=minloc, maxstart=maxloc)
			else:
				result = embed_motif(sequence, motif, maxerror=e, errorprob=ep, P=P, minstart=minloc, maxstart=maxloc, nmotifs=nmotifs)
			sequences.append((result,motif,i))
			# print result, motif, i
	if verbose:
			print('seq_gen: sequences generated...')
	write_to_fasta(directory+output, sequences, negative)
	if verbose:
			print('seq_gen: wrote sequences to file...')
	return (sequences,motif)

# generates a string from the alphabet 'ACGT' of the given length iteratively
def gen_seq(length,bgprob={'A':.25,'T':.25,'C':.25,'G':.25},template=None):
	if template is None:
		seq = ''
		for i in range(length):
			vals = [float(x) for x in bgprob.values()]
			m = sum(vals)
			r = random.uniform(0,m)
			total = 0.0
			for c in bgprob:
				if total + float(bgprob[c]) >= r:
					seq += c
					break
				total += float(bgprob[c])
		return seq
	else:
		seq = ''
		for i in range(length):
			if i == 0:
				vals = [float(x) for x in bgprob.values()]
				m = sum(vals)
				r = random.uniform(0,m)
				total = 0.0
				for c in bgprob:
					if total + float(bgprob[c]) >= r:
						seq += c
						break
					total += float(bgprob[c])
			else:
				m = sum([float(x) for x in template[seq[-1]].values()])
				r = random.uniform(0,m)
				total = 0.0
				for c in template[seq[-1]]:
					if total + float(template[seq[-1]][c]) >= r:
						seq += c
						break
					total += float(template[seq[-1]][c])
		return seq
	#return ''.join(random.SystemRandom().choice('ACGT') for _ in range(length))

# embeds motif in the given sequence
def embed_motif(seq, motif, maxerror=0, errorprob=.20, minstart=0, maxstart=-1, P=1,nmotifs=1):
	if random.random() < P and (maxstart > minstart or maxstart == -1):
		if maxstart < 0:
			maxstart = len(seq) - len(motif[0]) - 1

		# pick starting location for motif
		loc = [random.randint(minstart,maxstart)]

		# introduce artificial error in motif
		errors = [0]
		errorind = [[]]
		r = random.random()
		while r < errorprob and errors[0] < maxerror:
			# choose location for snp, if already snp pick another loc
			snp = random.randint(0,len(motif[0])-1)
			while snp in errorind[0]:
				snp = random.randint(0,len(motif[0])-1)

			# remove the correct bp from the options for replacement
			letters = 'ACGT'.replace(motif[0][snp],'')

			# replace bp for another random bp
			motif[0] = motif[0][:snp] \
				+ random.SystemRandom().choice(letters) \
				+ motif[0][snp+1:]
			errors[0] += 1
			errorind[0].append(snp)
			r = random.random()

		# embed modified motif in sequence
		newseq = seq[:loc[0]] + motif[0] + seq[loc[0] + len(motif[0]):]
		
		# sort snp locations in motif before returning
		errorind[0].sort()

		# try to add more motifs
		if nmotifs > 1 and random.random() < P:
			res = embed_motif(newseq, motif[1:],maxerror,errorprob,loc[0]+len(motif[0]),maxstart,P,nmotifs-1)
			newseq = res[0]
			loc += res[1]
			errors += res[3]
			errorind += res[4]
		return (newseq, loc, motif, errors, errorind)
	else:
		return (seq, [-1], None, [None], [None])

# writes a set of sequences to a fasta file
def write_to_fasta(filename, sequences, negative):
	with open(filename, 'w') as f:
		for sequence in sequences:
			if negative:
				f.write('>artificial'+str(sequence[2])+'-'+str(int(time.time()))+' NEGATIVE sequence'+'\n')
			else:
				f.write('>artificial'+str(sequence[2])+'-'+str(int(time.time()))+' seq w/ embedded motif:'+str(sequence[1])+' start_loc:'+str(sequence[0][1])+' errors:'+str(sequence[0][3])+' errorind:'+str(sequence[0][4])+'\n')
			linelen = 80
			line = sequence[0][0]
			for i in range(0, len(line), linelen):
				f.write(line[i:i+linelen]+'\n')
			f.write('\n\n')
	pass

def load_template(filename, order=1, save=True, verbose=False):
	model = build_markov(order)
	order_prob = {key:0 for key in model}

	threshold = 0
	seq_name = ""

	with open(filename, "r") as f:
		text = f.read().split("\n")

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

	final_markov = model
	for key in final_markov:
		subtotal = sum([final_markov[key][x] for x in final_markov[key]])
		for c in final_markov[key]:
			c_freq = final_markov[key][c]
			final_markov[key][c] = float(c_freq) / float(subtotal)

	if verbose:
		print('seq_gen: loaded template.')

	if save:
		fname = 'mkv_'+str(int(time.time()*100))+'.txt'
		save_markov(final_markov, fname)
		if verbose:
			print('seq_gen: saved template to '+fname)

	return final_markov

def save_markov(markov, filename):
	with open(filename, "w") as f:
		for k in markov:
			for c in markov[k]:
				f.write("P(%s|%s) = %f\n" % (c, k, markov[k][c]))



# loads markov model from file, wont work if file doesn't have all possibilities
def load_markov(filename):
	markov = dict()
	with open(filename, "r") as f:
		text = f.read().replace(" ","")
		text = text.replace("P(","")
		text = text.replace("|",",")
		text = text.replace(")=",",") # each line should now be c,k,v
		lines = text.split()

		for line in lines:
			print(line)
			if line == '':
				continue
			ckv = line.split(",")
			c = ckv[0]
			k = ckv[1]
			v = ckv[2]
			if k not in markov:
				markov[k] = dict()
			markov[k][c] = v
	return markov


# bgprob file format:
#    A=.25
#    T=.25
#    C=.25
#    G=.25
def load_bgprob(filename):
	with open(filename, 'r') as f:
		text = f.read().split('\n')

	prob = {'A':.25,'T':.25,'C':.25,'G':.25}
	for line in text:
		if '=' not in line:
			continue
		line = line.replace(' ','')
		if len(line) < 1:
			continue
		if line[0] == 'A':
			prob['A'] = float(line.split('=')[1])
		elif line[0] == 'T':
			prob['T'] = float(line.split('=')[1])
		elif line[0] == 'C':
			prob['C'] = float(line.split('=')[1])
		elif line[0] == 'G':
			prob['G'] = float(line.split('=')[1])
		else:
			print('seq_gen: unexpected line format in bgprob file')

	return prob



def build_markov(order=1):
	chars = ['A','G','C','T']

	if order == 1:
		model = {x:{x:0 for x in chars} for x in chars}
		return model

	perm = chars
	for i in range(1,order):
		perm = permutations(perm)

	return {x:{x:0 for x in chars} for x in perm}

def permutations(perm):
	chars = ['A','G','C','T']
	res = []
	for c in chars:
		res += [x+c for x in perm]
	return res

if __name__ == '__main__':
	main()

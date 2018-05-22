#! /usr/bin/python
# batch tester for synthetic data
# by: jake martinez (jrm98)

from seq_gen import generate
import subprocess

def test(name,w,N,n,e=0,ep=.65,P=1,minloc=0,maxloc=-1,directory='data/',negative=False,
	nmotifs=1,dyad=False):

	#auto names batch test files with prefix 'syn_batch_'
	(sequences,motifs) =generate(w,N,n,e=e,ep=ep,P=P,minloc=minloc,
		maxloc=maxloc, output='syn_batch_'+name+'.fasta',
		directory=directory,negative=negative,nmotifs=nmotifs,dyad=dyad)

	#auto generate complementary negative sequence
	# generate(w,N,n,output='syn_batch_'+name+'_neg.fasta',
	# 	directory=directory,negative=True)

	cmdLine = "python orange_pipeline.py data/syn_batch_"+name+".fasta "+"data/n_genes.fasta"
	subprocess.call([cmdLine], shell=True)

    #see if we found a motif that was generated
	res = parse_results()

	fout = open('results/syn_test.txt','a')
	fout.write('test'+name+'\t=> results: '+str(res)+'\n')
	fout.write('\t=> actual: '+str(motifs)+'\n\n')

	for found in res:
		if found in motifs:
			return True
	return False



def parse_results():
	try:
		with open('results/results.html','r') as f:
			lines = f.readlines()
			results = []
			for line in lines:
				if 'Motif:' in line:
					m = line.replace('Motif:','').replace(' ','').replace('\n','')
			 		results += [m]
			return results
	except Exception:
		print 'file not found...'
	return results

def editdist(string1,string2):
	pass

def cleanup():
	subprocess.call(["./cleanup_syn_batch.sh"], shell=True)
	pass

def main():
	print 'beginning synthetic batch test...'

	#test results
	res = []

	#parameters for synthetic data
	testname = ''
	numtests = 10
	width = 8
	numsequences = 100
	seqlength = 500
	maxerrors = 5
	errorprob = .55
	motifprob = 1
	nummotifs = 1
	dyadmode = False


	#varying motif width tests
	for x in range(numtests):
		res.append(test(testname+str(x),width,numsequences,seqlength,
			e=maxerrors,ep=errorprob,P=motifprob,nmotifs=nummotifs,
			dyad=dyadmode))

	width = 12
	for x in range(numtests):
		res.append(test(testname+str(x),width,numsequences,seqlength,
			e=maxerrors,ep=errorprob,P=motifprob,nmotifs=nummotifs,
			dyad=dyadmode))

	width = 16
	for x in range(numtests):
		res.append(test(testname+str(x),width,numsequences,seqlength,
			e=maxerrors,ep=errorprob,P=motifprob,nmotifs=nummotifs,
			dyad=dyadmode))

	width = 12
	motifprob = .65
	for x in range(numtests):
		res.append(test(testname+str(x),width,numsequences,seqlength,
			e=maxerrors,ep=errorprob,P=motifprob,nmotifs=nummotifs,
			dyad=dyadmode))

	width = 12
	motifprob = .45
	for x in range(numtests):
		res.append(test(testname+str(x),width,numsequences,seqlength,
			e=maxerrors,ep=errorprob,P=motifprob,nmotifs=nummotifs,
			dyad=dyadmode))
	
	width = 12
	motifprob = .25
	for x in range(numtests):
		res.append(test(testname+str(x),width,numsequences,seqlength,
			e=maxerrors,ep=errorprob,P=motifprob,nmotifs=nummotifs,
			dyad=dyadmode))

	print(res)
	pass

if __name__ == '__main__':
	main()
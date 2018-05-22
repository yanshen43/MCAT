#! /usr/bin/python
# batch tester for synthetic data
# by: jake martinez (jrm98)

from seq_gen import generate
from statfile import StatFile
import subprocess, os, time, sys
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

def accuracy(found,motifs):
	if found in motifs:
		return 1.0
	# for m in motifs:

	# possible options:
	#    - 2/3/4-mer counts
	#    - 1-(edit distance/motif length)
	return 0.0

def test(name,w=12,N=15,n=500,e=0,ep=.65,P=1,minloc=0,maxloc=-1,
	directory='test/data/',negative=False,nmotifs=1,dyad=False,threshold=2, log=None):

	seqfilename = 'testdata_'+name+'.fasta'
	negfilename = 'testdata_'+name+'_neg.fasta'

	if log is not None:
		log.write('[%s]: Using sequences file "%s"\n' % (time.strftime('%H:%M:%S'),seqfilename) )
		log.flush()

	#generate sequences file
	(sequences,motifs) = generate(w,N,n,e=e,ep=ep,P=P,minloc=minloc,
		maxloc=maxloc, output=seqfilename,
		directory=directory,negative=negative,nmotifs=nmotifs,dyad=dyad)

	#generate complementary negative sequences file
	(nsequences,_) = generate(w,N,n,output=negfilename,
		directory=directory,negative=True)

	#run pipeline for test data
	cmdLine = "python orange_pipeline.py "+directory+seqfilename+" "+directory+negfilename+" -s test/stat > test/logs/output.log 2> test/logs/output_error.log"

	if log is not None:
		log.write('[%s]: Beginning command "%s"\n' % (time.strftime('%H:%M:%S'),cmdLine) )
		log.flush()

	subprocess.call([cmdLine], shell=True)

	if log is not None:
		log.write('[%s]: Looking through statdir\n' % time.strftime('%H:%M:%S'))
		log.flush()

	statdir='test/stat' # path to stat directory
	statfiles = sorted([ f for f in os.listdir(statdir) if f.startswith('stat_')])

	if len(statfiles) <= 0:
		return None

	sf_filename = str(statfiles[0])

	sf_id = sf_filename.replace('stat','').replace('.','').replace('_','')

	if log is not None:
		log.write('[%s]: Found StatFile "%s"\n' % (time.strftime('%H:%M:%S'),sf_filename) )
		log.flush()

	sf = StatFile('test/stat/'+sf_filename)

	return sf

	# judge how well motif was found based on top three results
	# for i in range(1,4):
	# 	found = sf.get_by('rank',i)
	# 	acc = accuracy(found['motif'],motifs)
	# 	if acc > 0.50:
	# 		return acc / (2**(i-1)) # adjust accuracy score based on rank

	# return 0.0

def graph_res(res,name):
	points = []
	for (x,sf) in res:
		sf.load()
		ys = sf.sort_by('pval')[:20] # get 20 rows with best pvals

		for y in ys:
			points.append( (x,y) )

	fig = plt.figure(figsize=(7,7))
	ax = fig.add_subplot(111)

	ax.scatter(points)
	ax.set_xlabel(name)
	ax.set_ylabel('p-Value')

	fig.savefig('test/img/%s.png'%name)


def test_width(motifprob, guessed_width, interval=2, 
	start=6, stop=20, numtests=10, log=None):
	res = []

	for width in range(start, stop, interval):
		for x in range(numtests):
			res.append( (width,test('width'+str(width)+'-'+str(x),w=width,P=motifprob, log=log)) )

	graph_res(res,'width')

	return res

def test_mprob(guessed_width, interval=1, granularity=0.1, 
	start=0.1, count=10, numtests=10, log=None):
	res = []
	vals = [start+(x*granularity) for x in range(count)]

	for motifprob in vals:
		for x in range(numtests):
			res.append( (motifprob,test('mprob'+str(motifprob)+'-'++str(x),w=guessed_width,P=motifprob, log=log)) )

	graph_res(res,'mprob')

	return res

def test_eprob(max_e, guessed_width, interval=1, granularity=0.1, 
	start=0.1, count=10, numtests=10, log=None):
	res = []
	vals = [start+(x*granularity) for x in range(count)]

	for eprob in vals:
		for x in range(numtests):
			res.append( (eprob,test('eprob'+str(eprob)+'-'++str(x),w=guessed_width,e=max_e,ep=eprob, log=log)) )
	
	graph_res(res,'eprob')

	return res

def test_mloc(mloc, guessed_width, motifprob=.9, interval=2, 
	start=6, stop=20, numtests=10, log=None):
	res = []

	for mloc in range(start, stop, interval):
		for x in range(numtests):
			res.append( (mloc,test('mloc'+str(mloc)+'-'++str(x),w=guessed_width,P=motifprob, log=log)) )
	
	graph_res(res,'mloc')

	return res

def test_nmotifs():
	pass

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
	subprocess.call(["./test/cleanup_syn_batch.sh"], shell=True)
	pass

def main():
	print 'mcat_test: beginning synthetic batch test...'

	#test results
	res = []

	#parameters for synthetic data
	testname = ''
	numtests = 10
	w = 8
	N = 100
	n = 500
	e = 5
	eprob = .55
	P = 1
	nmotifs = 1
	dyadmode = False

	directory = 'test/data'

	#check for batch test data directory, create if doesn't exist
	if not os.path.exists(os.path.join(os.getcwd(),directory)):
		sys.stdout.write('mcat_test: creating test data directory...\n')
		os.makedirs(directory)
	
	if not os.path.exists(os.path.join(os.getcwd(),'test/img')):
		sys.stdout.write('mcat_test: creating test img directory...\n')
		os.makedirs('test/img')

	if not os.path.exists(os.path.join(os.getcwd(),'test/logs')):
		sys.stdout.write('mcat_test: creating test logs directory...\n')
		os.makedirs('test/logs')

	if not os.path.exists(os.path.join(os.getcwd(),'test/stat')):
		sys.stdout.write('mcat_test: creating test logs directory...\n')
		os.makedirs('test/stat')

	with open('test/logs/test.log','w') as log:
		#motif width sensitivity
		sys.stdout.write('mcat_test: running motif width sensitivity tests...')
		sys.stdout.flush()
		begin = time.time()
		res += test_width(.9, 12, log=log)
		end = time.time()
		sys.stdout.write('Done (%fs).\n' % (end-begin))
		sys.stdout.flush()

		#motif probability sensitivity
		sys.stdout.write('mcat_test: running motif probability sensitivity tests...')
		sys.stdout.flush()
		begin = time.time()
		res += test_mprob(12, log=log)
		end = time.time()
		sys.stdout.write('Done (%fs).\n' % (end-begin))
		sys.stdout.flush()

		#motif error sensitivity
		sys.stdout.write('mcat_test: running motif error sensitivity tests...')
		sys.stdout.flush()
		begin = time.time()
		res += test_eprob(5, 12, log=log)
		end = time.time()
		sys.stdout.write('Done (%fs).\n' % (end-begin))
		sys.stdout.flush()

		#motif location sensitivity
		# sys.stdout.write('mcat_test: running motif location sensitivity tests...')
		# sys.stdout.flush()
		# begin = time.time()
		# res += test_mloc(12)
		# end = time.time()
		# sys.stdout.write('Done (%fs).\n' % (end-begin))
		# sys.stdout.flush()

	pass

if __name__ == '__main__':
	try:
		main()
	except Exception as err:
		with open('log.txt','w') as f:
			f.write(str(err))

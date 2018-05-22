#! /usr/bin/python
# batch tester for real data
# by: jake martinez (jrm98)

import subprocess

def test(posseq, negseq='data/n_genes.fasta'):
	cmdLine = "python orange_pipeline.py "+posseq+" "+negseq
	subprocess.call([cmdLine], shell=True)

	f = open('results/results.html','r')
	html = f.read()
	f.close()
	copy = open('results/res/'+get_filename(posseq)+'.html','w')
	copy.write(html)
	copy.close()

	#see if we found a motif that was generated
	return parse_results()



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

def get_filename(s):
    index = s.rfind('/')
    end = s.rfind('.')
    if index < 0:
        return s
    else:
        return s[index+1:end]

def main():
	print('beginning test...')

	#test results
	res1 = []
	res2 = []

	#parameters for synthetic data
	species = ['dm','hs','pneumocystis','scerevisiae','scryophilus',
		'sjaponicus','soctosporus','sp']

	genes = ['mad1','mad2','mad3','bub1','bub3','mph1']

	for s in species:
		res1 += [test('data/'+s+'.fasta')]
	print('by species: '+str(res1))

	for g in genes:
		res2 += [test('data/'+g+'.fasta')]
	print('by gene: '+str(res2))
	
	pass

if __name__ == '__main__':
	main()
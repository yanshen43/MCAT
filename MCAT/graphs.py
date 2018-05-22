# an interface for storing and accessing statistical data for motifs
# by: jake martinez (jrm98)

from statfile import StatFile
import statistics
import argparse
from scipy import stats
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

def pval_vs_atcg_ratio(sf):
	ratio = dict()
	for motif in pvals.keys():
		count = 0
		total = 0
		for c in motif:
			total += 1
			if c in ['A','T']:
				count += 1

		if total > 0:
			ratio[motif] = float( float(count)/float(total) )
		else:
			ratio[motif] = 0

	print ratio
	print pvals

	sf.add_vals('ratio',ratio)

	df = sf.copy()

	print df

	fig = plt.figure(figsize=(7,7))
	ax = fig.add_subplot(111)

	ax.scatter(df['pval'],df['ratio'])
	ax.set_xlabel('p-Value')
	ax.set_ylabel('A/T to C/G Ratio')

	fig.savefig('stat/pval_vs_atcg_ratio.png')

	pass

def pval_vs_at_ratio(sf):
	ratio = dict()
	for motif in pvals.keys():
		count = 0
		total = 0
		for c in motif:
			if c in ['A','T']:
				total += 1
				if c == 'A':
					count += 1

		if total > 0:
			ratio[motif] = float( float(count)/float(total) )
		else:
			ratio[motif] = 0

	print ratio
	print pvals

	sf.add_vals('ratio',ratio)

	df = sf.copy()

	print df

	fig = plt.figure(figsize=(7,7))
	ax = fig.add_subplot(111)

	ax.scatter(df['pval'],df['ratio'])
	ax.set_xlabel('p-Value')
	ax.set_ylabel('A to T Ratio')

	fig.savefig('stat/pval_vs_at_ratio.png')

	pass

def pval_vs_cg_ratio(sf):
	ratio = dict()
	for motif in pvals.keys():
		count = 0
		total = 0
		for c in motif:
			if c in ['C','G']:
				total += 1
				if c == 'C':
					count += 1

		if total > 0:
			ratio[motif] = float( float(count)/float(total) )
		else:
			ratio[motif] = 0

	print ratio
	print pvals

	sf.add_vals('ratio',ratio)

	df = sf.copy()

	print df

	fig = plt.figure(figsize=(7,7))
	ax = fig.add_subplot(111)

	ax.scatter(df['pval'],df['ratio'])
	ax.set_xlabel('p-Value')
	ax.set_ylabel('C to G Ratio')

	fig.savefig('stat/pval_vs_cg_ratio.png')


def compare_cols(sf,col1,col2):
	df = sf.copy()

	fig = plt.figure(figsize=(7,7))
	ax = fig.add_subplot(111)

	ax.scatter(df[col1],df[col2])
	ax.set_xlabel(col1)
	ax.set_ylabel(col2)

	# compute spearman rank correlation
	print('graphs.py: computing spearman rank correlation...')
	corr = stats.spearmanr(df[col1],df[col2])
	print('   '+str(corr))

	filename = 'stat/%s_vs_%s.png'%(col1,col2)
	print('graphs.py: generating graph and writing to "'+str(filename)+'"')
	fig.savefig(filename)

def main():
	opt = argparse.ArgumentParser()
	opt.add_argument('statfile', default='', help='motif width')
	args = opt.parse_args()

	sf = StatFile(args.statfile)
	sf.load()
	
	# print('graphs.py: comparing pval to [A,T]:[C,G] ratio...')
	# pval_vs_atcg_ratio(sf)

	# print('graphs.py: comparing pval to A:T ratio...')
	# pval_vs_at_ratio(sf)

	# print('graphs.py: comparing pval to C:G ratio...')
	# pval_vs_cg_ratio(sf)

	print('graphs.py: comparing pval to qval')
	compare_cols(sf,'pval','qval')

	print('graphs.py: comparing pval to zscore')
	compare_cols(sf,'pval','zscore')

	# print('graphs.py: comparing pval to rank')
	# compare_cols(sf,'pval','rank')

if __name__ == '__main__':
	main()


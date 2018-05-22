# an interface for storing and accessing statistical data for motifs
# by: jake martinez (jrm98)

import pandas as pd
import numpy as np

class StatFile(object):
	"""interface for storing statistical data for motifs"""
	def __init__(self, filename):
		self._df = pd.DataFrame(columns=['rank','pval','zscore','qval'])
		self.filename = filename

	def __str__(self):
		return str(self._df)

	def __getitem__(self):
		return self._df

	def save(self):
		try:
			self._df.to_csv(self.filename)
			return True
		except Exception as err:
			print('statfile.py: %s' % err)
			return False

	def load(self):
		try:
			self._df = pd.read_csv(self.filename, index_col=0)
			return True
		except Exception as err:
			print('statfile.py: %s' % err)
			return False

	def copy(self):
		return self._df.copy()

	def add_vals(self,colname,vals):
		if colname not in self._df.columns:
			self._df[colname] = np.nan

		for k in vals:
			if k in self._df.index:
				self._df[colname][k] = vals[k]
			else:
				self.add_data(k,colname,vals[k])
		pass

	def get_motif(self,motif):
		return self._df.loc[motif]

	def get_all(self):
		return self._df

	def get_by(self,col,val):
		return self._df.loc[self._df[col] == val]

	def sort_by(self,col,ascending=True):
		return self._df.sort_values(col,axis=0,ascending=ascending)

	def add_data(self,motif,colname,value):
		if colname not in self._df.columns:
			self._df[colname] = np.nan

		if motif in self._df.index:
			self._df[colname][motif] = value
		else:
			row = {colname:value}
			if colname not in ['rank','pval','zscore','qval']:
				new_df = pd.DataFrame([row], index=[motif], columns=['rank','pval','zscore','qval',colname])
			else:
				new_df = pd.DataFrame([row], index=[motif], columns=['rank','pval','zscore','qval'])
			self._df = self._df.append(new_df)

	def add_motif(self,motif,rank=None,pval=None,zscore=None,qval=None):
		row = {'rank':rank,'pval':pval,'zscore':zscore,'qval':qval}
		if motif in self._df.index:
			print('statfile.py: motif already exists in DataFrame')
			return
		new_df = pd.DataFrame([row], index=[motif], columns=['rank','pval','zscore','qval'])
		self._df = self._df.append(new_df)

	def update_motif(self,motif,rank=None,pval=None,zscore=None,qval=None):
		if motif in self._df.index:
			if rank is not None:
				self._df['rank'][motif] = rank
			if pval is not None:
				self._df['pval'][motif] = pval
			if zscore is not None:
				self._df['zscore'][motif] = zscore
			if qval is not None:
				self._df['qval'][motif] = qval
		else:
			print('statfile.py: motif does not exist in DataFrame')

if __name__ == '__main__':
	import statistics
	sf = StatFile('stat/stat.csv')

	if sf.load():
		print('statfile: loaded file successfully')

	pvals = statistics.pval_all(8,'data/test.fasta')
	sf.add_vals('pval',pvals)

	pvals2 = statistics.pval_all(8,'data/test.fasta',order=2)
	sf.add_vals('pval_2',pvals2)

	qvals = statistics.qval_all(8,'data/test.fasta',pvals=pvals)
	sf.add_vals('qval',qvals)

	zscores = statistics.zscore_all(8,'data/test.fasta')
	sf.add_vals('zscore',zscores)

	print sf.sort_by('pval').iloc[0]
	print sf.sort_by('zscore',ascending=False).iloc[0]

	if sf.save():
		print('statfile: saved file successfully')

#! /user/bin/python
# tests the results from MCAT, XXMotifs..
#by: christy coghlan (ccogs)

import os
import numpy as np
from itertools import islice

#from parseTestResults.py
def match(s1, s2):
    return max([sum([1 if c1 == c2 else 0 for c1, c2 in zip(s1, s2[p:]+s2[:p])])
        for p in xrange(len(s2))])

#from syn_test.py
#filename is .html results file for mcat
def get_mcat_results(filename):  
	try:
		with open(filename,'r') as f:
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


#get results from xxmotif run
def get_xxmotif_results(filename):
    abs_filename = os.path.abspath(filename)
    try:
        with open(abs_filename, 'r') as f:
            lines = f.readlines()
            results = []
            for line in lines:
                if '=>' in line:
                    m = line.replace('=>', '').replace(' ','').replace('\n', '')
                    results += [m]
            return results
    except Exception as err:
        print err
    return results

#get results from xxmotif pwm data
def get_xxmotif_results_pwm(filename):
	abs_filename = os.path.abspath(filename)
	print filename
 	try:
        	with open(abs_filename, 'r') as f:
            		motifs = []

			while True:
        			next_n_lines = list(islice(f, 6))
        			if not next_n_lines:
            				break
		#		print next_n_lines
        	    		motifs += [get_motif_from_pwm(next_n_lines)]
		#		print motifs
            		return motifs
    	except Exception as err:
        	print err
    	#return []

#given a pwm, return the motif
#using a simple max algorithm
def get_motif_from_pwm(lines):
	arr = []
	for idx, line in enumerate(lines[1:5]):
		line = line.rstrip()
		#print line, idx
        	split = line.split("\t")
		rowdata = map(float, split)
		arr += [rowdata]

	#get max columns
   	#print arr
	#print ["{0:0.2f}".format(i) for i in arr]
   	arr_numpy = np.array(arr)
   	maximums_columns = np.argmax(arr_numpy, axis=0)
   
	#print maximums_columns	
	ans = ""
	for c in maximums_columns:
		if c == 0:
			ans += 'A'
		elif c == 1:
			ans += 'C'
		elif c == 2:
			ans += 'G'
		elif c == 3:
			ans += 'T'
		else:
			ans += "INVALID"

	#print ans
	return ans
	#ans = [tup[1] for tup in maximums_columns]
	#return ans

#score results of a test
def score_results(motifs, found, tool_name):
    count = 0
    found_list = []
    not_found_list = []
    for m in motifs:
	for res in found:
		if match(res, m.upper()) >= len(m)*.75:
        	    	count = 1
            		found_list.append(m)
        	else:
            		not_found_list.append(m)

    out = ""
    out = 'tool ' + tool_name + ' found motifs: ' + str(found) + '\n'
    out += 'actual was: ' + str(motifs)
    out += '\ncorrect: ' + str(count) + ' of ' + str(len(motifs))
    print out
    return (count, len(motifs))
            

#gets motifs from a synthetic file
def get_motifs_from_syn_gen(filename):
	#use this if you're using syn_gen to generate	
	abs_filename = os.path.abspath(filename)

	try:
		with open(abs_filename,'r') as f:
			lines = f.readlines()
			motifs = []
			first_line = lines[0]
			motifs_part = first_line.split('[')[1]
			
			motifs_part = motifs_part.split(']')[0]
			motifs_part = motifs_part.replace("'", "")
			
		
			motifs  = motifs_part.split(', ')	
			return motifs
	except Exception as err:
		print err
	return motifs
    
#parses out the name of the motif from the filename
def get_motifs_from_filename(filename):
	motif = filename.replace('.fasta', '').replace(' ', '').replace('.', '')
	return [motif]

def main():
    print 'beginning scoring of different pipelines'

    input_directory = 'test_data/BuddingYeastData2/pos/'
    mcat_results_directory = ''
    mcat_suffix = '.mcat.out'
    xxmotif_results_directory = 'test_results/xmotif_yeast_motiffiles_2/'
    xxmotif_suffix = '_pos_hf_MotifFile.txt'
    xxmotif_pwm_directory = 'test_results/xmotif_yeast_results_with_neg2/home/ccogs/gasp/test_data/BuddingYeastData2/pos/'
    xxmotif_pwm_suffix = '_pos_hf.pwm'
    xxmotif_correct = 0
    total_motifs = 0

    #iterate over input files
    for filename in os.listdir(input_directory):
	motifs = get_motifs_from_filename( filename)
        #mcat_results = get_mcat_results(mcat_results_directory + '/' + filename + mcat_suffix)

	#print motifs
	passing_filename = filename.replace('.fasta', '')
        #xxmotif_results = get_xxmotif_results(xxmotif_results_directory + passing_filename + xxmotif_suffix)
	xxmotif_results = get_xxmotif_results_pwm(xxmotif_pwm_directory + passing_filename + "/" + passing_filename + xxmotif_pwm_suffix)
	#print "got motifs: " + str(motifs)

	(xx_correct, total) = score_results(motifs, xxmotif_results, 'xxmotif')
        #score_results(motifs, mcat_results, 'mcat')                                        
	
	xxmotif_correct += xx_correct
	total_motifs += total

    print "total correct for xxmotif = " + str(xxmotif_correct) + " of " + str(total_motifs)
    
     


if __name__ == '__main__':
	main()

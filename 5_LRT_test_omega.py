#!/usr/bin/env python

import os
import sys
import re
import argparse
import statsmodels.api as sm
import pandas as pd
import scipy
from statsmodels.sandbox.stats.multicomp import multipletests

def lrt_test(input_dir, null_ext, alt_ext, df, alpha):
	null2lnl = dict()
	alternative2lnl = dict()
	for file in os.listdir(input_dir):
		if file.endswith(alt_ext):
			file_path = os.path.join(input_dir, file)
			null_file = file_path.replace(alt_ext, null_ext)
			core_gene = file.replace(alt_ext, "")

			null_dict = dict()
			with open(null_file) as null_file:
				for line in null_file:
					line = line.rstrip()
					line = re.sub("[\t ]+", "\t", line)
					tabs = line.split("\t")

					if "lnL" in line:
						lnl_value = float(tabs [4])
						null2lnl[core_gene] = lnl_value
					else:
						pass
			null_file.close()	

			with open(file_path) as file_path:
				for line in file_path:
					line = line.rstrip()
					line = re.sub("[\t ]+", "\t", line)
					tabs = line.split("\t")

					if "lnL" in line:
						lnl_value = float(tabs [4])
						alternative2lnl[core_gene] = lnl_value
					else:
						pass

			file_path.close()
			
	###Doing test
	all_p_vals = list()
	all_core_genes = list()
	all_lr_stats = list()
	for key, value in alternative2lnl.items():
		lnl_alt = value
		lnl_null = null2lnl[key]
		LR_statistic = -2*(lnl_null-lnl_alt)
		p_val = scipy.stats.chi2.sf(LR_statistic, df)
		all_p_vals.append(p_val)
		all_core_genes.append(key)
		all_lr_stats.append(LR_statistic)

	##Correcting p-vals
	p_adjusted = multipletests(all_p_vals, alpha=alpha, method='bonferroni')
	p_adjusted_vals = p_adjusted[1]

	print("Core gene" + "\t" + "LR_stat"+ "\t" + "P value" + "\t" + "Corrected P value" + "\t" + "Result")
	count = 0
	for x in p_adjusted_vals:
		core = all_core_genes[count]
		p_val = all_p_vals[count]
		lr_stat = all_lr_stats[count]
		if x<alpha:
			print(core + "\t" + str(lr_stat) +"\t" + str(p_val) +"\t" + str(x) + "\t" + "*")
		else:
			print(core + "\t" + str(lr_stat) +"\t" + str(p_val) +"\t" + str(x) + "\t" + "NS")
		count = count + 1


def main(argv=None):
	args_parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description="INFO:\nThis script will perform a Likelihood-ratio test, which will assesses the goodness of fit of two competing statistical models.\nP-values will be correctec through the bonferroni method.", epilog='*******************************************************************\n\n*******************************************************************\n\nMake sure you cite our book chapter!')
	args_parser.add_argument('-i', '--input', required=True, help='Input folder where codeml output files are located.')
	args_parser.add_argument('-n', '--null', required=True, help='Extension of codeml output files that resulted from model = 0 and NSsites = 0 (null hypothesis).')
	args_parser.add_argument('-a', '--alternative', required=True, help='Extension of codeml output files that resulted from model = 1 and NSsites = 0 in our example (alternative hypothesis). Moel and NSsites may change depending on your question.')
	args_parser.add_argument('-df', '--df', required=True, help='Degrees of freedom for LR test.')
	args_parser.add_argument('-alpha', '--alpha', required=False, default=int(0.05), help='Alpha value for P-values correction.')
	args_parser = args_parser.parse_args()

	#Setting up parameters
	input_dir = args_parser.input
	null = args_parser.null
	alternative = args_parser.alternative
	df = int(args_parser.df)
	alpha = float(args_parser.alpha)
	
	lrt_test(input_dir, null, alternative, df, alpha)

if __name__ == '__main__':
	status = main()

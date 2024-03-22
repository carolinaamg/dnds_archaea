#!/usr/bin/env python

import os, sys, re, subprocess, shlex, argparse

def predict_proteins(input_dir):
	outfolder = input_dir
	for i in os.listdir(input_dir):
		if i.endswith(".fna"): 
			fasta = os.path.join(input_dir, i)
			core = re.sub(".fna", "", i)

			prot = os.path.join(outfolder, core+".faa")
			nucl = os.path.join(outfolder, core+".genes.fna")

			cmd = "prodigal -i "+ fasta +" -a "+ prot +" -d "+ nucl
			print (i)
			cmd2 = shlex.split(cmd)
			subprocess.call(cmd2, stdout=open("prodigal.out", "w"), stderr=open("prodigal.err", "w"))

	print("DONE.")

def main(argv=None):
	args_parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description="INFO:\nThis script will run Prodigal on a set of genome files. Make sure Prodigal is installed.\n\nIMPORTANT:\nYour input files must have an .fna extension.", epilog='*******************************************************************\n\n*******************************************************************\n\nMake sure you cite Prodigal and our book chapter!')
	args_parser.add_argument('-i', '--input', required=True, help='Input folder where genome files are located with ".fna" extension.')
	args_parser = args_parser.parse_args()

	#Setting up parameters
	input_dir = args_parser.input
	
	predict_proteins(input_dir)

if __name__ == '__main__':
	status = main()

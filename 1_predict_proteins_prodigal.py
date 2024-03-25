#!/usr/bin/env python

import os, sys, re, subprocess, shlex, argparse

def predict_proteins(input_dir,outfolder, quiet):
	for i in os.listdir(input_dir):
		if i.endswith(".fna"): 
			fasta = os.path.join(input_dir, i)
			core = re.sub(".fna", "", i)

			prot = os.path.join(outfolder, core+".faa")
			nucl = os.path.join(outfolder, core+".genes.fna")

			if quiet == 1:
				cmd = "prodigal -i "+ fasta +" -a "+ prot +" -d "+ nucl
			else:
				cmd = "prodigal -i "+ fasta +" -a "+ prot +" -d "+ nucl + " -q"
			print (cmd)
			cmd2 = shlex.split(cmd)
			subprocess.call(cmd2, stdout=open("prodigal_log.out", "w"), stderr=open("prodigal.err", "w"))

	print("DONE.")

def main(argv=None):
	args_parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description="INFO:\nThis script will run Prodigal on a set of genome files. Make sure Prodigal is installed.\n\nIMPORTANT:\nYour input files must have an .fna extension.", epilog='*******************************************************************\n\n*******************************************************************\n\nMake sure you cite Prodigal and our book chapter!')
	args_parser.add_argument('-i', '--input', required=True, help='Input folder where genome files are located with ".fna" extension.')
	args_parser.add_argument('-o', '--output', required=True, help='Output folder where your output files will go.We suggest you use the same folder as your input to keep your files more organized.')
	args_parser.add_argument('-q', '--quiet', required=False, default=int(1), help='Run Prodigal quietly? yes? use 0, no? use 1. Default is 1.')
	args_parser = args_parser.parse_args()

	#Setting up parameters
	input_dir = args_parser.input
	outfolder = args_parser.output
	quiet = int(args_parser.quiet)
	
	predict_proteins(input_dir,outfolder, quiet)

if __name__ == '__main__':
	status = main()

#!/usr/bin/env python

import os, sys, re, argparse
from Bio.Phylo.PAML import codeml

def run_codeml(input_dir, tree, extension, model, nssites, clean):
	for file in os.listdir(input_dir):
		if file.endswith(".pal2nal"):
			alignment = os.path.join(input_dir,file)
			codeml_output = re.sub(".pal2nal", extension, alignment)
			core_gene = file.replace(".pal2nal","")
			print("Your input files are: ", alignment, "and", tree)
			
			#Let's run codeml!
			print("Running codeml for core gene:", core_gene)
			cml = codeml.Codeml()
			cml.alignment = alignment
			cml.tree = tree
			cml.out_file = codeml_output
			cml.working_dir = input_dir

	 		#Setting options
			cml.set_options(noisy=1) #Limited information on screen 
			cml.set_options(verbose=1) #Details of output file
			cml.set_options(seqtype=1) #Analysis on codons
			cml.set_options(ndata = 1) #Number of alignments
			cml.set_options(icode = 0) #Universal genetic code
			cml.set_options(cleandata=clean) #Remove positions with ambiguous data (might be bad if too many positions are)		
			cml.set_options(runmode=0) #Analysis will use tree topology 
			cml.set_options(model = model) #Compute an omega value for each branch		
			cml.set_options(NSsites = nssites)#One omega for all sites
			cml.set_options(CodonFreq=1) #Mutation-selection model with observed codon frequencies used as estimates. This model accounts for the mutational bias and selection affecting codon usage,
			cml.set_options(clock=0) #Clock assumption is not used (unrooted phylogeny must be used)
			cml.set_options(fix_omega=0)#Estimate omega (not used for seqtype 1)
			cml.set_options(omega=0.5)#Starting omega value (not used for seqtype 1)

			cml.print_options()#Print the options on the screen


			result = cml.run(verbose=True)
			print("Omega for core gene", core_gene, "succesfully calculated.")

	print("Run finished!")

def main(argv=None):
	args_parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description="INFO: \n script to estimate dN/dS a.k.a omega values on codon alignemnts.", epilog='*******************************************************************\n\n*******************************************************************\n\nMake sure you cite PAML and our book chapter!')
	args_parser.add_argument('-i', '--input', required=True, help='Input folder of codon alignments (ending in .pal2nal, obtained using PAL2NAL).')
	args_parser.add_argument('-t', '--tree', required=True, help='Input phylogenetic tree.')
	args_parser.add_argument('-e', '--extension', required=True, help='Extension of output file (e.g., .cml.out).')
	args_parser.add_argument('-m', '--model', required=False, default=int(0), help='Model for omega estimation. Default value is 0, meaning that all branches have the same rate.')
	args_parser.add_argument('-n', '--nssites', required=False, default=int(0), help='NSsites option available to estimate omega. Default value is 0, meaning that codons in the alignment have the same rate.')
	args_parser.add_argument('-c', '--clean', required=False, default=int(0), help='Do you expect your data to have a lot of gaps and/or ambiguous sites? yes? use -c 1, no? use default value -c 0.')
	args_parser = args_parser.parse_args()

	#Setting up parameters to run omega
	input_dir = args_parser.input
	tree = args_parser.tree
	extension = str(args_parser.extension)
	model = int(args_parser.model)
	nssites = args_parser.nssites
	clean = args_parser.clean
	
	run_codeml(input_dir, tree, extension, model, nssites, clean)

if __name__ == '__main__':
	status = main()

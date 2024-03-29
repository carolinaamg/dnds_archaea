#!/usr/bin/env python

import os, sys, subprocess, shlex, argparse
from collections import defaultdict
from Bio import SeqIO
from Bio import Phylo

def gene_trees(core_genome, genome_folder, ingroup_branch, quiet):
	###Getting a nucleotide file for each core gene
	aminoacid_seqs = defaultdict(list)
	protein2genome = dict()
	for file in os.listdir(core_genome):
		if file.endswith(".faa"):
			file_path = os.path.join(core_genome, file)
			core_gene = file.replace(".faa","")

			for record in SeqIO.parse(file_path, "fasta"):
				seq_name = str(record.id).split("&")
				protein_name = seq_name[1]
				aminoacid_seqs[core_gene].append(protein_name)
				genome_name = seq_name[0].replace(".faa","")
				protein2genome[protein_name] = genome_name

	###Finding nuc sequences that are part of the core genome
	all_genomes_seqs = dict()
	for file in os.listdir(genome_folder):
		if file.endswith(".genes.fna"):
			file_path = os.path.join(genome_folder, file)
			
			for record in SeqIO.parse(file_path, "fasta"):
				seq_name = str(record.id)
				all_genomes_seqs[seq_name] = record

	###Linking core gene and nucleotide sequences
	for key, value in aminoacid_seqs.items():
		core_seqs = list()
		output = os.path.join(core_genome, key + ".fna")
		for i in value:
			genome_name = protein2genome[i]
			i_seq = all_genomes_seqs[i]
			i_seq.id = genome_name
			i_seq.description = ""
			core_seqs.append(i_seq)
		SeqIO.write(core_seqs, output, "fasta")
		print("Nucleotide sequences for core gene", key, "are ready!")

	###Getting nucleotide alignments
	for file in os.listdir(core_genome):
		if file.endswith(".fna"):
			file_path = os.path.join(core_genome, file)
			output = file_path.replace(".fna",".fna.aln")
			
			if quiet == 1:
				cmd = "mafft --auto "+ file_path 
			else:
				cmd = "mafft --auto --quiet "+ file_path 
			cmd2 = shlex.split(cmd)
			print("Running alignment: ",cmd)
			subprocess.call(cmd2, stdout=open(output, "w"), stderr=open("mafft_log.txt", "w"))

	###Getting gene trees from our freshly made nucleotide alignments :D
	for file in os.listdir(core_genome):
		if file.endswith(".fna.aln"):
			family_name = file.replace(".fna.aln","")
			input_aln = os.path.join(core_genome, file)
			output_tree = input_aln.replace(".fna.aln", ".fasttree")

			cmd = "fasttree -nt "+ input_aln 
			cmd2 = shlex.split(cmd)
			print("Running your tree with FastTree --> ",cmd)
			subprocess.call(cmd2, stdout=open(output_tree, "w"), stderr=open("fasttree_log.txt", "a"))

			print("Modifying your tree to have it ready for CODEML...")

			treefile = Phylo.read(output_tree, "newick")
			output_tree_mod = input_aln.replace(".fna.aln", ".temporary")
			output_tree_nwk = open(input_aln.replace(".fna.aln", ".nwk"), "w")

			total_tips = 0
			for tip in treefile.get_terminals():
				total_tips = total_tips + 1
				if str(tip.name) == ingroup_branch:
					tip.name = ingroup_branch + " #1"
				else:
					pass

			Phylo.write(treefile, output_tree_mod, "newick")

			output_tree_nwk.write(str(total_tips)+"  1\n")
			with open(output_tree_mod) as temporary:
				for line in temporary:
					line = line.rstrip()
					line = line.replace("'","")
					output_tree_nwk.write(line +"\n")		
			temporary.close()
			output_tree_nwk.close()

			os.remove(output_tree_mod)
			print("Your output tree for CODEML is", family_name+ ".nwk")


def main(argv=None):
	args_parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description="INFO:\nThis script will build a gene tree for each core gene. These gene trees will be used by CODEML to estimate distance between genomes. Make sure you have MAFFT and fasttree installed!", epilog='*******************************************************************\n\n*******************************************************************\n\nMake sure you cite MAFFT, Fasttree, and our book chapter!')
	args_parser.add_argument('-c', '--core', required=True, help='Input folder where core genome files are located/Corecruncher output files.')
	args_parser.add_argument('-g', '--genomes', required=True, help='Folder where genome files (.fna) and Prodigal output files are located.')
	args_parser.add_argument('-ib', '--ingroup_branch', required=True, help='What is your ingroup branch? Make sure you write the name of the branch exactly as it appears in your tree, including special characters.')	
	args_parser.add_argument('-q', '--quiet', required=False, default=int(1), help='Run MAFFT quietly? yes? use 0, no? use 1. Default is 1.')
	args_parser = args_parser.parse_args()

	#Setting up parameters
	core_genome = args_parser.core
	genome_folder = args_parser.genomes
	ingroup_branch = str(args_parser.ingroup_branch)
	quiet = int(args_parser.quiet)
	
	gene_trees(core_genome, genome_folder, ingroup_branch, quiet)

if __name__ == '__main__':
	status = main()

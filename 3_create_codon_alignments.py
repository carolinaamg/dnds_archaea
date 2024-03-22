#!/usr/bin/env python

import os, sys, subprocess, shlex, argparse
from collections import defaultdict
from Bio import SeqIO

def codon_alignment(core_genome, genome_folder, quiet):
	core2genome = dict()
	core2seqs = defaultdict(list)
	###Getting amino acid alignments
	for file in os.listdir(core_genome):
		if file.endswith(".faa"):
			core_gene = file.replace(".faa","")
			file_path = os.path.join(core_genome, file)
			output = file_path.replace(".faa",".aln.prot")

			for record in SeqIO.parse(file_path, "fasta"):
				seq_name = str(record.id).split("&")
				protein_name = seq_name[1]
				genome_name = seq_name[0].replace(".faa", "")
				core2genome[protein_name] = genome_name
				core2seqs[core_gene].append(protein_name)
				
			if quiet == 1:
				cmd = "mafft --auto "+ file_path 
			else:
				cmd = "mafft --auto --quiet "+ file_path 
			cmd2 = shlex.split(cmd)
			print("Running: ",cmd)
			subprocess.call(cmd2, stdout=open(output, "a"), stderr=open("mafft_log.txt", "a"))

	for file in os.listdir(core_genome):
		if file.endswith(".aln.prot"):
			file_path = os.path.join(core_genome, file)
			output = file_path.replace(".aln.prot",".aln.seq")

			all_seqs = list()
			for record in SeqIO.parse(file_path, "fasta"):
				name = str(record.id).split("&")
				genome_name = name [0].replace(".faa","")
				record.id = genome_name
				record.description = ""
				all_seqs.append(record)
			SeqIO.write(all_seqs, output, "fasta")
			os.remove(file_path)

	###Finding nuc sequences that are part of the core genome
	all_genomes_seqs = dict()
	for file in os.listdir(genome_folder):
		if file.endswith(".genes.fna"):
			file_path = os.path.join(genome_folder, file)
			
			for record in SeqIO.parse(file_path, "fasta"):
				seq_name = str(record.id)
				all_genomes_seqs[seq_name] = record

	###Getting nucleotide sequences for each core gene
	for key, value in core2seqs.items():
		all_seqs = list()
		output = os.path.join(core_genome, key + ".genes.fna")
		for i in value:
			seq_record = all_genomes_seqs[i]
			seq_record.description = ""
			seq_record.id = str(core2genome[i])
			all_seqs.append(seq_record)

		SeqIO.write(all_seqs, output, "fasta")

	###Getting codon alignments with PAL2NAL
	for file in os.listdir(core_genome):
		if file.endswith(".aln.seq"):
			input_alignment = os.path.join(core_genome, file)
			genes_seq = input_alignment.replace(".aln.seq",".genes.fna")
			output_pal2nal = input_alignment.replace(".aln.seq",".pal2nal")

			cmd = "perl pal2nal.pl " + input_alignment + " " + genes_seq + " -output paml -nogap"
			print ("Running: ",cmd)
			cmd2 = shlex.split(cmd)
			subprocess.call(cmd2, stdout = open(output_pal2nal, "w"), stderr = open("pal2nal_log.txt", "w"))

def main(argv=None):
	args_parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description="INFO:\nThis script will build individual codon core-genome alignments with the extension .pal2nal. Make sure you have MAFFT install, and edit Line 72 of this script to add the absolute path of the pal2nal binary.", epilog='*******************************************************************\n\n*******************************************************************\n\nMake sure you cite MAFFT, PAL2NAL, and our book chapter!')
	args_parser.add_argument('-c', '--core', required=True, help='Input folder where core genome files are located/Corecruncher output files. Your output files will be located here!')
	args_parser.add_argument('-g', '--genomes', required=True, help='Folder where genome files (.fna) and Prodigal output files are located.')
	args_parser.add_argument('-q', '--quiet', required=False, default=int(1), help='Run MAFFT quietly? yes? use 0, no? use 1. Default is 1.')
	args_parser = args_parser.parse_args()

	#Setting up parameters
	core_genome = args_parser.core
	genome_folder = args_parser.genomes
	quiet = int(args_parser.quiet)
	
	codon_alignment(core_genome, genome_folder, quiet)

if __name__ == '__main__':
	status = main()

#!/usr/bin/env python

import os, sys, subprocess, shlex, argparse
from collections import defaultdict
from Bio import SeqIO

def run_program(core_genome, genome_folder, output_alignment):
	###Getting a nucleotide file for each core gene
	aminoacid_seqs = defaultdict(list)
	genome_list = list()
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
				if genome_name not in genome_list:
					genome_list.append(genome_name)
				else:
					pass

	###Finding nuc sequences that are part of the core genome
	all_genomes_seqs = dict()
	for file in os.listdir(genome_folder):
		if file.endswith(".genes.fna"):
			file_path = os.path.join(genome_folder, file)
			
			for record in SeqIO.parse(file_path, "fasta"):
				seq_name = str(record.id)
				all_genomes_seqs[seq_name] = record

	###Linking core genes and nucleotide sequences
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
		print("Nucleotide sequences for core gene", key, "is ready!")

	###Getting nucleotide alignments
	for file in os.listdir(core_genome):
		if file.endswith(".fna"):
			file_path = os.path.join(core_genome, file)
			output = file_path.replace(".fna",".aln")
			cmd = "/overflow/bobaylab/carolina/mafft-7.505-with-extensions/scripts/mafft --auto "+ file_path 
			cmd2 = shlex.split(cmd)
			print("Running: ",cmd)
			subprocess.call(cmd2, stdout=open(output, "a"), stderr=open("log_file.txt", "a"))

	###Concatenating nucleotide sequences
	genome2seqs = defaultdict(list)
	for file in os.listdir(core_genome):
		if file.endswith(".aln"):
			core_seqs = dict()
			path_file = os.path.join(core_genome, file)

			for record in SeqIO.parse(path_file, "fasta"):
				genome_name = str(record.id)
				core_seqs[genome_name] = record
				sequence = list(record.seq)
				seq_len = len(sequence)


			for genome in genome_list:
				if genome in core_seqs:
					core_record = core_seqs[genome]
					genome2seqs[genome].append(str(core_record.seq))
				else:
					core_record = "-"*seq_len
					genome2seqs[genome].append(core_record)

	output_alignment = open(output_alignment, "w")
	###Getting concatenated alignment
	for key, value in genome2seqs.items():
		all_seqs = "".join(value)
		output_alignment.write(">"+ key + "\n")
		output_alignment.write(all_seqs + "\n")
	output_alignment.close()

	print("Concatenated alignment", output, "is ready!")

def main(argv=None):
	args_parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description="INFO:\nThis script will build a nucleotide core-genome concatenated alignment for phylogenetic reconstruction.", epilog='*******************************************************************\n\n*******************************************************************\n\nMake sure you cite MAFFT and our book chapter!')
	args_parser.add_argument('-c', '--core', required=True, help='Input folder where core genome files are located/Corecruncher output files.')
	args_parser.add_argument('-g', '--genomes', required=True, help='Folder where genome files (.fna) and Prodigal output files are located.')
	args_parser.add_argument('-o', '--output', required=True, help='Your output concatenated alignment.')
	args_parser = args_parser.parse_args()

	#Setting up parameters
	core_genome = args_parser.core
	genome_folder = args_parser.genomes
	output = args_parser.output
	
	run_program(core_genome, genome_folder, output)
	return 0

if __name__ == '__main__':
	status = main()
	sys.exit(status)
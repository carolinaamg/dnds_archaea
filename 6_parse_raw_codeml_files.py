#!/usr/bin/env python

import argparse
import os
from collections import defaultdict
import sys


def parse_dnds(input_folder, output, ext):
    # Get list of all infiles.
    files = os.listdir(input_folder)
    files = [f for f in files if f.endswith(ext)]

    # Throw error if 0 files found.
    if len(files) == 0:
        print('No files found with extension:', ext, file=sys.stderr)
        sys.exit(1)

    # Make output directory, unless it exists.
    if not os.path.exists(output):
        os.makedirs(output)

    poly_pos = defaultdict(set)
    summary_tab = defaultdict(dict)
    branch_tab = []

    trees = dict()

    branch_header = None
    just_printed_sites = False
    just_parsed_counts = False
    next_line_is_poly_breakdown = False
    tree_next = False
    poly_breakdown_length = None

    for f in files:

        branch_observed = False

        # Get input name, which is filename without extension.
        input_name = f.replace(ext, '')

        summary_tab[input_name]['num_compare'] = 0
        summary_tab[input_name]['summed_N'] = 0
        summary_tab[input_name]['summed_S'] = 0

        with open(os.path.join(input_folder, f)) as file:
            lines = file.readlines()
            for line in lines:
                line_split = line.rstrip().split()

                if len(line_split) == 0:
                    continue

                if next_line_is_poly_breakdown:
                    if poly_breakdown_length is not None and len(line_split) == poly_breakdown_length:
                        all_pos = ''.join(line_split[1:])
                        seq_poly_pos = set([i for i, x in enumerate(all_pos) if x != '.'])
                        poly_pos[input_name] = poly_pos[input_name].union(seq_poly_pos)
                    else:
                        next_line_is_poly_breakdown = False

                if line.startswith('Printing out site pattern counts'):
                    just_printed_sites = True
                elif just_printed_sites:
                    just_printed_sites = False
                    summary_tab[input_name]['num_seq'] = line_split[0]
                    summary_tab[input_name]['num_sites'] = line_split[1]
                    just_parsed_counts = True
                elif just_parsed_counts:
                    next_line_is_poly_breakdown = True
                    just_parsed_counts = False
                    poly_breakdown_length = len(line_split)
                elif line_split[0] == 'kappa':
                    summary_tab[input_name]['kappa'] = line_split[3]
                elif line_split[0] == 'omega':
                    summary_tab[input_name]['omega'] = line_split[3]
                elif line.startswith('tree length for dN:'):
                    summary_tab[input_name]['total_dN'] = line_split[-1]
                elif line.startswith('tree length for dS:'):
                    summary_tab[input_name]['total_dS'] = line_split[-1]
                if line.startswith('tree length = '):
                    summary_tab[input_name]['tree_length'] = line_split[-1]
                    tree_next = True
                elif tree_next:
                    trees[input_name] = line.rstrip()
                    tree_next = False
                elif line_split[0] == 'branch':
                    if branch_header is None:
                        branch_header = '\t'.join(['locus'] + line_split)
                    # Get mapping of which column index corresponds to 'N' and 'S'.
                    N_index = line_split.index('N')
                    S_index = line_split.index('S')
                    branch_num_col = len(line_split)
                    branch_observed = True
                elif branch_observed and len(line_split) == branch_num_col:
                    branch_tab.append('\t'.join([input_name] + line_split))
                    summary_tab[input_name]['num_compare'] += 1.0
                    summary_tab[input_name]['summed_N'] += float(line_split[N_index])
                    summary_tab[input_name]['summed_S'] += float(line_split[S_index])

    # Check that all info parsed per locus.
    loci_w_missing_info = set()
    for locus in summary_tab:
        if 'kappa' not in summary_tab[locus]:
            loci_w_missing_info.add(locus)
        elif 'omega' not in summary_tab[locus]:
            loci_w_missing_info.add(locus)
        elif summary_tab[locus]['num_compare'] == 0:
            loci_w_missing_info.add(locus)

    if len(loci_w_missing_info) > 0:
        print('Loci with missing info: ' + ', '.join(sorted(list(loci_w_missing_info))), file=sys.stderr)

    # Write output.
    with open(os.path.join(output, 'summary_tab.tsv'), 'w') as outfile:
        summary_header = '\t'.join(['locus', 'kappa', 'omega', 'total_dN', 'total_dS', 'mean_N', 'mean_S',
                                    'num_polymorphic_sites'])
        outfile.write(summary_header + '\n')

        for locus in sorted(summary_tab.keys()):

            # Skip if the output files was incomplete.
            if locus in loci_w_missing_info:
                continue

            mean_N = str(summary_tab[locus]['summed_N'] / summary_tab[locus]['num_compare'])
            mean_S = str(summary_tab[locus]['summed_S'] / summary_tab[locus]['num_compare'])

            num_polymorphisms = str(len(poly_pos[locus]))

            summary_line = '\t'.join([locus,
                                      summary_tab[locus]['kappa'],
                                      summary_tab[locus]['omega'],
                                      summary_tab[locus]['total_dN'],
                                      summary_tab[locus]['total_dS'],
                                      mean_N,
                                      mean_S,
                                      num_polymorphisms])
            outfile.write(summary_line + '\n')

    with open(os.path.join(output, 'branch_tab.tsv'), 'w') as branchfile:
        branchfile.write(branch_header + '\n')
        branchfile.write('\n'.join(branch_tab))

    # Make trees output folder, unless it exists.
    if not os.path.exists(os.path.join(output, 'trees')):
        os.makedirs(os.path.join(output, 'trees'))

    for locus in trees.keys():
        with open(os.path.join(output, 'trees', locus + '.tre'), 'w') as treefile:
            treefile.write(trees[locus])


def main(argv=None):
    args_parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="INFO:\nThis script will parse raw CODEML output files for further analysis.",
        epilog='*******************************************************************\n\n*******************************************************************\n\nMake sure you cite our book chapter!')
    args_parser.add_argument('-i', '--input', required=True, help='Input folder where codeml output files are located.')
    args_parser.add_argument('-o', '--output', required=True, help='Folder where output files will go.')
    args_parser.add_argument('-e', '--ext', required=False, default=".cml.out", help='Extension of codeml output files. Default is .cml.out.')
    args_parser = args_parser.parse_args()

    # Setting up parameters.
    input_folder = args_parser.input
    output = args_parser.output
    extension = args_parser.ext

    parse_dnds(input_folder, output, extension)


if __name__ == '__main__':
    status = main()

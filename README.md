## A workflow for estimating the ratio of non-synonymous to synonymous substitutions in core protein-coding genes

### Introduction

This repository provides a step-by-step guide to estimate the dN/dS ratio on protein-coding genes of archeal genomes. This framework is also applicable to bacteria.

Moreover, this repository provides example input files to go from genome files to the final statistical analyses to compare two models.

Input example files (genomes and a formatted tree) can be found in the following Zenodo repository: https://zenodo.org/records/10854407

*This repository is part of the book chapter: "Title", where we provide a detailed summary of the scripts shown above and a step-by-step tutorial to estimate dN/dS with the example genome files.*


### Installation

This repository can be downloaded by running:
```
git clone git@github.com:carolinaamg/dnds_archaea.git
```

The dependencies needed to run this workflow are described in our chapter. To install these dependencies in a _conda_ environment, you can use the following commands.

_Optional_: We recommend using _[mamba](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html)_ below for efficiency when installing packages in _conda_ environments. To do so, you will need to follow the _mamba_ [installation instructions]([https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html]) and then simply swap in `mamba` for all instances of `conda` below.

To create a new conda environment called `dnds_workflow` and to install most dependencies in this environment:
```
conda create --name dnds_workflow \
    anaconda::python \
    anaconda::statsmodels \
    bioconda::mafft \
    bioconda::pal2nal \
    bioconda::paml \
    bioconda::prodigal \
    bioconda::raxml \
    conda-forge::biopython
```

To activate this environment:
```
conda activate dnds_workflow
```

_[CoreCruncher](https://github.com/lbobay/CoreCruncher)_ needs to be installed separately.

Ddownload _CoreCruncher_ and make the scripts are executable. Note taht you need to specify the full path when running CoreCruncher.
```
git clone git@github.com:lbobay/CoreCruncher.git

cd CoreCruncher

chmod 775 *.py
```

Then download [_USEARCH_ v6.1](https://drive5.com/usearch/download.html), which is required by _CoreCruncher_. Note that you can use the free version only for non-commerical use (please see the CoreCruncher documentation on how to alternatively use `BLAST`).

You will need to also add the _USEARCH_ binary to your _$PATH_, and rename the binary (indicated by the example placeholder `usearch6.1.XX` below) to be `usearch61`, which is what _CoreCruncher_ expects.
```
ln -s $PWD/usearch6.1.XX /path/to/envs/dnds_workflow/bin/usearch61
```


### Tutorial

The following commands will illustrate how to run this workflow on the _Methanosarina barkeri_ vs. _Methanosarcina mazei_ comparison described in our chapter. This workflow will be conducted in `tutorial_working`.

#### 0. Download the genome FASTAs and species tree (in newick format), and extract the genomes.

```
mkdir tutorial_working

cd tutorial_working

wget https://zenodo.org/records/10854407/files/Methanosarcina_genus_genomes.tar.gz?download=1 \
	-O Methanosarcina_genus_genomes.tar.gz

wget https://zenodo.org/records/10854407/files/test_tree.nwk?download=1 \
	-O test_tree.nwk

tar -xzvf Methanosarcina_genus_genomes.tar.gz
```

#### 1. Call core genes and build alignment

First, run `1_predict_proteins_prodigal.py` to wrap [Prodigal](https://github.com/hyattpd/Prodigal) to call genes across these genomes:

```
mkdir prodigal_out

1_predict_proteins_prodigal.py \
	-i Methanosarcina_genus_genomes \
	-o prodigal_out
```

Then run `CoreCruncher` to identify core genes:

```
python /path/to/corecruncher_master.py \
	-in prodigal_out \
	-out Methanosarcina_genus_genomes_out \
	-ext .faa \
	-freq 100 \
	-score 80 \
	-length 80 \
```

Note that we specified the reference genome's protein set to use as pivot genome, to ensure the analysis returns identical results to ours (otherwise the first genome alphabetically will be used).


#### 2. Build core gene alignments and phylogenetic trees.

Generate core gene alignments and phylogenies with MAFFT and FastTree:
```
2_generate_gene_trees_core_genome.py \
	-c Methanosarcina_genus_genomes_out/core \
	-g Methanosarcina_genus_genomes \
	-ib ingroup_branch \
	--quiet 1
```

#### 3. Create codon alignment for each core gene.
```
3_create_codon_alignments.py \
	-c Methanosarcina_genus_genomes_out/core \
	--quiet 1
```


#### 4. Run CodeML on all call genes to estimate omega under two models

Note that this requires specifying the "input branch" lineage.

Also, omega must estimated with a single omega across all genomes, and also under a scenario when the two species have differing omega values.

Note that the below commands run CODEML serially on each gene (i.e., one at a time). Because of this, there are intermediate files created in the output directory with the same name for each job -- make sure you do not run the below commands at the same time, as this would cause these files to be overwritten and could create confusing bugs!

```
4_estimate_omega_core_genes.py \
	-i Methanosarcina_genus_genomes_out/core \
	-e .null \
	-m 0 \
	-n 0 \
	-c 0
```

After the above finishes, also run:
```
4_estimate_omega_core_genes.py \
	-i Methanosarcina_genus_genomes_out/core \
	-e .alt \
	-m 2 \
	-n 0 \
	-c 0
```

The above commands will produce CODEML output files for each gene family in `corecruncher_out/core`. Output files corresponding to the null model fits will end in `.null`, while those corresponding to fitting  the two models will be `.alt`.


#### 5. Run likelihood-ratio test to test for differences in omega between _Methanosarina barkeri_ and _Methanosarcina mazei_.


```
5_LRT_test_omega.py \
	-i Methanosarcina_genus_genomes_out/core \
	-n .null \
	-a .alt \
	-df 1 \
	-alpha 0.05 \
	> null_vs_alt_LRT.tsv
```

The output table (`null_vs_alt_LRT.tsv`) will look like this:
```
Core gene	LR_stat	P value	Corrected P value	Result
fam632	0.0005879999998796848	0.9806542235362545	1.0	NS
fam2350	0.0738099999998667	0.785868112457999	1.0	NS
fam2196	2.790219999999863	0.09484120606210841	1.0	NS
fam2125	0.07075199999962933	0.7902449997159229	1.0	NS
fam1062	0.0535520000012184	0.8169939273679201	1.0	NS
```

Significant core genes indicate those where separate omega values for the two tested lineages fits better than a single overall omega.

This step will not be appropriate for every analysis, as it is hard-coded for this particular comparison. However, the code could be altered to conduct similar analyses with different datasets.


#### 6. Generate overall summaries.

The raw CODEML output can be difficult to parse. To get the output in a more user-friendly format, you can run these two commands:

```
6_parse_raw_codeml_files.py \
	-i corecruncher_out/core/ \
	-o null_codeml_summary \
	-e .null

6_parse_raw_codeml_files.py \
	-i corecruncher_out/core/ \
	-o alt_codeml_summary \
	-e .alt
```

Note that this script also will not work for any arbitrary CODEML output (e.g., if more than two omega values are fit), but can be the basis for further adjustments.

This will create three outputs for the genes fit with the null and alternative models, respectively: a table giving the per-branch parameter estimates, the trees for each gene, and a table of overall summary values. These will be output within the new folders `null_codeml_summary` and `alt_codeml_summary`.


The overall summary table (`summary_tab.tsv`) is primarily of interest, and will look like this:
```
locus   kappa   omega   total_dN        total_dS        mean_N  mean_S  num_polymorphic_sites
fam100  2.45035 0.06719 0.0748  1.1139  655.0   245.0   148
fam1004 2.64685 0.03962 0.0528  1.3335  417.20000000000033      152.80000000000013      86
fam1005 2.45667 0.07525 0.1358  1.8040  712.4000000000004       262.6000000000002       198
fam1006 4.37718 0.09050 0.1140  1.2595  254.0   97.0    66
```

Each row is a separate gene, and the columns correspodn to the estimated kappa and omega parameters, as well as the observed dN and dS substitution rates (summed over all branches). The mean number of sites where non-synonymous (mean_N) and synonymous (mean_S) substitutions can occur across all sequences is also indicated. Last, "num_polymorphic_sites" indicates the total number of independent sites across all input sequences that vary, which can be a useful metric for filtering out genes with insufficient variation.

The output for the alternative model fit will include a separate column for each omega value (omega1 and omega2).

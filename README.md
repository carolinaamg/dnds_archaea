## A workflow for estimating the ratio of non-synonymous to synonymous substitutions in core protein-coding genes

This repository provides a step-by-step guide to estimate the dN/dS ratio on protein-coding genes of archeal genomes. This framework is applicable to bacteria.

Moreover, this repository provides example input files to go from genome files to the final statistical analyses to compare two models.

Input example files (genomes) can be found in the following Zenodo repository: https://zenodo.org/records/10854407

Once downloaded, your .tar.gz can be unzipped as follows:

``
tar -xzvf Methanosarcina_genus_genomes.tar.gz
``

and

``
tar -xzvf Methanosarcina_mazei_genomes.tar.gz
``

_In our tutorial, the resulting folders will be your input folders._

#### This repository is part of the book chapter: "Title", where we provide a detailed summary of the scripts shown above and a step-by-step tutorial to estimate dN/dS with the example genome files.


# A workflow for estimating the ratio of non-synonymous to synonymous substitutions in core protein-coding genes

This is a step-by-step guide to estimate the dN/dS ratio on core protein-coding genes across archeal genomes, as described in the chapter *_XXXXX_*.

This workflow starts with raw genome files and ends with the a statistical analyses to compare two evolutionary models.

The example genomes analyzed in this chapter are also available for readers to run the tutorial themselves, although this general workflow is applicable to prokaryotes in general, including bacteria.


## Installation

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
    anaconda::statsmodels
    bioconda::mafft \
    bioconda::pal2nal \
    bioconda::paml \
    bioconda::prodigal \
    bioconda::raxml \
    conda-forge::biopython \
```

To activate this environment:
```
conda activate dnds_workflow
```

_[CoreCruncher](https://github.com/lbobay/CoreCruncher)_ needs to be installed separately.

To download _CoreCruncher_:
```
git clone git@github.com:lbobay/CoreCruncher.git
```

You can add the _CoreCruncher_ scripts to your $PATH by making symbolic links to the scripts in the conda environment `bin` folder:
```
cd CoreCruncher

for SCRIPT in *.py; do
   ln -s $PWD/$SCRIPT /path/to/envs/dnds_workflow/bin/$SCRIPT
done
```

Then download [_USEARCH_ v6.1](https://drive5.com/usearch/download.html), which is required by _CoreCruncher_. Note that you can use the free version only for non-commerical use (please see the CoreCruncher documentation on how to alternatively use `BLAST`).

You will need to also add the _USEARCH_ binary to your _$PATH_, as above, and also rename the binary (indicated by the example placeholder `usearch6.1.XX` below) to be `usearch61`, which is what _CoreCruncher_ expects.
```
ln -s $PWD/usearch6.1.XX /path/to/envs/dnds_workflow/bin/usearch61
```



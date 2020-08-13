# Viral Deep-Sequencing Pipeline

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.17-brightgreen.svg)](https://snakemake.bitbucket.io)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

The `viral-deepseq` pipeline is a generalizable framework for analyzing next-generation sequencing data from RNA viruses. The workflow is flexible, allowing for multiple virus types and analysis configurations. In addition to its sequence analysis capabilities, the pipeline also functions as a toolkit with powerful scripts for analyzing intra-host genetic diversity.   

## Authors

* [Will Hannon](https://www.linkedin.com/in/williamhannon/)
* [Jesse Bloom](https://www.fredhutch.org/en/faculty-lab-directory/bloom-jesse.html)

## Getting Started 

First, clone the repository to your desired location. 

```
git clone https://github.com/jbloomlab/viral-deepseq.git
```

To get started running this pipeline, you will need `conda` as implemented in either `anaconda` or `miniconda`. To install either of these, follow the instructions [here](https://docs.conda.io/projects/continuumio-conda/en/latest/user-guide/install/index.html#regular-installation).

Once you've cloned the repository and installed conda, you can create the environment needed to run the pipeline. To do this, run the following command. 

```
conda env create --file environment.yml; conda activate viral-deepseq
```

## Running Analysis

To configure the analysis, you will need a table with the [NCBI SRA](https://www.ncbi.nlm.nih.gov/sra) accessions, the Library Layout (paired-end or single-end reads), the name of the virus corresponding to the location of its genome in the [`config file`](/config/config.yml), and the identity of the host organism if there are contaminating reads. See below for the correct format. 

| Run         | LibraryLayout | Virus | Host  |
|-------------|---------------|-------|-------|
| SRR11549941 | SINGLE        | SARS2 | human |
| SRR11549942 | SINGLE        | SARS2 | human |
| SRR11140746 | PAIRED        | SARS2 | none  |
| SRR11140748 | PAIRED        | SARS2 | none  |

To run the analysis locally, you can use the following command. I would not recommend this, because the computational time will be extensive if you don't have many cores. 

```
snakemake --use-conda --conda-prefix ./env --cores 4
```

If you plan to run the analysis on Fred Hutch `rhino` or `gizmo` servers, you can submit the pipeline to `Slurm` with the following command. 

```
sbatch run_analysis.bash
```
# SARS-CoV-2 Transmission Bottleneck

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.17-brightgreen.svg)](https://snakemake.bitbucket.io)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

*This analysis was conducted using the `viral-deepseq` pipeline template.*

## Contributors

* [Will Hannon](https://www.linkedin.com/in/williamhannon/)
* [Jesse Bloom](https://www.fredhutch.org/en/faculty-lab-directory/bloom-jesse.html)

## Analysis Overview

This analysis aims to see if there is a substantial transmission of intra-host minor variants between individuals in a large cluster of infections. This information could help characterize the inter-host transmission bottleneck in SARS-CoV-2 with the highest resolution achieved so far. 

## Running Analysis

The analysis is configured in the `./config` directory. Reference genomes and associated gene annotations are specified in the [`config file`](/config/config.yml), and the location of the included samples is found in the `samples.csv` file.  

To run the analysis locally, you can use the following command. I would not recommend this, because the computational time will be extensive if you don't have many cores. 

```
snakemake --use-conda --conda-prefix ./env --cores 4
```

If you plan to run the analysis on Fred Hutch `rhino` or `gizmo` servers, you can submit the pipeline to `Slurm` with the following command. 

```
sbatch run_analysis.bash
```
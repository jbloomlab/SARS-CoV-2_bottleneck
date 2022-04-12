# Narrow transmission bottlenecks and limited within-host viral diversity during a SARS-CoV-2 outbreak on a fishing boat

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.17-brightgreen.svg)](https://snakemake.bitbucket.io)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/287312783.svg)](https://zenodo.org/badge/latestdoi/287312783)

## Contributors

* [Will Hannon](https://www.linkedin.com/in/williamhannon/)
* [Jesse Bloom](https://www.fredhutch.org/en/faculty-lab-directory/bloom-jesse.html)

## Authors
William W. Hannon<sup>1,2*</sup>, Pavitra Roychoudhury<sup>4,5*</sup>, Hong Xie<sup>5</sup>, Lasata Shrestha<sup>5</sup>, Amin Addetia<sup>1,5</sup>, Keith R. Jerome<sup>4,5</sup>, Alexander L. Greninger,<sup>4,5</sup> Jesse D. Bloom<sup>2,3,6,#</sup>
 
<sup>1</sup> Molecular and Cellular Biology Graduate Program, University of Washington, Seattle, WA 98109
<sup>2</sup> Basic Sciences and Computational Biology, Fred Hutchinson Cancer Research Center, Seattle, WA 98109
<sup>3</sup> Department of Genome Sciences, University of Washington, Seattle, WA 98109
<sup>4</sup> Vaccine and Infectious Disease Division, Fred Hutchinson Cancer Research Center, Seattle, Washington, USA
<sup>5</sup> Department of Laboratory Medicine and Pathology, University of Washington School of Medicine, Seattle, Washington, USA
<sup>6</sup> Howard Hughes Medical Institute, Seattle, WA 98109

## Abstract

The long-term evolution of viruses is ultimately due to viral mutants that arise within infected individuals and transmit to other individuals. Here we use deep sequencing to investigate the transmission of viral genetic variation among individuals during a SARS-CoV-2 outbreak that infected the vast majority of crew members on a fishing boat. We deep-sequenced nasal swabs to characterize the within-host viral population of infected crew members, using experimental duplicates and strict computational filters to ensure accurate variant calling. We find that within-host viral diversity is low in infected crew members. The mutations that did fix in some crew members during the outbreak are not observed at detectable frequencies in any of the sampled crew members in which they are not fixed, suggesting viral evolution involves occasional fixation of low-frequency mutations during transmission rather than persistent maintenance of within-host viral diversity. Overall, our results show that strong transmission bottlenecks dominate viral evolution even during a superspreading event with a very high attack rate. 

## Running Analysis

The analysis is configured in the `./config` directory. Reference genomes and associated gene annotations are specified in the [`config file`](/config/config.yml), and the path to the sequencing samples is found in the `samples.csv` file.  

To run the analysis locally, you can use the following command. I would not recommend this, because the computational time will be extensive if you don't have many CPUs. 

```
snakemake --use-conda --conda-prefix ./env --cores 4
```

If you plan to run the analysis on Fred Hutch `rhino`, you can submit the pipeline to `Slurm` with the following command. 

```
sbatch run_analysis.bash
```

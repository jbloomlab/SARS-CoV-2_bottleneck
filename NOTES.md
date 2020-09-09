## Notes

This document contains notes and useful information for the project – `SARS-CoV-2_Bottleneck`.

## Pre-print Notes

#### Cohort Overview

There was an outbreak on a fishing boat in Seattle. Individuals (120/122) were screened for SARS-CoV-2 RNA by RT-qPCR and seroconversion by the Abbot IgG assay both before departure and after return. Although all were negative, an outbreak occured nonetheless. A total of 104/122 crew members had viral RNA or seroconverted upon return to shore (84% Attack Rate). 

#### Data Overview 

From the 104 reportedly infected crew memebers there are 39 metagenomic RNA-sequencing runs avaliable. The data is avaliable publically at NCBI BioProject `PRJNA610428`. Metadata for the identity of samples is avaliable in Supplementary Table 2.

The metadata only contains the GISAID-ID for the consensus and the strain. I might be able to use the strain to figure out the SRR# using the metadata on the [SRA run selector](https://www.ncbi.nlm.nih.gov/Traces/study/?query_key=2&WebEnv=MCID_5f344f743a76fab9ab22e42b&o=acc_s%3Aa). 

#### Experimental Overview 

> Libraries were sequenced on a 1x75 bp Illumina NextSeq run. A median of 509,551 sequencing reads were obtained for each sample.

#### Pipeline Overview

In the paper, the authors use a [custom genome calling pipeline](https://github.com/proychou/hCoV19) to construct consensus genome from the samples for use in phlyogenetic analyses. Briefly, they trimmed adaptors with `BBDuk`, aligned reads to the SARS-CoV-2 reference genome (Wuhan-Hu-1 NC_045512.2) with `Bowtie2`, filtered viral reads with `BBDuk`, and assembled new consensus genomes with `SPAdes`. 

#### Results Overview

*This is a quick overview of the relevant results from the section "Confirmation of outbreak with whole genome sequencing"*

> Metagenomic recovery of 39 SARS-CoV-2 whole genomes from the outbreak indicated a major single outbreak clade (FastTree support value: 1.00) covering 38 isolates that differed by a median of one nucleotide across the genome (range 0-5). Sixteen of these isolates shared completely identical sequence.

## Project Notes

The goal of this preliminary analysis is to identify homoplastic, low-frequency variants that might reveal something about the size of the inter-patient transmission bottleneck for SARS-CoV-2. 

There are three useful papers to go over closely when thinking about how to construct this analysis: 

1. [Shared SARS-CoV-2 diversity suggests localised transmission of minority variants](https://doi.org/10.1101/2020.05.28.118992)
2. [Population Bottlenecks and Intra-host Evolution during Human-to-Human Transmission of SARS-CoV-2](https://doi.org/10.1101/2020.06.26.173203)
3. [Stochastic processes constrain the within and between host evolution of influenza virus](https://doi.org/10.7554/eLife.35962)

Duplicates by re-sequencing: 

For 24 samples (with Ct < 20),  Pavitra re-prepped and re-sequenced them using the same shotgun metagenomics approach as before. 

This resulted in 24 total samples (each with duplicates: the original sequencing on the SRA, and the new sequencing). 

One of the samples is missing a corresponding run: SpID 10110 1 10
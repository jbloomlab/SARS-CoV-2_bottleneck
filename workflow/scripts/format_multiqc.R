## ---------------------------
##
## Script name: `rule format_multiqc`
##
## Purpose of script: Takes multiqc tables and converts them into a single 
## formatted table for further analysis.
##
## Author: Will Hannon
##
## Date Created: 2021-01-15
##
## Copyright (c) Will Hannon, 2020
##
## Email: wwh22@uw.edu
##
## TODO: Figure out how to get/calculate the % of filtered reads from BBTools
##
## ---------------------------

require(tidyverse)

## ==== Get file paths from snakemake object ==== ##

# TODO: Convert these to Snakemake objects for integration
samtools.stats = "results/qc/multiqc/multiqc_report_data/multiqc_samtools_stats.txt"
fastqc.stats = "results/qc/multiqc/multiqc_report_data/multiqc_fastqc.txt"
picard.stats = "results/qc/multiqc/multiqc_report_data/multiqc_picard_dups.txt"
alignment.stats = "results/qc/multiqc/multiqc_report_data/mqc_samtools_alignment_plot_1.txt"
sources = "results/qc/multiqc/multiqc_report_data/multiqc_sources.txt"

## ==== Process the multiqc data files ==== ##

# Data for merged and unmerged BAMs from `samtools stats`
samtools.df = read_tsv(samtools.stats) %>% 
  # Get the metadata from the sample names
  separate(Sample, into = c("sample", "replicate", "library", "aligner"), fill = "right") %>% 
  # Fix the case where merged samples don't have the correct columns 
  mutate(library = if_else(library == "BWA", "merged", library)) %>% 
  # Select only useful columns
  select(sample, replicate, library,
         raw_total_sequences, reads_mapped,
         reads_unmapped, total_length, bases_mapped,
         mismatches, average_length, error_rate,
         average_quality, reads_mapped_percent,
         reads_unmapped_percent, reads_QC_failed_percent)


# Fastqc data for only the trimmed reads
fastqc.df = read_tsv(fastqc.stats) %>% 
  # Get the metadata from the sample names
  separate(Sample, into = c("sample", "replicate", "library"), fill = "right") %>% 
  # Select the useful columns
  select(sample, replicate, library,
         fastqc_raw_total_sequences = "Total Sequences",
         poor_quality_sequences = "Sequences flagged as poor quality",
         per_base_n_content, per_sequence_quality_scores)
  

# Data from Picard MarkDuplicates
picard.df = read_tsv(picard.stats) %>% 
  # Get the metadata from the sample names
  separate(Sample, into = c("sample", "replicate", "library", "aligner"), fill = "right") %>% 
  # Select the useful columns
  select(sample, replicate, library,
         duplicates = "UNPAIRED_READ_DUPLICATES",
         percent_duplicated = "PERCENT_DUPLICATION")

# Join everything together
qc.df = left_join(samtools.df, full_join(picard.df, fastqc.df, by = c("sample", "replicate", "library")), by = c("sample", "replicate", "library"))

## ==== Export the results ==== ##

write_csv(qc.df, snakemake@output[[1]])

## ==== End of script ==== ##
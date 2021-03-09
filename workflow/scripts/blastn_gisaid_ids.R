## ---------------------------
##
## Script name: `blastn_gisaid_ids``
##
## Purpose of script: Get a sequences from clade 20C to make a blastn database
##
## Author: Will Hannon
##
## Date Created: 2021-03-04
##
## Copyright (c) Will Hannon, 2021
## Email: wwh22@uw.edu
##
## ---------------------------

require(tidyverse)

## ---------------------------

# IMPORTANT: this is randomized, so set the seed for reproduciblity 
set.seed(20210120)

# Set the max collection date - currently the middle of June (samples were collected by May 29th)
max.collection.date = "2020-06-15"

# Process the data downloaded from gisaid. This contains all metadata as of 2021-02-22
GISAID.df = read_tsv("../../config/gisaid/metadata_2021-03-04_10-31.tsv") %>% 
  # Convert the collection date into a date object
  mutate(collection_date = as.Date(date, format = "%Y-%m-%d")) %>% 
  # Remove anything without a conforming collection date
  filter(!is.na(collection_date)) %>% 
  # Filter sequences from before the end of May (as per the time the boat samples were collected)
  # Also, filter out sequences that are substatially different from the length of the reference
  filter(collection_date <= max.collection.date & host == "Human" & (length > 29000 & length < 30000))

# Process the data from Addieta, et. al., 2020 to exclude these from the selected genomes
Addieta.seqs.df = read_csv("../../config/addieta-et-al_supplemental-table-2.csv")
Addieta.seqs = Addieta.seqs.df$`GISAID Accession Number`

# All of the samples excluding those used in the paper. 
all.clades.df = GISAID.df[which(!GISAID.df$gisaid_epi_isl %in% Addieta.seqs),]
addieta.seqs.df = GISAID.df[which(GISAID.df$gisaid_epi_isl %in% Addieta.seqs),] %>% 
  mutate(study = "internal")

# All of the sequences that meet quality standards in clade 20C for blastn database (excludes the boat seqs)
write.csv(select(filter(all.clades.df, Nextstrain_clade == '20C'), gisaid_epi_isl, strain), "../../config/gisaid/bastn_gisaid_ids.csv", row.names = F)

# The full metadata for the boat samples
write.csv(addieta.seqs.df, "../../config/gisaid/boat_gisaid_ids.csv", row.names = F)


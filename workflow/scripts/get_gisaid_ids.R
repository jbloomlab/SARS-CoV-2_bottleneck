## ---------------------------
##
## Script name: `get_gisaid_ids``
##
## Purpose of script: Get a random samples of 100 genomes from each clade before a specific date
##
## Author: Will Hannon
##
## Date Created: 2021-02-23
##
## Copyright (c) Will Hannon, 2021
## Email: wwh22@uw.edu
##
## ---------------------------

require(tidyverse)

## ---------------------------

# IMPORTANT: this is randomized, so set the seed for reproduciblity 
set.seed(20210120)

# How many random genomes to samples
n_genomes = 100

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

# Similar strains from BlastN
similar.strains = read_csv("../../config/gisaid/smilar_to_boat_genomes.csv")$strain

# All of the samples excluding those used in the paper. 
all.clades.df = GISAID.df[which(!GISAID.df$gisaid_epi_isl %in% Addieta.seqs),]
all.clades.df = all.clades.df[which(!all.clades.df$strain %in% similar.strains),]

# Get the clades to sample (Any clade represented by more than 1000 sequences globablly at that time)
clades.to.samples = all.clades.df %>% 
  filter(!is.na(Nextstrain_clade)) %>% 
  group_by(Nextstrain_clade) %>% 
  count() %>% 
  filter(n > 1000) %>% 
  pull(Nextstrain_clade)

# Randomly select 100 sequences from each of the major clades; 19A, 19B, 20A, 20B, 20C, 20D
selected.clades.df = data.frame()

for (i in 1:length(clades.to.samples)) {
  
  selection = all.clades.df %>% 
    filter(Nextstrain_clade == clades.to.samples[i]) %>% 
    sample_n(., n_genomes)
  
  selected.clades.df = rbind(selected.clades.df, selection)
  
}

# The GISAID ids that contain a representation of each clade and the boat samples
write.csv(selected.clades.df, "../../config/gisaid/phylogeny_gisaid_ids.csv", row.names = F)








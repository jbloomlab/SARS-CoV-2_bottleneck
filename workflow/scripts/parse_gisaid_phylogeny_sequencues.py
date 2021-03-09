"""
Parse all of the sequences from GISAID to assembly a multifasta file
for alignment and subsequent tree generation. 
"""
__author__ = "Will Hannon"
__copyright__ = "Copyright 2021 Will Hannon"
__email__ = "wwh22@uw.edu"
__license__ = "MIT"

from Bio import SeqIO
import pandas as pd

# File name for metadata and fasta
all_gisaid_sequences = "../../config/gisaid/sequences_2021-03-04_08-34.fasta"
metadata_path = "../../config/gisaid/metadata_2021-03-04_10-31.tsv"
representative_strains_path = "../../config/gisaid/phylogeny_gisaid_ids.csv"
boat_strains_path = "../../config/gisaid/boat_gisaid_ids.csv"
similar_strains_path = "../../config/gisaid/smilar_to_boat_genomes.csv"


# Import the metadata as a csv
metadata_df = pd.read_csv(metadata_path, sep='\t', low_memory=False) 
boat_samples_df = pd.read_csv(boat_strains_path) 
similar_samples_df = pd.read_csv(similar_strains_path) 

# List of strains for each category
samples = set(pd.read_csv(representative_strains_path).strain.to_list()) 

# Make a list to save which samples end up in the final fasta
strain_list = []

# Filter the sequences that have low complexity
for record in SeqIO.parse(all_gisaid_sequences, 'fasta'):
    
    # Get the records
    header = record.id
    genome = str(record.seq)

    # Check if it's an acceptable sample
    if header in samples:

        # Check how many N's (less than 5% N)
        if genome.count('N')/len(genome) < 0.05: 
            strain_list.append(header)

included_df = metadata_df[metadata_df.strain.isin(strain_list)].groupby('Nextstrain_clade').apply(pd.DataFrame.sample, n=30, random_state=1).drop(columns = "Nextstrain_clade").reset_index().drop(columns = "level_1")

final_samples_df = pd.concat([included_df, boat_samples_df , similar_samples_df ])

with open('../../config/gisaid/phylogeny.fasta', 'w') as outfasta:

    for record in SeqIO.parse(all_gisaid_sequences, 'fasta'):
        # Get the records
        header = record.id
        genome = str(record.seq)

        # Check if it's an acceptable sample
        if header in set(final_samples_df.strain):
            outfasta.write(">" + header + "\n" + genome + "\n")

final_samples_df.to_csv('../../config/gisaid/phylogeny_metadata.csv', index=False)

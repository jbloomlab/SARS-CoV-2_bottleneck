"""
Parse all of the GISAID sequences to get those that can be included in a 
BLASTN database. This code also gets all of the fasta sequences from the
boat to use as a query against the blast database. 
"""
__author__ = "Will Hannon"
__copyright__ = "Copyright 2021 Will Hannon"
__email__ = "wwh22@uw.edu"
__license__ = "MIT"

from Bio import SeqIO
import pandas as pd


# File name for metadata and fasta
all_gisaid_sequences = "../../config/gisaid/sequences_2021-03-04_08-34.fasta"
samples_to_make_database = "../../config/gisaid/bastn_gisaid_ids.csv"
boat_samples = "../../config/gisaid/boat_gisaid_ids.csv"

# Import the metadata as a csv
clade_20C_samples = set(pd.read_csv(samples_to_make_database).strain.to_list())
boat_ids = set(pd.read_csv(boat_samples).strain.to_list()) 

with open('../../config/gisaid/blastn_covid_database.fasta', 'w') as outfasta:
    # Iterate over the records in the multifasta file
        for record in SeqIO.parse(all_gisaid_sequences, 'fasta'):
            # Get the records
            header = record.id
            genome = str(record.seq)

            # Check if it's an acceptable sample
            if header in clade_20C_samples:

                # Check how many N's (less than 5% N)
                if genome.count('N')/len(genome) < 0.05: 
                    outfasta.write(">" + header + "\n" + genome + "\n")


with open('../../config/gisaid/boat_genomes.fasta', 'w') as outfasta:
    # Iterate over the records in the multifasta file
        for record in SeqIO.parse(all_gisaid_sequences, 'fasta'):
            # Get the records
            header = record.id
            genome = str(record.seq)

            # Check if it's an acceptable sample
            if header in boat_ids:
                outfasta.write(">" + header + "\n" + genome + "\n")


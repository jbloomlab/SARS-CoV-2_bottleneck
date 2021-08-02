"""
The purpose of this script is to make a high-confidence consensus sequence
from two independent sequencing replicates of the same patient. This 
approach only takes SNPs into account.

The approach makes a consensus sequence from each replicate by using 
the majority vote of an allele at each position. There is a filter for
bases with a quality below a specified value. There is also a filter for 
minimum coverage. Anythin position with less than 100X reads will be called
as an N. After the individual consensus sequences are generated, the
missing gaps (Ns) due to low coverage are filled in with high coverage bases
from the matched replicate. Finally, a global mask is made by combining
regions covered by Ns in all samples.
"""

__author__ = "Will Hannon"
__copyright__ = "Copyright 2020 Will Hannon"
__email__ = "wwh22@uw.edu"
__license__ = "MIT"

import pysam # count variant alleles from BAM
import pandas as pd # data frames
import numpy as np # arrays
import os # interacting with files
from Bio import SeqIO # reading fasta format
import re # regular expressions
from collections import Counter # count the frequeny of bases

## === Functions === ##
def get_ref_sequence(refpath):
    """
    Function to read in the reference sequence as a list.
    """
    return [base.upper() for base in list(SeqIO.parse(refpath, "fasta"))[0].seq]


def get_sample_info(filepath):
    """
    Get the sample and the replicate from the path.
    """
    metadata = re.split("\\.|_", os.path.basename(filepath))
    return (metadata[0], metadata[1])


def check_read(read):
    """
    Helper function to decide what reads should
    be keep when parsing alignment file with `pysam`. 
    """
    # Exclude Quality Failures
    if read.is_qcfail:
        return False
    # Exclude Secondary Mappings
    if read.is_secondary:
        return False
    # Exclude Unmapped Reads
    if read.is_unmapped:
        return False
    else:
        return True

    
def build_consensus_seq(
                bampath, 
                callback_function = check_read, 
                contig = "NC_045512.2", 
                refpath = "../../config/ref/SARS2.fa", 
                minimum_QUAL = 25, 
                minimum_COV = 100):
    """
    Read in BAM file and convert to a consensus fasta the length of the
    reference genome using the `pysam` command `count_coverage`. 
    """
     # String to hold the final consensus seq
    consensus = ""
    
    # Open alignment with pysam
    with pysam.AlignmentFile(bampath, "rb") as bamfile:
        # Get a dataframe of the counts
        count_df = pd.DataFrame.from_dict({base:counts for base, counts in zip("ACGT", bamfile.count_coverage(contig = contig, read_callback=callback_function, quality_threshold=minimum_QUAL))})
        # Add the depth at each position
        count_df['DP'] = count_df.sum(axis = 1)
        # Add the position 
        count_df['POS'] = count_df.index + 1
        # Add the reference allele
        count_df['REF'] = get_ref_sequence(refpath)
        # Convert counts to frequency 
        count_df.iloc[:,0:4] = count_df.iloc[:,0:4].div(count_df.DP, axis = 0)
        # Handle any NaNs created by dividing by 0 coverage
        count_df = count_df.fillna(0)
        
    # Iterate over each row and determine the conesnus sequence and build consensus string
    for index, row in count_df.iterrows(): 
        if row.DP < minimum_COV:
            consensus += "N"
        else: 
            consensus += max([(row["A"], "A"), (row["C"], "C"), (row["G"], "G"), (row["T"], "T")], key=lambda x:x[0])[1]

    return consensus


def merge_consensus_seq(consensus_1, consensus_2, sample, refpath):
    """
    Combine the consensus sequences between the two replicates using higher coverage
    reigons of each replicate to fill in gaps. 
    """
    
    # Get the reference sequence. 
    ref = get_ref_sequence(refpath)
    
    # Hold the merged consensus from the two replicates. 
    merged_consensus = ""
    
    # List to hold positions where the two replicates disagree.
    discrepencies = []
    
    # Iterate over ever position and compare depth and allele to merge.
    for pos, bases in enumerate(zip(consensus_1, consensus_2)): 

        if len(set(bases)) == 1 and 'N' not in bases:
            merged_consensus += bases[0]
        elif len(set(bases)) == 1 and 'N' in bases:
            merged_consensus += ref[pos]
        elif len(set(bases)) > 1:
            if 'N' in bases and bases[0] == 'N':
                merged_consensus += bases[1]
            elif 'N' in bases and bases[1] == 'N':
                merged_consensus += bases[0]
            else:
                merged_consensus += ref[pos]
       
    return merged_consensus
    
    

def main(): 
    """
    Main function to run the analysis.
    """
    
    ## ======== Input data ======== ##
    
    # Get the path list for the input BAM files from snakemake rule.
    bampaths = snakemake.input.bam
    # Get the path to the output multi-fasta file. 
    outpath = str(snakemake.output)
    # Get the path to the reference genome.
    refpath = str(snakemake.input.genome)
    
    ## ======== Generate individual consensus sequences ======== ##
    
    # Iterate over the samples and make a rough consensus sequence. 
    consensus_dict = {get_sample_info(bampath):build_consensus_seq(bampath, refpath=refpath) for bampath in bampaths}
    
    ## ======== Merge individual consensus sequences ======== ##
    
    samples = list({name_tup[0] for name_tup in consensus_dict.keys()})
    
    polished_consensus_dict = {sample:merge_consensus_seq(consensus_dict[(sample, '1')], consensus_dict[(sample, '2')], sample, refpath) for sample in samples}
    
    ## ======== Write out to a fasta file ======== ##
    
    with open(outpath, "w") as outfile:
        for sample, consensus in polished_consensus_dict.items():
            outfile.write(">" + sample + "\n" + consensus + "\n")

if __name__ == '__main__': 
    main()

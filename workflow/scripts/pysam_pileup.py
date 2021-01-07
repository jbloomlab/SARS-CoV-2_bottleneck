"""
This module contains functions for calling single nucleotide 
polymorphisms from Illumina short read sequencing data in the 
form of indexed and sorted BAM files. 

The main library that facillitates this is `pysam`.

Importantly, this module's functions differ from other
variants callers in that there is no statistical filtering
of variants. All variants are identified and reported consistent
with heuristic cutoffs. 
"""
__author__ = "Will Hannon"
__copyright__ = "Copyright 2020 Will Hannon"
__email__ = "wwh22@uw.edu"
__license__ = "MIT"

## === Import Libraries === ##
import pysam #count variant alleles from BAM
import pandas as pd #data frames
import numpy as np #arrays
import os #interacting with files
from Bio import SeqIO #reading fasta format
import re #regular expressions

## === Functions === ##

def translate(codon):
    """
    Translate a three letter DNA string into 
    a one letter amino acid code. 

    Parameters
    ----------
    codon : str
        three letter DNA sequence

    Returns
    -------
    str
        one letter amino acid code

    Raises
    ------
    AssertionError
        error if codon sequence is invalid
        
    """
    
    table = { 
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                  
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
        'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*', 
        'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W', 
    } 
    
    assert codon in table.keys(), "Not a valid codon sequence."
    
    return table[codon]


def mutate(codon, alt, index):
    """
    Replace (mutate) a base in a codon with an 
    alternate base. 
        
    Parameters
    ----------
    codon : str
        three letter DNA sequence
        
    alt : str
        alternative base
        
    index : int
        index of the alt base in codon (0|1|2). 

    Returns
    -------
    str
        codon with alternative base

    Raises
    ------
    AssertionError
        error if index is not valid (0|1|2)
        
    AssertionError
        error if base is not valid (A|T|C|G)
        
    """
    
    assert index in [0,1,2], "Not a valid index."
    
    assert alt in ["A", "T", "C", "G"], "Not a valid base."
    
    return "".join([alt if i == index else b for i,b in enumerate(codon)])


def check_read(read):
    """
    Helper function to decide what reads should
    be keep when parsing alignment file with `pysam`. 

    Parameters
    ----------
    read : AlignedSegment
        read from alignment file parsed with `pysam`.

    Returns
    -------
    bool
        True/False if read should be included
        
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
        
    
def build_af_df(filepath, 
                callback_function = check_read, 
                ref = "NC_045512.2", 
                ref_path = "../../config/ref/SARS2.fa", 
                minimum_AF = 0.01, 
                minimum_qual = 25):
    """
    Read in BAM file and convert to a dataframe containing the frequency of 
    of any bases present at a given position in the reference genome using the
    `pysam` command `count_coverage`. 

    Parameters
    ----------
    filepath : str
        path to the bam file to be parsed
        
    callback_function : function
        function that decides which reads to keep/exclude
        
    ref : str
        name of the contig to count coverage over
        
    ref_path : str
        path to the reference genome as fasta
        
    minimum_AF : float
        minimum allele frequency to include allele in dataframe
        
    minimum_qual : int
        minimum QUAL score to count read at a position.

    Returns
    -------
    Pandas.DataFrame
       Data Frame containing the bases represented at each positon in the genome
        
    """

    with pysam.AlignmentFile(filepath, "rb") as bamfile:
        
        # The count_coverage method counts the occurances of each base at each position. It excludes reads based on the callback function
        count_df = pd.DataFrame.from_dict({base:counts for base, counts in zip("ACGT", bamfile.count_coverage(contig = ref, read_callback=callback_function, quality_threshold=minimum_qual))})
        
        # Add the depth at each position
        count_df['DP'] = count_df.sum(axis = 1)
        
        # Add the position 
        count_df['POS'] = count_df.index + 1
        
        # Add the reference allele
        count_df['REF'] = [base.upper() for base in list(SeqIO.parse(ref_path, "fasta"))[0].seq]
        
        # convert counts to frequency 
        count_df.iloc[:,0:4] = count_df.iloc[:,0:4].div(count_df.DP, axis = 0)
        
        # handle any NaNs created by dividing by 0 coverage
        count_df = count_df.fillna(0)
        
        # Melt the data frame to a longer ('tidy') form
        count_df = pd.melt(count_df, 
                           id_vars=['POS', 'DP', 'REF'],
                           value_vars=[base for base in 'ATGC'],
                           value_name='AF',
                           var_name='ALT')
        
        # Filter out anything less than minimum allele freq
        count_df = count_df[count_df['AF'] >= minimum_AF]
        
        # TRUE/FALSE if it's a SNP
        count_df['SNP'] = np.where(count_df['ALT'] != count_df['REF'], True, False)
        
        # Is a base consensus or not.
        count_df['CONS'] = count_df['AF'].map(lambda x: x >= 0.5)
    
        # Sort by position and reset the index
        return count_df.sort_values('POS').reset_index(drop=True)
     

def main():
    """Main function to call pysam pileup scipt.
    """
    # Get the path for the input BAM file from snakemake rule.
    inputpath = str(snakemake.input.bam)

    # Get the path to the output csv file from snakemake rule.
    outpath = str(snakemake.output)

    # Get the path to the reference genome from snakemake rule. 
    ref_genome = [base.upper() for base in list(SeqIO.parse(str(snakemake.input.genome), "fasta"))[0].seq]

    # Get the sample name, i.e. accession, from the path.
    accession = os.path.basename(inputpath)

    # Build the count df
    # TODO: Currently hardcoded for COVID Wuhan-1 reference build.
    count_df = build_af_df(inputpath, 
                callback_function=check_read, 
                ref = "NC_045512.2", 
                ref_path = str(snakemake.input.genome), 
                minimum_AF = 0.01, 
                minimum_qual = int(snakemake.params.score))

    # Add accession and day to data frame
    count_df.insert(0, 'ACCESSION', accession)

    # Exporting to csv for final analysis in `R`
    count_df.to_csv(outpath, index=False)


if __name__ == '__main__':
    main()


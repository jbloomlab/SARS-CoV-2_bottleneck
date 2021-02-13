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
import annotate_coding_changes as annotate

## === Functions === ##

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
        
    
def identify_snps(filepath, 
                callback_function = check_read, 
                ref = "NC_045512.2", 
                ref_path = "../../config/ref/SARS2.fa", 
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
         
    minimum_qual : int
        minimum QUAL score to count base observation at a position.

    Returns
    -------
    Pandas.DataFrame
       Data Frame containing the bases represented at each positon in the genome
        
    """

    with pysam.AlignmentFile(filepath, "rb") as bamfile:
        
        # The count_coverage method counts the occurances of each base at each position. 
        # It excludes reads based on the callback function.
        count_df = pd.DataFrame.from_dict({base:counts for base, counts in zip("ACGT",
                                                                               bamfile.count_coverage(contig=ref,
                                                                                                      read_callback=callback_function,
                                                                                                      quality_threshold=minimum_qual))})
        # Add the depth at each position
        count_df['DP'] = count_df.sum(axis = 1)
        # Add the position (1-indexed)
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
        # Remove anything with 0 coverage.
        count_df = count_df[count_df['AF'] > 0]
        # TRUE/FALSE if it's a SNP
        count_df['SNP'] = np.where(count_df['ALT'] != count_df['REF'], True, False)
        # TRUE/FALSE if it's a consensus base.
        count_df['CONS'] = count_df['AF'].map(lambda x: x >= 0.5)
        # Add the number of times a given allele is observed.
        count_df['OBSV'] = count_df.DP * count_df.AF
        
        # Filter to only get SNPs
        count_df = count_df.loc[count_df['SNP']]

        # Iterate over the pileup to get info on read position, base quality, and strand bias
        snp_set = {pair for pair in zip(count_df.POS, count_df.ALT)}
        pos_set = {pos for pos, alt in snp_set}
        # Dictionary to store information for each alternative allele
        allele_dict = dict()
        # Get the pileup column for a specific region
        for pileupcolumn in bamfile.pileup(maxdepth = 0, stepper = 'nofilter', flag_filter = 0, min_base_quality = minimum_qual):
            # Only iterate through reads that have alternative alleles
            if pileupcolumn.reference_pos + 1 not in pos_set:
                continue
            # Iterate over every alignment in the pileup column. 
            for pileupread in pileupcolumn.pileups:
                # Check if the read is valid and can be parsed
                if check_read(pileupread.alignment) and not pileupread.is_del and not pileupread.is_refskip:
                    # Save the 1-indexed position
                    pos = pileupcolumn.reference_pos + 1
                    # Position in the read
                    readpos = pileupread.query_position
                    # Base at this position in the read
                    base = pileupread.alignment.query_sequence[readpos]
                    # Base quality
                    qual = pileupread.alignment.query_qualities[readpos]
                    # Orientation
                    if pileupread.alignment.is_reverse:
                        orientation = 1
                    else:
                        orientation = 0
                    # Check if this read has a SNP  
                    if (pos, base) in snp_set: 
                        # Add this SNP to a dict to hold the read postions, qualities, and orientation
                        if (pos, base) not in allele_dict.keys():
                            allele_dict[(pos, base)] = ([readpos], [qual], [orientation])
                        else: 
                            allele_dict[(pos, base)][0].append(readpos)
                            allele_dict[(pos, base)][1].append(qual)
                            allele_dict[(pos, base)][2].append(orientation)
        # List to hold data for readposition, quality, and orientation
        filter_data = []
        # Get the average or ratio of each value as list of tuples to add back to main dict.
        for key, value in allele_dict.items(): 
            pos = key[0]
            alt = key[1]
            mean_readpos = sum(int(x) for x in value[0])/len(value[0])
            mean_qual = sum(int(x) for x in value[1])/len(value[1])
            strand_ratio = sum(int(x) for x in value[2])/len(value[2])
            filter_data.append((pos, alt, mean_readpos, mean_qual, strand_ratio))
        # Save dict with pileup info
        pileup_df = pd.DataFrame(filter_data, columns=('POS', 'ALT', 'MEAN_READPOS', 'MEAN_QUAL', 'STRAND_RATIO'))
        
        # Merge the pileup dict and the snp dict. 
        return pd.merge(count_df, pileup_df,  how='left', left_on=['POS','ALT'], right_on = ['POS','ALT']).sort_values('POS').reset_index(drop=True)

    
def coverage_filter(snp_df, threshold = 10): 
    """
    Filter rows that don't meet the heuristic 
    threshold of coverage specified by the user. 
    
    default is set to 10X the reciprocal of the 
    allele frequency.
    """
    return snp_df.assign(COV_FILTER = lambda df: df.DP >= (1/df.AF * threshold))


def observation_filter(snp_df, threshold = 2): 
    """
    Filter for the minimum number of reads where a SNP was observed.
    
    The default is 2 reads containing a SNP.
    """
    return snp_df.assign(OBSV_FILTER = lambda df: df.OBSV >= threshold)


def read_position_filter(snp_df, avg_read_length = 71, percentile = 0.9):
    """
    Filter out SNPs that at the extreme start or end of reads. 
    """
    maxavg = avg_read_length*percentile
    minavg = avg_read_length*(1-percentile)
    return snp_df.assign(READPOS_FILTER = lambda df: (df.MEAN_READPOS >= minavg) & (df.MEAN_READPOS <= maxavg))


def strand_bias_filter(snp_df, threshold = 0.9):
    """
    Filter out SNPs with a strand bias
    """
    return snp_df.assign(STRANDBIAS_FILTER = lambda df: (df.STRAND_RATIO >= (1-threshold)) & (df.STRAND_RATIO <= threshold))


def mask(snp_df, mask = [x for x in range(29860, 29904)]):
    """
    Define a mask of position to exclude in the analysis. 
    
    These could be low complexity, close to the start or 
    end of the genome, etc.. 
    """
    return snp_df.loc[~snp_df.POS.isin(mask)]


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
    accession = os.path.basename(inputpath).split(".")[0]

    # Build the snp df
    snp_df = identify_snps(inputpath, 
                callback_function=check_read, 
                ref = "NC_045512.2", 
                ref_path = str(snakemake.input.genome), 
                minimum_qual = int(snakemake.params.score))

    # Add accession and day to data frame
    snp_df.insert(0, 'ACCESSION', accession)

    # Filter based on snakemake params
    strand_bias_threshold = 0.9
    read_position_threshold = 0.9
    avg_read_length = 71
    observation_threshold = 2
    coverage_threshold = 10

    snp_df = coverage_filter(snp_df, threshold = coverage_threshold)

    snp_df = observation_filter(snp_df, threshold = observation_threshold)

    snp_df = read_position_filter(snp_df, avg_read_length = avg_read_length, percentile = read_position_threshold)

    snp_df = strand_bias_filter(snp_df, threshold = strand_bias_threshold)

    snp_df = mask(snp_df)

    # TODO: Annotate coding changes
    # Load in the information from genbank
    cds_dict, genome = annotate.parse_genbank(email = "wwh22@uw.edu", ref_id = "NC_045512.2")

    snp_df[['CODON_POS', 'AA_CHANGE', 'EFFECT', 'GENE']] = snp_df.apply(lambda df: annotate.annotate_effect(cds_dict, genome, (df.POS, df.ALT)), axis=1, result_type='expand')

    # Exporting to csv for final analysis in `R`
    snp_df.to_csv(outpath, index=False)


if __name__ == '__main__':
    main()


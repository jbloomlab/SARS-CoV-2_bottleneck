"""
This script downloads fastq files for viral deep-sequencing Snakemake workflows. 

It does this in one of two ways:
    (1) Moving files from a safe location on a local machine into the working repo.
    (2) Downloading files from the NCBI SRA using `fasterq-dump`.

The output is a gzipped fastq file. 
"""
_author__ = "Will Hannon"
__copyright__ = "Copyright 2020, Will Hannon"
__email__ = "wwh22@uw.edu"
__license__ = "MIT"

import os
import pandas as pd
from snakemake.shell import shell

# Get the list of all samples from the config file.
SAMPLE_DF = pd.read_csv(snakemake.config['samples']['file'])

def get_source(name, library):
    """
    Using the sample name, determine whether the fastq file
    is on the SRA or located in your file system. Do this
    by referencing the sample.csv file specified in the 
    config.
    """
    # Get the source, must be 'local' or 'SRA'
    try: 
        return SAMPLE_DF.loc[(SAMPLE_DF.Run == name) & (SAMPLE_DF.Library == int(library))].Source.item()
    except ValueError:
        print(f"The sample: {name} has more than one value for source.")


def get_layout(name, library):
    """
    Use the sample name and library to determine what the 
    library layout is ("PAIRED" or "SINGLE"). 

    Get the information from the samples dataframe 
    specified in the config file.
    """
    # Get the layout, must be 'PAIRED' or 'SINGLE'
    try: 
        return SAMPLE_DF.loc[(SAMPLE_DF.Run == name) & (SAMPLE_DF.Library == int(library))].LibraryLayout.item()
    except ValueError:
        print(f"The sample: {name} has more than one value for LibraryLayout.")


def download_local(name, library, layout):
    """
    If the samples are in the local file system, 
    get them, zip them, and finish.
    """
    # Path to samples on filesystem is a string of a list.
    path = SAMPLE_DF.loc[(SAMPLE_DF.Run == name) & (SAMPLE_DF.Library == int(library))].Path.item()
    path = path.strip('][').split(', ')
    if layout == "PAIRED": 
        assert len(path) == 2, "There is only one file for paired reads, is it interleaved?"
        # Sometimes the files are zipped and othertimes not.
        if path[0].endswith(".gz"):
            shell("cp {path[0]} {snakemake.output[0]}; cp {path[1]} {snakemake.output[1]}")  
        else:
            shell("gzip -c {path[0]} > {snakemake.output[0]}; gzip -c {path[1]} > {snakemake.output[1]}")
    elif layout == "SINGLE": 
        assert len(path) == 1, "There shouldn't be more than one file for single-ended reads."
        # Check if runs are zipped
        if path[0].endswith(".gz"):
            shell("cp {path[0]} {snakemake.output[0]}")  
        else:
            shell("gzip -c {path[0]} > {snakemake.output[0]}")
    else: 
        assert False, "The layout must be either PAIRED or SINGLE" 


# TODO: Will only work for paired files at present. Needs updating.
def download_sra(name, library, attempts = 3):
    """
    If the samples are on the SRA, download them, 
    try if it fails up the specified number of attempts.
    Output is zipped fastq file.
    """
    while attempt < attempts:
        try: 
            shell("fasterq-dump {snakemake.wildcards.accession} --outdir {snakemake.params} --temp {snakemake.params} --threads {snakemake.threads} -f; gzip -c {snakemake.params}/*_1.fastq > {snakemake.output[0]}; gzip -c {snakemake.params}/*_2.fastq > {snakemake.output[1]}")  
            break 
        except:
            attempt += 1


def main():
    """
    Main function to run the script.
    """
    # The sample name should refer to the accession attribute
    name = snakemake.wildcards.accession
    # The library should refer to when it was prepared (will be merged downstream)
    library = snakemake.wildcards.library
    # Get the source for this run (i.e., SRA or local)
    source = get_source(name, library)
    # Get the library layout (paried or single ended)
    layout = get_layout(name, library)
    # Download the files
    if source == 'local':
        download_local(name, library, layout)
    elif source == 'SRA':
        download_sra(name, library)
    else:
        assert False, "Source must be either local or SRA."


if __name__ == '__main__':
    main()

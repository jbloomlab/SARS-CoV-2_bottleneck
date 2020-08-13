__author__ = "Will Hannon"
__copyright__ = "Copyright 2020, Will Hannon"
__email__ = "wwh22@uw.edu"
__license__ = "MIT"

import os
import pandas as pd
from snakemake.shell import shell

# Determine if the sample is local or public/remote
source = pd.read_csv(snakemake.config['samples']['file']).set_index('Run').at[snakemake.wildcards.accession, 'Source']
shell("touch {source}")
# If the sample needs to be downloaded from the NCBI SRA
if source == "public":
    # Number of the attempt if the download fails 
    attempt = 0
    # Attempt at three times
    while attempt < 3:
    # Try the download
        try:
            shell("fasterq-dump {snakemake.wildcards.accession} --outdir {snakemake.params.outdir} --temp {snakemake.params.outdir} --threads {snakemake.threads} -f &> {snakemake.log}")  
            break              
        # If the download fails increment the attempt and try again
        except: 
            attempt += 1

# If the sample exists somewhere locally
elif source == "local":
    # Get the path to the files
    path = pd.read_csv(snakemake.config['samples']['file']).set_index('Run').at[snakemake.wildcards.accession, 'Path'].strip('][').split(', ')
    # Check if the run is paired or single-ended. 
    layout = pd.read_csv(snakemake.config['samples']['file']).set_index('Run').at[snakemake.wildcards.accession, 'LibraryLayout']
    # Two approaches depending on the layout
    if layout == "PAIRED":
        # Make sure there are two paths 
        assert len(path) == 2
        # Check if runs are zipped
        if path[0].endswith(".gz"):
            # Make directory and move into the correct direcotry with naming scheme 
            shell("cp {path[0]} {snakemake.output[0]}.gz; cp {path[1]} {snakemake.output[1]}.gz") 
            # Unzip the *.gz files
            shell("gunzip {snakemake.output[0]}.gz; gunzip {snakemake.output[1]}.gz")
        else:
            # Make directory and move into the correct direcotry with naming scheme 
            shell("cp {path[0]} {snakemake.output[0]}; cp {path[1]} {snakemake.output[1]}")  
    else: 
        # Make sure there is only one path
        assert len(path) == 1
        # Check if runs are zipped
        if path[0].endswith(".gz"):
            # Make directory and move into the correct direcotry with naming scheme 
            shell("cp {path[0]} {snakemake.output[0]}.gz")
            # USnzip the *.gz files
            shell("gunzip {snakemake.output[0]}.gz")
        else:
            # Make directory and move into the correct direcotry with naming scheme 
            shell("cp {path[0]} {snakemake.output}")  

### ======= Python utilities for running the pipeline ======= ###
#
# Author: Will Hannon 
# Email: wwh22@uw.edu
# Date: 10/30/2020
#

def single_ended(Run):
    """
    Returns `True` if the Library Layout is single-end reads.
    """
    # Get the Library Layout
    layout = pd.read_csv(config['samples']['file']).set_index('Run').at[Run, 'LibraryLayout']

    if not isinstance(layout, str):
        layout = layout.to_list()[0]
    # Check the layout
    if layout == "SINGLE":
        return True
    else:
        return False


def aggregate_csv(csv_list, out_csv):
    """
    This function aggregates lists of CSVs from
    the input of a snakemake rule. It's intended
    to be used in any aggregation rules throughout 
    the pipeline. 

    Try/Except to deal with empty CSVs.
    """
    valid_csv_list = []
    for csv in csv_list:
        try:
            pd.read_csv(csv)
            valid_csv_list.append(csv)
        except: 
            pass
    df = pd.concat(map(pd.read_csv, valid_csv_list))
    df.to_csv(out_csv[0], index = False)


def get_avaliable_fastqs(wildcards):
    """
    This function fills in the avaliable fastqs depending on 
    the library layout for the rule `trim_adapters_(pe | se)`. 
    """
    # If the layout is single-ended.
    if single_ended(wildcards.accession):
        # Return the target files.
        return expand(join(config["fastq_dir"], "{accession}-{library}", "{accession}-{library}.fastq.gz"), accession=wildcards.accession, library=wildcards.library)
    # Otherwise the layout is assumed to be paired-ended. 
    return expand([join(config["fastq_dir"], "{accession}-{library}", "{accession}-{library}_1.fastq.gz"),
                   join(config["fastq_dir"], "{accession}-{library}", "{accession}-{library}_2.fastq.gz")], accession=wildcards.accession, library=wildcards.library)


def get_avaliable_trimmed_fastqs(wildcards):
    """
    This function fills in the avaliable trimmed fastqs 
    depending on the library layout for the rule `filter_reads_(pe | se)`. 
    """
    # If the layout is single-ended.
    if single_ended(wildcards.accession):
        # Return the target files.
        return expand(join(config["trim_dir"], "{accession}-{library}", "{accession}-{library}.trimmed.fastq.gz"), accession=wildcards.accession, library=wildcards.library)
    # Otherwise the layout is assumed to be paired-ended. 
    return expand([join(config["trim_dir"], "{accession}-{library}", "{accession}-{library}_1.trimmed.fastq.gz"),
                   join(config["trim_dir"], "{accession}-{library}", "{accession}-{library}_2.trimmed.fastq.gz")], accession=wildcards.accession, library=wildcards.library)
     

def get_avaliable_filtered_fastqs(wildcards):
    """
    This function fills in the avaliable filtered fastq's 
    depending on the library layout. 
    """
    # If the layout is single-ended.
    if single_ended(wildcards.accession):
        # Return the target files.
        return expand(join(config["filter_dir"], "{accession}-{library}", "{accession}-{library}.filtered.fastq.gz"), accession=wildcards.accession, library=wildcards.library)
    # Otherwise the layout is assumed to be paired-ended. 
    return expand([join(config["filter_dir"], "{accession}-{library}", "{accession}-{library}_1.filtered.fastq.gz"),
                   join(config["filter_dir"], "{accession}-{library}", "{accession}-{library}_2.filtered.fastq.gz")], accession=wildcards.accession, library=wildcards.library)
     

# TODO: Will break if the same accession needs multiple viral genomes
def get_genome(wildcards, index = "samtools"):
    """ Function to get the correct genome for variant calling, indexed with samtools. 
    """
    # Read the metadata into Pandas dataframe
    samples_df = pd.read_csv(config['samples']['file'])
    # Get the viral genome for a given accession
    virusname = samples_df.loc[samples_df.Run == wildcards.accession, ["Virus"]].values.flatten().tolist()[0]
    # Return the viral genome file
    if index == "samtools": 
        return expand(join(config['index_dir']['samtools'], '{virus}.fa'), virus=virusname)
    elif index == "BWA":
        return expand(join(config['index_dir']['bwa'], '{virus}.fa'), virus=virusname)


def get_organism(wildcards): 
    """ Determine if organism is host or virus depending on wildcard.
    """
    # Read in the sample dataframe
    samples_df = pd.read_csv(config['samples']['file'])
    # Get the host and virus names as a list
    virus = samples_df.Virus.tolist()
    host = samples_df.Host.tolist()
    # Figure out if the organism is host or virus 
    if wildcards.organism in host:
        # Invert grep
        return "-v"
    # Don't invert grep
    return ""


def get_avaliable_bams(wildcards):
    """ 
    This function determines which filtered bam files to generate
    based on the accession.
    """
    # Read the metadata into Pandas dataframe
    samples_df = pd.read_csv(config['samples']['file'])
    # Get the viral genome for a given accession
    library = np.unique(samples_df.loc[samples_df.Run == wildcards.accession, ["Library"]].values.flatten()).tolist()

    return expand(join(config['align_dir'], "{{aligner}}", "{{accession}}-{library}", "{{accession}}-{library}.{{aligner}}.sorted.bam"), library=library)




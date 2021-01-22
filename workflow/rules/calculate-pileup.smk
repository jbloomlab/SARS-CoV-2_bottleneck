### ======= Identify variants by calculating and parsing a pileup file ======= ###
#
# Author: Will Hannon 
# Email: wwh22@uw.edu
# Date: 11/02/2020
#

localrules: aggregate_pileup

rule get_pileup:
    """ Calculate pileup statistics and process with python/pysam.
    """
    input: bam=join(config['align_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.sorted.marked.merged.bam"),
           bai=join(config['align_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.sorted.marked.merged.bam.bai"),        
           genome=get_genome
    output: join(config['pileup_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.pysam.pileup.csv")
    params: score=config['BQ']
    conda: "../envs/pysam.yml"
    script: "../scripts/pysam_pileup.py"

    
rule aggregate_pileup:
    """ Aggregate variants called with python and pysam.. 
    """
    input: expand(join(config['pileup_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.pysam.pileup.csv"), accession=pd.read_csv(config['samples']['file'])['Run'], aligner=['BWA'])
    output: join(config['pileup_dir'], "pysam-variants.csv")
    run: aggregate_csv(input, output)


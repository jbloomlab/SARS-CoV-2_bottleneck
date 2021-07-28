### ======= Generate consensus sequences from aligned viral reads. ======= ###
#
# Author: Will Hannon 
# Email: wwh22@uw.edu
# Date: 01/15/2021
#

rule make_consensus:
    input: bam=expand(join(config['align_dir'], "{aligner}", "{accession}", "{accession}.{spid}_{replicate}.sorted.merged.bam"), replicate = [1,2]), 
           bai=expand(join(config['align_dir'], "{aligner}", "{accession}", "{accession}.{spid}_{replicate}.sorted.merged.bam.bai"), replicate = [1,2]), 
           genome = lambda wildcards: get_genome(wildcards) 
    output: join(config['consensus_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.consensus.fa")
    conda: "../envs/pysam.yml"
    script: "../scripts/make-consensus-sequence.py"


rule join_consensus: 
    """
    Join the consensus sequences from bcftools 
    for multiple sequence alignment with MAFFT.
    """
    input: expand(join(config['consensus_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.consensus.fa"), accession=list(set(pd.read_csv(config['samples']['file'])['Run'])), aligner=['BWA'])
    output: join(config['consensus_dir'], "merged.consenus.fa")
    shell: "cat {input} > {output}"


rule align_consensus:
    """
    Perform multiple sequence alignment with 
    MAFFT. 
    """
    input: join(config['consensus_dir'], "merged.consenus.fa")
    output: join(config['consensus_dir'], "aligned_consensus.fa")
    conda: '../envs/consensus.yml'
    shell: "mafft {input} > {output}"


### ======= Generate consensus sequences from aligned viral reads. ======= ###
#
# Author: Will Hannon 
# Email: wwh22@uw.edu
# Date: 01/15/2021
#

rule make_consensus:
    input: bam=expand(join(config['align_dir'], "{{aligner}}", "{{spid}}_{replicate}", "{{spid}}_{replicate}.{{aligner}}.sorted.merged.bam"), replicate = [1,2]), 
           bai=expand(join(config['align_dir'], "{{aligner}}", "{{spid}}_{replicate}", "{{spid}}_{replicate}.{{aligner}}.sorted.merged.bam.bai"), replicate = [1,2]), 
           genome = join(config['index_dir']['samtools'], "SARS2.fa") 
    output: join(config['consensus_dir'], "{aligner}", "{spid}", "{spid}.{aligner}.consensus.fa")
    conda: "../envs/pysam.yml"
    script: "../scripts/make-consensus-sequence.py"


rule join_consensus: 
    """
    Join the consensus sequences from bcftools 
    for multiple sequence alignment with MAFFT.
    """
    input: expand(join(config['consensus_dir'], "{aligner}", "{spid}", "{spid}.{aligner}.consensus.fa"), spid=list(set(pd.read_csv(config['samples']['file'])['spID'])), aligner=['BWA'])
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


rule iqtree_boat_phylogeny:
    """Make a phylogeny of the high-quality samples from the boat."""
    input: join(config['consensus_dir'], "aligned_consensus.fa")
    output: join(config['consensus_dir'], "aligned_consensus.fa.iqtree")
    conda: "../envs/iqtree.yml"
    shell: "iqtree -s {input}"

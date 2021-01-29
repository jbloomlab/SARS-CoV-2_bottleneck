### ======= Generate consensus sequences and make a phylogeny. ======= ###
#
# Author: Will Hannon 
# Email: wwh22@uw.edu
# Date: 01/29/2021
#

rule make_consensus:
    input: bam=expand(join(config['align_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.sorted.marked.merged.bam"), accession=list(set(pd.read_csv(config['samples']['file'])['Run'])), aligner=['BWA']),
           bai=expand(join(config['align_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.sorted.marked.merged.bam.bai"), accession=list(set(pd.read_csv(config['samples']['file'])['Run'])), aligner=['BWA']), 
           genome=join(config['index_dir']['samtools'], 'SARS2.fa')
    output: join(config['consensus_dir'], "merged.consenus.fa")
    conda: "../envs/pysam.yml"
    script: "../scripts/make-consensus-sequence.py"

rule iqtree_phylogeny:
    input: join(config['consensus_dir'], "merged.consenus.fa")
    output: join(config['consensus_dir'], "merged.consenus.fa.iqtree")
    conda: "../envs/iqtree.yml"
    shell: "iqtree -s {input}"

### ======= Generate consensus sequences from aligned viral reads. ======= ###
#
# Author: Will Hannon 
# Email: wwh22@uw.edu
# Date: 01/15/2021
#

rule calculate_pileup:
    """ 
    Calculate the pileup of bases at every position in virus genome.
    Only calculates bases with Phred scaled quality score higher than 30.
    """
    input: bam=join(config['align_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.sorted.marked.merged.bam"),
           bai=join(config['align_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.sorted.marked.merged.bam.bai"),        
           genome=get_genome
    output: join(config['pileup_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.mpileup.txt")
    params: score=config['BQ']
    conda: '../envs/samtools.yml'
    shell: "samtools mpileup -d 0 -E --excl-flags UNMAP,SECONDARY,QCFAIL -q {params.score} -Q {params.score} -f {input.genome} {input.bam} -O -s --reverse-del -a -o {output}"


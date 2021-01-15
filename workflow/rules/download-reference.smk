### ======= Download and format reference genomes and annotations ======= ###
#
# Author: Will Hannon 
# Email: wwh22@uw.edu
# Date: 10/30/2020
#

rule get_ref:
    """ Download the reference genome fasta.
    """
    output: join(config['ref_dir'], '{genome}.fa')
    params: ftp = lambda wildcards: config[wildcards.genome]['ref']
    wildcard_constraints: genome="[^.]+"
    shell: 'wget -O - {params.ftp} | gunzip -c > {output}'


rule get_gtf:
    """ Download the gtf format annotation. 
    """
    output: join(config['gtf_dir'], '{genome}.gtf')
    params: ftp = lambda wildcards: config[wildcards.genome]['gtf']
    wildcard_constraints: genome="[^.]+"
    shell: 'wget -O - {params.ftp} | gunzip -c > {output}'



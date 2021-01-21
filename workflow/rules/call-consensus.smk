### ======= Generate consensus sequences from aligned viral reads. ======= ###
#
# Author: Will Hannon 
# Email: wwh22@uw.edu
# Date: 01/15/2021
#

# TODO: Add in a rule to aggregate the filtered fastqs. 

rule tadpole_error_correct: 
    """
    Using `tadpole.sh` from the BBmap suite,
    correct errors in the reads at the ends
    and in the middle of the reads. Use a kmer
    that at least less than 1/2 the length of the 
    read, ideally 1/3 read length. 
    """
    input: get_avaliable_filtered_fastqs
    output: join(config["consensus_dir"], "{accession}-{library}", "{accession}-{library}.ec.fastq.gz")
    threads: config['threads']['max_cpu']
    params: kmer = config['kmer']['error_correct']
    conda: '../envs/filter.yml'
    shell: "tadpole.sh in={input.reads} out={output} mode=correct k={params.kmer}"


rule tadpole_extend_reads:
    """
    Using `tadpole.sh` from the BBmap suite, 
    extend the error-corrected reads to longer
    reads that can be assembled into contigs. 
    """
    input: join(config["consensus_dir"], "{accession}-{library}", "{accession}-{library}.ec.fastq.gz")
    output: join(config["consensus_dir"], "{accession}-{library}", "{accession}-{library}.extended.fastq.gz")
    threads: config['threads']['max_cpu']
    params: kmer = config['kmer']['extend_reads']
    conda: '../envs/filter.yml'
    shell: "tadpole.sh in={input} out={output} mode=extend k={params.kmer} el=50 er=50"


rule tadpole_assemble_reads:
    """
    Using `tadpole.sh` from the BBmap suite, 
    assemble the error corrected and extended
    reads into longer contigs which can then 
    be polished into a full assembly.
    """
    input: join(config["consensus_dir"], "{accession}-{library}", "{accession}-{library}.extended.fastq.gz")
    output: join(config["consensus_dir"], "{accession}-{library}", "{accession}-{library}.contigs.fa")
    threads: config['threads']['max_cpu']
    params: kmer = config['kmer']['assemble_contigs']
    conda: '../envs/filter.yml'
    shell: "tadpole.sh in={input} out={output} k={params.kmer}"

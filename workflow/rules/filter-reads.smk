### ======= Filter processed fastq files by kmers with BBduk ======= ###
#
# Author: Will Hannon 
# Email: wwh22@uw.edu
# Date: 10/30/2020
#

rule filter_reads_se:
    """ Filter out viral reads using kmer filtering w/ BBduk Single-Ended.
    """
    input: 
        reads=get_avaliable_trimmed_fastqs,
        genome=get_BWA_ref_genome
    output: 
        matched=join(config['filter_dir'], "{accession}", "{accession}.filtered.fastq.gz"),
        unmatched=join(config['filter_dir'], "{accession}", "{accession}.unfiltered.fastq.gz"),
        stats=join(config['filter_dir'], "{accession}", "{accession}.filter.stats.txt")
    threads: config['threads']['max_cpu']
    conda: '../envs/filter.yml'
    shell:
        """
        bbduk.sh -Xmx80g \
            in={input.reads} \
            out={output.unmatched} \
            outm={output.matched} \
            ref={input.genome} \
            k=31 \
            hdist=2 \
            stats={output.stats} \
            overwrite=TRUE \
            t={threads}
        """

rule filter_reads_pe:
    """ Filter out viral reads using kmer filtering w/ BBduk Paired-end.
    """
    input: 
        reads=get_avaliable_trimmed_fastqs,
        genome=get_BWA_ref_genome
    output: 
        matched=[join(config['filter_dir'], "{accession}", "{accession}_1.filtered.fastq.gz"), 
                 join(config['filter_dir'], "{accession}", "{accession}_2.filtered.fastq.gz")],
        unmatched=[join(config['filter_dir'], "{accession}", "{accession}_1.unfiltered.fastq.gz"), 
                   join(config['filter_dir'], "{accession}", "{accession}_2.unfiltered.fastq.gz")],
        stats=join(config['filter_dir'], "{accession}", "{accession}.filter.stats.txt")
    threads: config['threads']['max_cpu']
    conda: '../envs/filter.yml'
    shell:
        """
        bbduk.sh -Xmx80g \
            in1={input.reads[0]} \
            in2={input.reads[1]} \
            out1={output.unmatched[0]} \
            out2={output.unmatched[1]} \
            outm1={output.matched[0]} \
            outm2={output.matched[1]} \
            ref={input.genome} \
            k=31 \
            hdist=2 \
            stats={output.stats} \
            overwrite=TRUE \
            t={threads}
        """
### ======= Short read alignment of filtered reads to normal or hybrid genome ======= ###
#
# Author: Will Hannon 
# Email: wwh22@uw.edu
# Date: 11/02/2020
#

rule bwa_align:
    """ 
    Perform short read alignment of the filtered 
    viral readas with `bwa-mem`.
    
    Sort the aligned reads with samtools sort.
    """
    input: 
        reads = get_avaliable_filtered_fastqs,
        genome = lambda wildcards: get_genome(wildcards, index = "BWA") 
    output: join(config['align_dir'], "BWA", "{accession}-{library}", "{accession}-{library}.BWA.sorted.bam")
    threads: config['threads']['max_cpu']
    params: tmp=join(config['align_dir'], "BWA", "{accession}-{library}", "{accession}-{library}.final.tmp")
    conda: '../envs/align.yml'
    shell: 
        """
        bwa mem -t {threads} \
            {input.genome} \
            {input.reads} | \
            samtools sort -o {output} -T {params.tmp} - 
        """


rule mark_duplicates:
    """ 
    This rule uses `Picard` to mark/remove probable PCR duplicates.
    Only removes duplicates if `true` is specified in the config file, 
    otherwise this tool only marks them.
    """
    input: bam=join(config['align_dir'], "{aligner}", "{accession}-{library}", "{accession}-{library}.{aligner}.sorted.bam")
    output: bam=join(config['align_dir'], "{aligner}", "{accession}-{library}", "{accession}-{library}.{aligner}.sorted.marked.bam"),
            metrics=join(config['qc_dir'], "{accession}-{library}", "Picard", "{accession}-{library}.{aligner}.metrics.txt")
    params: remove_duplicates=config['remove_duplicates']
    conda: '../envs/picard.yml'
    shell:
        """
        picard MarkDuplicates REMOVE_DUPLICATES={params.remove_duplicates} INPUT={input.bam} \
            OUTPUT={output.bam} METRICS_FILE={output.metrics} 
        """


rule merge_bam:
    """
    Merge the BAM files from the same sample but different runs. 
    """
    input: get_avaliable_bams
    output: join(config['align_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.sorted.marked.merged.bam")
    conda: '../envs/samtools.yml'
    shell: "samtools merge {output} {input}" 


rule index_bam:
    """ Index mapped, sorted, merged, and marked bam files. 
    """
    input: join(config['align_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.sorted.marked.merged.bam")
    output: join(config['align_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.sorted.marked.merged.bam.bai")
    conda: '../envs/samtools.yml'
    shell: "samtools index {input}"





    

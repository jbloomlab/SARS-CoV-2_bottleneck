### ======= Short read alignment ======= ###

rule bwa_mem:
    """ Perform mapping with `bwa-mem`.
    """
    input: 
        reads=get_avaliable_trimmed_fastqs,
        genome=get_BWA_ref_genome
    output: join(config['align_dir'], "BWA", "{accession}", "{accession}.BWA.bam")
    threads: config['max_cpu']
    params: tmp=join(config['align_dir'], "BWA", "{accession}", "{accession}.final.tmp")
    log: "logs/bwa_mem/{accession}.log"
    conda: '../envs/align.yml'
    shell: 
        """
        bwa mem -t {threads} \
            {input.genome} \
            {input.reads} | \
            samtools sort -o {output} -T {params.tmp} - &> {log}
        """

rule star_align:
    """ Perform mapping with `STAR` against host/viral genomes.
    """
    input: 
        reads=get_avaliable_trimmed_fastqs,
        genome=get_STAR_ref_genome
    output: join(config['align_dir'], "STAR", "{accession}", "{accession}.STAR.bam")
    params:
        prefix=join(config['align_dir'], "STAR", "{accession}", "{accession}.STAR."),
        intermediate_bam=join(config['align_dir'], "STAR", "{accession}", "{accession}.STAR.Aligned.sortedByCoord.out.bam"),
        logtmp=join(config['align_dir'], "STAR", "{accession}", "{accession}.STAR.Log.final.out")
    log: join(config['qc_dir'], "{accession}", "STAR", "{accession}.STAR.Log.final.out")
    conda: '../envs/align.yml'    
    threads: config['max_cpu']
    shell:
        """
        STAR --genomeDir {input.genome} \
            --runThreadN {threads} \
            --readFilesIn {input.reads} \
            --outSAMtype BAM SortedByCoordinate \
            --twopassMode Basic \
            --limitBAMsortRAM 100000000000 \
            --outFileNamePrefix {params.prefix} \
            --outSAMunmapped Within KeepPairs

        # Rename bam files to target
        mv {params.intermediate_bam} {output}

        # Move log files
        mv {params.logtmp} {log} 
        """


rule samtools_filter:
    """ Filter reads that mapped to either virus or host from `STAR` and `BWA`. 
    """
    input: join(config['align_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.bam"),
    output: bam=join(config['split_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.{organism}.bam"),
            index=join(config['split_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.{organism}.bam.bai")
    params: 
        invert=lambda wildcards: "true" if (wildcards.organism != "virus") else "false",
        contig=lambda wildcards: config[pd.read_csv(config['samples']['file']).set_index('Run').at[wildcards.accession, 'Virus']]['contig']
    conda: '../envs/samtools.yml'
    shell:
        """
        samtools view -H {input} > {output}
        
        if {params.invert}; then

            samtools view {input} \
                | awk '$3!= "{params.contig}"' \
                >> {output.bam} 

            samtools sort {output.bam} -o {output.bam}
            samtools index {output.bam}

        else

            samtools view {input} \
                | awk '$3== "{params.contig}"' \
                >> {output.bam} 

            samtools sort {output.bam} -o {output.bam}
            samtools index {output.bam}

        fi
        """

     
rule remove_duplicates:
    """ This rule uses `Picard` to mark/remove probable PCR duplicates.
    """
    input: 
        bam=join(config['split_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.{organism}.bam")
    output:
        bam=join(config['dedup_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.{organism}.bam"),
        metrics=join(config['dedup_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.{organism}.metrics.txt")
    conda: '../envs/picard.yml'
    shell:
        """
        picard MarkDuplicates REMOVE_DUPLICATES=true INPUT={input.bam} \
            OUTPUT={output.bam} METRICS_FILE={output.metrics} 
        """

rule index_bam:
    """ Sort and index mapped reads from `STAR` and `BWA`. 
    """
    input: join(config['dedup_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.{organism}.bam"),
    output: join(config['dedup_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.{organism}.bam.bai")
    conda: '../envs/samtools.yml'
    shell:
        """
        samtools sort {input} -o {input}
        samtools index {input} 
        """

##========##

rule bedtools_coverage:
    """ 
    """
    input: bam = join(config['split_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.virus.bam"),
           bai = join(config['split_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.virus.bam.bai"),
           genome = join(config['index_dir']['samtools'], "SARS2.fa.fai")
    output: join(config['split_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.bedgraph")
    params: binsize = 50
    conda: '../envs/samtools.yml'
    shell:
        """
        awk 'FS=OFS="\t"{{print $1, 0, $2}}' {input.genome} \
            | bedtools makewindows -b - -w {params.binsize} \
            | bedtools coverage -a stdin -b {input.bam} -mean \
            > {output}
            
        sed -i "s/$/\t{wildcards.accession}/" {output} 
                       
        """


rule merge_coverage:
    """ 
    """
    input: expand(join(config['split_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.bedgraph"),accession=pd.read_csv(config['samples']['file'])['Run'], aligner=['BWA'])
    output: join(config['split_dir'], "merged.bedgraph")
    params: header=temp(join(config['split_dir'], "merged.bedgraph.tmp"))
    shell:
        """
        cat {input} > {params.header}
        awk 'BEGIN{{print "Name\tStart\tStop\tDepth\tAccession"}}1' {params.header} > {output}
        """


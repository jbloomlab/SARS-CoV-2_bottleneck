### ======= Quality Control and Reporting ======= ###

rule fastqc: 
    """ Generate a QC report for unfiltered reads. 
    """
    input: get_avaliable_fastqs
    output: directory(join(config['qc_dir'], '{accession}/fastqc'))
    conda: '../envs/qc.yml'
    log: "logs/fastqc/{accession}.log"
    shell: "mkdir -p {output}; fastqc {input} --outdir {output} &> {log}"


rule samtools_stats:
    """ Calculate bam stats with samtools
    """
    input: 
        join(config['dedup_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.virus.bam"),
        join(config['dedup_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.virus.bam.bai")
    output: join(config['qc_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.virus.bam.stats")
    log: "logs/samtools_stats/{accession}.{aligner}.log"
    conda: '../envs/samtools.yml'
    shell: "samtools stats {input[0]} > {output} &> {log}"


rule multiqc:
    """
    Collate all QC reports for all samples into a
    single report. 
    """
    input: 
        qc=expand([join(config['qc_dir'], '{accession}/fastqc'),
                   join(config['qc_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.virus.bam.stats"),
                   join(config['qc_dir'], "{accession}", "STAR", "{accession}.STAR.Log.final.out")],
                   accession=pd.read_csv(config['samples']['file'])['Run'], aligner=['BWA','STAR'])
    output: directory(join(config['qc_dir'], 'multiqc'))
    conda: '../envs/qc.yml'
    params: 
        basename="multiqc_report",
        indir=config['qc_dir']
    shell:
        """
        multiqc {params.indir} \
            -f -o {output} \
            -n {params.basename} 
        """

### ======= Quality control and reporting with Multiqc ======= ###
#
# Author: Will Hannon 
# Email: wwh22@uw.edu
# Date: 10/31/2020
#

rule fastqc_raw: 
    """ Generate a QC report for unfiltered reads. 
    """
    input: get_avaliable_fastqs
    output: directory(join(config['qc_dir'], "{accession}-{library}", "fastqc-raw"))
    conda: '../envs/qc.yml'
    shell: "mkdir -p {output}; fastqc {input} --outdir {output}"


rule fastqc_trim: 
    """ Generate a QC report for trimmed reads. 
    """
    input: get_avaliable_trimmed_fastqs 
    output: directory(join(config['qc_dir'], "{accession}-{library}", "fastqc-trim"))
    conda: '../envs/qc.yml'
    shell: "mkdir -p {output}; fastqc {input} --outdir {output}"


rule samtools_stats:
    """ Calculate bam stats with samtools.
    """
    input: join(config['align_dir'], "{aligner}", "{accession}-{library}", "{accession}-{library}.{aligner}.sorted.marked.bam"),
           join(config['align_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.sorted.marked.merged.bam")
    output: join(config['qc_dir'], "{accession}-{library}", "{aligner}", "{accession}-{library}.{aligner}.bam.stats"),
            join(config['qc_dir'], "{accession}", "{aligner}", "{accession}.{aligner}.merged.bam.stats")
    conda: '../envs/samtools.yml'
    shell: 
        """
        samtools stats {input[0]} > {output[0]}
        samtools stats {input[1]} > {output[1]}
        """

# TODO: fix the paths to reflect the correct inputs. This might require a dedicated function like those in `common.smk`
rule multiqc:
    """
    Collate all QC reports for all samples into a
    single report. 
    """
    input: 
        qc=expand([join(config['qc_dir'], "{accession}", "fastqc"),
                   join(config['qc_dir'], "{accession}", "{aligner}", "{accession}.{aligner}.sorted.bam.stats"),
                   join(config['qc_dir'], "{accession}", "{aligner}", "{accession}.{aligner}.virus.sorted.bam.stats"),                  
                   join(config['qc_dir'], "{accession}", "Picard", "{accession}.{aligner}.virus.metrics.txt")],
                   accession=pd.read_csv(config['samples']['file'])['Run'], aligner=['BWA'])
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

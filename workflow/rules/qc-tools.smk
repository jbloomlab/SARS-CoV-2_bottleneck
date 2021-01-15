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


rule samtools_unmerged_stats:
    """ 
    Calculate basic statistics about the unmerged BAM files.
    """
    input: join(config['align_dir'], "{aligner}", "{accession}-{library}", "{accession}-{library}.{aligner}.sorted.marked.bam"),
    output: join(config['qc_dir'], "{accession}-{library}", "{aligner}", "{accession}-{library}.{aligner}.bam.stats"),
    conda: '../envs/samtools.yml'
    shell: 
        """
        samtools stats {input} > {output}
        """


rule samtools_merged_stats:
    """ 
    Calculate basic statistics about the merged BAM files.
    """
    input: join(config['align_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.sorted.marked.merged.bam")
    output: join(config['qc_dir'], "{accession}", "{aligner}", "{accession}.{aligner}.merged.bam.stats")
    conda: '../envs/samtools.yml'
    shell: 
        """
        samtools stats {input} > {output}
        """


def get_multiqc_targets():
    """
    This function generates all of the target files
    as input to the multiqc rule. This needs a 
    dedicated function because there are an uneven 
    number of extra runs per sample.
    """
    samples_df = pd.read_csv(config['samples']['file'])
    input_set = set()
    for index, row in samples_df.iterrows():
        # Fastqc for untrimmed and unfiltered library
        input_set.add(join(config['qc_dir'], f"{row['Run']}-{row['Library']}", "fastqc-raw"))
        # Fastqc post-trimming and basic quality filtering
        input_set.add(join(config['qc_dir'], f"{row['Run']}-{row['Library']}", "fastqc-trim"))
        # Trimmed read reports from fastp
        input_set.add(join(config['qc_dir'], f"{row['Run']}-{row['Library']}", "fastp", f"{row['Run']}-{row['Library']}.fastp.se.json"))
        # Filtered reads statistics from BBduk
        input_set.add(join(config['qc_dir'], f"{row['Run']}-{row['Library']}", "BBduk", f"{row['Run']}-{row['Library']}.filter.se.stats"))
        # Samtools stats for unmerged alignment
        input_set.add(join(config['qc_dir'], f"{row['Run']}-{row['Library']}", "BWA", f"{row['Run']}-{row['Library']}.BWA.bam.stats"))
        # Samtools stats for merged alignment
        input_set.add(join(config['qc_dir'], f"{row['Run']}", "BWA", f"{row['Run']}.BWA.merged.bam.stats"))
    
    return list(input_set)


rule multiqc:
    """
    Gather all relevant library and alignment stats into a single place
    that can be incorportated into the Snakemake report.
    """
    input: qc = get_multiqc_targets()
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

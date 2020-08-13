### ======= Import runs from the NCBI SRA ======= ###

rule get_fastq_se:
    """ Single-end implementation of `sra-toolkit`. 
    """
    output: join(config['sra_dir'], "{accession}", "{accession}.fastq")
    params: outdir=join(config['sra_dir'], "{accession}")
    log: "logs/sra_download/{accession}.log"
    threads: 8
    conda: '../envs/sra.yml'
    script: '../scripts/sra_download.py'

rule get_fastq_pe:
    """
    Single-end implementation of `sra-toolkit`. 
    """
    output: 
        join(config['sra_dir'], "{accession}", "{accession}_1.fastq"),
        join(config['sra_dir'], "{accession}", "{accession}_2.fastq")
    params: outdir=join(config['sra_dir'], "{accession}")
    log: "logs/sra_download/{accession}.log"
    threads: 8
    conda: '../envs/sra.yml'
    script: '../scripts/sra_download.py'


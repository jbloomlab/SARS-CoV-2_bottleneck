### ======= Index reference genomes in specified format ======= ###
#
# Author: Will Hannon 
# Email: wwh22@uw.edu
# Date: 10/30/2020
#

rule bwa_index:
    """ Index the genome with `BWA` before mapping.
    """
    input: join(config['ref_dir'], '{genome}.fa')
    output: join(config['index_dir']['bwa'], '{genome}.fa')
    conda: '../envs/align.yml'
    params: algorithm="bwtsw"
    shell:
        """
        cp {input} {output}
        bwa index -a {params.algorithm} {output}
        """

rule star_index:
    """ Index the virus & host genome with `STAR` before mapping.
    """
    input: 
        ref=join(config['ref_dir'], '{genome}.fa'),
        gtf=join(config['gtf_dir'], '{genome}.gtf')
    output: directory(join(config['index_dir']['star'], "{genome}"))
    params: 
        Nbases=SAindexNbases,
        sjdb=sjdbGTFfile
    threads: config['max_cpu']
    conda: '../envs/align.yml'
    shell:
        """
        # Make the output directory specified in the config file.
        mkdir -p {output}

        # run STAR generate genome command
        STAR --runThreadN {threads} \
            --runMode genomeGenerate \
            --genomeDir {output} \
            --genomeFastaFiles {input.ref} \
            --genomeSAindexNbases {params.Nbases} \
            --limitGenomeGenerateRAM 180000000000 \
            {params.sjdb} 
        """

rule samtools_index:
    """ Index genome with `samtools` for `BSQR`. 
    """
    input: join(config['ref_dir'], '{genome}.fa')
    output: 
        fa=join(config['index_dir']['samtools'], '{genome}.fa'),
        idx=join(config['index_dir']['samtools'], '{genome}.fa.fai'),
        idxdict=join(config['index_dir']['samtools'], '{genome}.fa.dict')
    conda: '../envs/samtools.yml'
    shell:
        """
        cp {input} {output.fa}
        samtools faidx {output.fa}
        samtools dict -o {output.idxdict} {output.fa}
        """


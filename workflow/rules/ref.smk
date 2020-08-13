### ======= Assemble reference genome ======= ###

rule get_ref:
    """ Download a reference genome.
    """
    output: join(config['ref_dir'], '{genome}.fa')
    params: ftp = lambda wildcards: config[wildcards.genome]['ref']
    wildcard_constraints: genome="[^.]+"
    shell: 'wget -O - {params.ftp} | gunzip -c > {output}'


rule get_gtf:
    """ Download a gtf annotation. 
    """
    output: join(config['gtf_dir'], '{genome}.gtf')
    params: ftp = lambda wildcards: config[wildcards.genome]['gtf']
    wildcard_constraints: genome="[^.]+"
    shell: 'wget -O - {params.ftp} | gunzip -c > {output}'


rule cat_ref:
    """ Concatenate virus/host genomes.
    """
    input: 
        virus_ref=join(config['ref_dir'], '{virus}.fa'),
        host_ref=join(config['ref_dir'], '{host}.fa')
    output: join(config['ref_dir'], '{virus}.{host}.fa')
    shell: "cat {input.virus_ref} {input.host_ref} > {output}"


rule cat_gtf:
    """ 
    Concatenate virus/host gtf files. This is 
    primarily for `STAR` which requires custom
    GTF files, so this is not downstream of
    `get_ref`. See run instructions for details. 
    """
    input: 
        virus_gtf=join(config['gtf_dir'], '{virus}.gtf'),
        host_gtf=join(config['gtf_dir'], '{host}.gtf')
    output: join(config['gtf_dir'], '{virus}.{host}.gtf')
    shell: "cat {input.virus_gtf} {input.host_gtf} > {output}"


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


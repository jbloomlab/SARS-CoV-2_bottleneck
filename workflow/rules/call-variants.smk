### ======= Identify variants using vatiant calling/annotation programs ======= ###
#
# Author: Will Hannon 
# Email: wwh22@uw.edu
# Date: 10/19/2020
#

rule samtools_mpileup:
    """ 
    Calculate the pileup of bases at every position in virus genome.
    """
    input: bam=join(config['align_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.sorted.merged.bam"),
           bai=join(config['align_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.sorted.merged.bam.bai"),        
           genome=get_genome
    output: join(config['pileup_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.mpileup.txt")
    params: score=config['BQ'],
            maxdepth=config['maxdepth']
    conda: '../envs/samtools.yml'
    shell: "samtools mpileup -d {params.maxdepth} -E --excl-flags UNMAP,SECONDARY,QCFAIL -q {params.score} -Q {params.score} -f {input.genome} {input.bam} -O -s --reverse-del -a -o {output}"


rule ivar_calling:
    """
    Variant calling with iVar. Benchmarked for viruses, seems like 
    a good alternative to Varscan.
    """
    input:
           bam=join(config['align_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.sorted.merged.bam"),
           bai=join(config['align_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.sorted.merged.bam.bai"),        
           gff=join(config['gff_dir'], "MeVChiTok.gff"),
           genome=get_genome
    output: txt=join(config['variant_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.ivar.ann.txt"),
            tsv=temp(join(config['variant_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.ivar.ann.tsv"))
    params: 
        quality=config['BQ'],
        maxdepth=config['maxdepth'],
        minimum_coverage=config['min_coverage'],
        prefix=join(config['variant_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.ivar.ann")
    conda: '../envs/variant.yml'
    threads: config['threads']['max_cpu']
    shell: 
        """
        samtools mpileup -aa -A --excl-flags UNMAP,SECONDARY,QCFAIL -d {params.maxdepth} -B -Q 0 {input.bam} | ivar variants -p {params.prefix} -q {params.quality} -t 0.01 -r {input.genome} -g {input.gff}
        cp {output.tsv} {output.txt}
        """


rule varscan_calling:
    """ SNP calling with Varscan. Parameters are controlled from the config file.  
    """
    input: 
        varscan=join(config['tools'], "VarScan.v2.4.0.jar"),
        pileup=join(config['pileup_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.mpileup.txt")
    output: 
        variants=join(config['variant_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.varscan.vcf"),
    params:
        minimum_coverage=config['min_coverage'],
        minumum_supporting_reads=config['min_reads_supporting'],
        minimum_base_quality=config['BQ'],
        minimum_variant_freq=config['min_allele_frequency'],
        strand_filter=config['strand_bias_filter']
    conda: '../envs/variant.yml'    
    shell:
        """
        # Call SNPs and InDels using the mpileup file
        java -jar {input.varscan} \
            mpileup2cns {input.pileup} \
            --variants 1 \
            --output-vcf 1 \
            --min-coverage {params.minimum_coverage} \
            --min-reads2 {params.minumum_supporting_reads} \
            --min-avg-qual {params.minimum_base_quality} \
            --strand-filter {params.strand_filter} \
            --min-var-freq {params.minimum_variant_freq} > {output.variants}
        """



rule lofreq_calling:
    """ Call variants (SNP and indels) with the reference. Using only defaut filtering.
    """
    input: bam=join(config['align_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.sorted.merged.bam"),
           bai=join(config['align_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.sorted.merged.bam.bai"),        
           genome=get_genome
    output: join(config['variant_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.lofreq.vcf")
    params: 
        maxdepth=config['maxdepth']
    conda: '../envs/variant.yml'
    threads: config['threads']['max_cpu']
    shell:
        """
        if [[ $(samtools view {input.bam} | head -n 5) ]]; then
            lofreq call-parallel --pp-threads {threads} \
                -f {input.genome} \
                -N \
                -d {params.maxdepth} \
                {input.bam} \
                -o {output}
        else
            touch {output}
        fi
        """ 


def get_virus(wildcards):
    """ Get the virus name from accession
    """
    # Get the sample dataframe
    sample_df = pd.read_csv(config['samples']['file'])
    # Get the correct wildcard
    genome = sample_df.loc[sample_df['Run'] == wildcards.accession, 'Virus'].iloc[0]
    # Generate the ouput file
    return genome


rule annotate_vcf:
    """
    Use the program SnpEff to annotate the effect of 
    mutations using the custom virus genome.
    """
    input:
        virusdir=lambda wildcards: expand(join(config['tools'], 'snpEff/data/{genome}'), genome=get_virus(wildcards)),
        vcf=join(config['variant_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.{caller}.vcf")
    output: join(config['variant_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.{caller}.ann.vcf")
    conda: '../envs/java.yml'
    params:
        snpEff=join(config['tools'], 'snpEff/snpEff.jar'),
        config=join(config['tools'], "snpEff/snpEff.config"),
        genome=get_virus
    shell:
        """
        java -jar {params.snpEff} -c {params.config} -noStats \
            -v {params.genome} {input.vcf} > \
            {output}
        """


rule vcf_to_table:
    """
    Convert varscan VCF files to tables for easy data
    analysis in R or Python. This will only work with 
    Varscan2 due to the program specific data fields. 
    """
    input: join(config['variant_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.{caller}.ann.vcf")
    output: join(config['variant_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.{caller}.ann.txt")
    conda: '../envs/variant.yml'    
    shell:
        """
        if [ -s {input} ]; then
            gatk VariantsToTable \
                -V {input} \
                -F CHROM -F POS -F QUAL -F REF -F ALT \
                -F DP -F AF -F FILTER -GF DP \
                -GF RD -GF FREQ -GF SDP -GF AD -F ANN \
                -O {output}
        else
            touch {output}
        fi
        """


rule add_metadata:
    """
    This rule adds in the metadata from the csv file
    that is used to run the experiment. 
    """
    input: join(config['variant_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.{caller}.ann.txt")
    output: join(config['variant_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.{caller}.ann.csv")
    params: metadata=config['samples']['file']
    conda: '../envs/r.yml'
    script: "../scripts/add_metadata.R"


rule aggregate_variants:
    """
    This rule aggregates all of the variants. 
    """
    input: expand([join(config['variant_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.{caller}.ann.csv")], accession=samples, aligner=['BWA'], caller=['varscan', 'lofreq'])
    output: join(config['variant_dir'], "variants.csv")
    run: aggregate_csv(input, output)



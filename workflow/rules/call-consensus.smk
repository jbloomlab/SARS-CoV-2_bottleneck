### ======= Generate consensus sequences from aligned viral reads. ======= ###
#
# Author: Will Hannon 
# Email: wwh22@uw.edu
# Date: 01/15/2021
#

def aggregate_fastqs(wildcards):
    """ 
    Get the fastq files from the same library for input to concatenation.
    """
    # Read the metadata into Pandas dataframe
    samples_df = pd.read_csv(config['samples']['file'])
    # Get the viral genome for a given accession
    library = np.unique(samples_df.loc[samples_df.Run == wildcards.accession, ["Library"]].values.flatten()).tolist()

    return expand(join(config['filter_dir'], "{{accession}}-{library}", "{{accession}}-{library}.filtered.fastq.gz"), library=library)


rule concatenate_fastqs: 
    """
    Join fastqs from the same library together 
    before de novo assembly of the consensus 
    sequence.
    """
    input: aggregate_fastqs
    output: join(config["consensus_dir"], "{accession}", "{accession}.merged.fastq.gz")
    shell: "cat {input} > {output}"


rule tadpole_error_correct: 
    """
    Using `tadpole.sh` from the BBmap suite,
    correct errors in the reads at the ends
    and in the middle of the reads. Use a kmer
    that at least less than 1/2 the length of the 
    read, ideally 1/3 read length. 
    """
    input: join(config["consensus_dir"], "{accession}", "{accession}.merged.fastq.gz")
    output: join(config["consensus_dir"], "{accession}", "{accession}.ec.fastq.gz")
    threads: config['threads']['max_cpu']
    params: kmer = config['kmer']['error_correct']
    conda: '../envs/filter.yml'
    shell: "tadpole.sh in={input} out={output} mode=correct k={params.kmer}"


rule tadpole_extend_reads:
    """
    Using `tadpole.sh` from the BBmap suite, 
    extend the error-corrected reads to longer
    reads that can be assembled into contigs. 
    """
    input: join(config["consensus_dir"], "{accession}", "{accession}.ec.fastq.gz")
    output: join(config["consensus_dir"], "{accession}", "{accession}.extended.fastq.gz")
    threads: config['threads']['max_cpu']
    params: kmer = config['kmer']['extend_reads']
    conda: '../envs/filter.yml'
    shell: "tadpole.sh in={input} out={output} mode=extend k={params.kmer} el=50 er=50"


rule tadpole_assemble_reads:
    """
    Using `tadpole.sh` from the BBmap suite, 
    assemble the error corrected and extended
    reads into longer contigs which can then 
    be polished into a full assembly.
    """
    input: join(config["consensus_dir"], "{accession}", "{accession}.extended.fastq.gz")
    output: join(config["consensus_dir"], "{accession}", "{accession}.contigs.fa")
    threads: config['threads']['max_cpu']
    params: kmer = config['kmer']['assemble_contigs']
    conda: '../envs/filter.yml'
    shell: "tadpole.sh in={input} out={output} k={params.kmer}"


rule bcftools_consensus:
    """
    Incorporate variants into the reference sequence from lofreq. 
    Lofreq is used here becuase it is fairly accurate and it 
    has indel re-alignemnt and indels. 
    """
    input: 
        vcf=join(config['variant_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.lofreq.vcf"),
        genome=get_genome
    output: join(config['consensus_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.consensus.fa")
    params: name=">{accession}.{aligner}"
    conda: '../envs/consensus.yml'
    shell:
        """
        bgzip -c {input.vcf} > {input.vcf}.gz
        tabix -p vcf {input.vcf}.gz
        cat {input.genome} | bcftools consensus {input.vcf}.gz > {output}
        sed -i "1s/.*/{params.name}/" {output}
        """


rule join_consensus: 
    """
    Join the consensus sequences from bcftools 
    for multiple sequence alignment with MAFFT.
    """
    input: expand(join(config['consensus_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.consensus.fa"), accession=list(set(pd.read_csv(config['samples']['file'])['Run'])), aligner=['BWA'])
    output: join(config['consensus_dir'], "joined_consensus.fa")
    shell: "cat {input} > {output}"


rule align_consensus:
    """
    Perform multiple sequence alignment with 
    MAFFT. 
    """
    input: join(config['consensus_dir'], "joined_consensus.fa")
    output: join(config['consensus_dir'], "aligned_consensus.fa")
    conda: '../envs/consensus.yml'
    shell: "mafft {input} > {output}"


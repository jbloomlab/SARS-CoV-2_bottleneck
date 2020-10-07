### ======= Calculate the coverage in bins over the genome ======= ###
#
# Author: Will Hannon 
# Email: wwh22@uw.edu
# Date: 10/31/2020
#

rule bedtools_coverage:
    """ Calculate the average coverage over bins in the genome.
    """
    input: bam=join(config['align_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.virus.sorted.marked.bam"),
           bai=join(config['align_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.virus.sorted.marked.bam.bai"),
           genome=get_genome
    output: join(config['coverage_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.bedgraph")
    params: binsize = config['bin_size']
    conda: '../envs/samtools.yml'
    shell:
        """
        awk 'FS=OFS="\t"{{print $1, 0, $2}}' {input.genome}.fai \
            | bedtools makewindows -b - -w {params.binsize} \
            | bedtools coverage -a stdin -b {input.bam} -mean \
            > {output}
            
        sed -i "s/$/\t{wildcards.accession}/" {output} 
        """


rule merge_coverage:
    """ Merge the coverage tables for all of the accessions into a single file.
    """
    input: expand(join(config['coverage_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.bedgraph"), accession=pd.read_csv(config['samples']['file'])['Run'], aligner=['BWA'])
    output: bedgraph=join(config['coverage_dir'], "merged.bedgraph"),
            header=temp(join(config['coverage_dir'], "merged.bedgraph.tmp"))
    shell:
        """
        cat {input} > {output.header}

        awk 'BEGIN{{print "Name\tStart\tStop\tDepth\tAccession"}}1' {output.header} > {output.bedgraph}
        """


rule samtools_depth:
    """ Calculate the depth over each position filtering by the phred base score. 
    """
    input: bam=join(config['align_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.virus.sorted.marked.bam"),
           bai=join(config['align_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.virus.sorted.marked.bam.bai")
    output: join(config['coverage_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.depth")
    params: score=config['BQ'],
            binsize=config['bin_size']
    conda: '../envs/samtools.yml'
    shell: 
        """
        samtools depth -q {params.score} -g DUP {input.bam} \
            | awk '{{sum+=$3}} NR%50==0 {{print $2-50 "\t" $2 "\t" sum/50 "\t"; sum=0}}' - \
            > {output}

        sed -i "s/$/\t{wildcards.accession}/" {output} 
        """


rule merge_depth:
    """ Merge the samtools depth tables for all of the accessions into a single file.
    """
    input: expand(join(config['coverage_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.depth"), accession=pd.read_csv(config['samples']['file'])['Run'], aligner=['BWA'])
    output: depth=join(config['coverage_dir'], "merged.depth"),
            header=temp(join(config['coverage_dir'], "merged.depth.tmp"))
    shell:
        """
        cat {input} > {output.header}

        awk 'BEGIN{{print "Start\tStop\tDepth\tAccession"}}1' {output.header} > {output.depth}
        """

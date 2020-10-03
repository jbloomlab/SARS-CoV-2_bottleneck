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


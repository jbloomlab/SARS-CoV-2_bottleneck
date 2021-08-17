### ======= Calculate the coverage in bins over the genome ======= ###
#
# Author: Will Hannon 
# Email: wwh22@uw.edu
# Date: 10/31/2020
#

rule samtools_depth:
    """ Calculate the depth over each position filtering by the phred base score. 
    """
    input: bam=join(config['align_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.sorted.merged.bam"),
           bai=join(config['align_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.sorted.merged.bam.bai")
    output: join(config['coverage_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.depth")
    params: score=config['BQ'],
            binsize=config['bin_size']
    conda: '../envs/samtools.yml'
    shell: 
        """
        samtools depth -m 0 -a -q {params.score} -g DUP {input.bam} \
            | awk '{{sum+=$3}} NR%50==0 {{print $2-50 "\t" $2 "\t" sum/50 "\t"; sum=0}}' - \
            > {output}

        sed -i "s/$/\t{wildcards.accession}/" {output} 
        """


rule merge_depth:
    """ Merge the samtools depth tables for all of the accessions into a single file.
    """
    input: expand(join(config['coverage_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.depth"), accession=samples, aligner=['BWA'])
    output: depth=join(config['coverage_dir'], "merged.depth"),
            header=temp(join(config['coverage_dir'], "merged.depth.tmp"))
    shell:
        """
        cat {input} > {output.header}

        awk 'BEGIN{{print "Start\tStop\tDepth\tAccession"}}1' {output.header} > {output.depth}
        """

rule average_depth:
    input: bam=join(config['align_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.sorted.merged.bam"),
           bai=join(config['align_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.sorted.merged.bam.bai")
    output: join(config['coverage_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.average.depth")
    params: score=config['BQ']
    conda: '../envs/samtools.yml'
    shell:
        """
        samtools depth -a -q {params.score} -g DUP {input.bam} \
            | awk '{{ if ($3 >= 100) sum+=1}} END {{print sum/NR*100}}' - \
            > {output}

        sed -i "s/$/\t{wildcards.accession}/" {output} 
        """

rule merge_average_depth:
    """ Merge the samtools depth tables for all of the accessions into a single file.
    """
    input: expand(join(config['coverage_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.average.depth"), accession=samples, aligner=['BWA'])
    output: depth=join(config['coverage_dir'], "merged.average.depth"),
            header=temp(join(config['coverage_dir'], "merged.average.depth.tmp"))
    shell:
        """
        cat {input} > {output.header}

        awk 'BEGIN{{print "Percent\tAccession"}}1' {output.header} > {output.depth}
        """

rule coverage_stats: 
    input: expand(join(config['align_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.sorted.merged.bam"), accession=samples, aligner=['BWA'])
    output: join(config['coverage_dir'], "coverage.stats")
    params: score=config['BQ']
    conda: '../envs/samtools.yml'
    shell:
        """
        echo "rname\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq\tfilename" > {output}
        
        list="{input}"

        for f in $list; do
            echo -e "$(samtools coverage -H -Q {params.score} --excl-flags UNMAP,SECONDARY,QCFAIL $f)\t$f" >> {output}
        done
        """


### ======= Pileup Test ======= ###

rule generate_pileup:
    """
    """
    input: 
        bam = join(config['split_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.virus.bam"),
        bai = join(config['split_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.virus.bam.bai"),
        genome = join(config['index_dir']['samtools'], "SARS2.fa")
    output: join(config['pileup_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.mpileup.txt")
    conda: '../envs/samtools.yml'
    shell: "samtools mpileup -d 0 -E -q 30 -Q 30 -f {input.genome} {input.bam} -O -s --reverse-del -a -o {output}"

rule parse_pileup:
    """
    """
    input: join(config['pileup_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.mpileup.txt")
    output: join(config['pileup_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.mpileup.csv")
    script: '../scripts/process_pileup.py'

rule parse_pileup_test:
    """
    """
    input: join(config['pileup_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.mpileup.txt")
    output: join(config['pileup_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.will.mpileup.csv")
    shell: 
        """
        python workflow/scripts/ParsePileup.py --input {input} --output {output} -si
        """

rule aggregate_pileup:
    """
    """
    input: expand(join(config['pileup_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.mpileup.csv"), accession=pd.read_csv(config['samples']['file'])['Run'], aligner=['BWA'])
    output: join(config['pileup_dir'], "raw-variants.csv")
    run:
        paths = []
        for f in input:
            try:
                pd.read_csv(f)
                paths.append(f)
            except:
                pass
        df = pd.concat(map(pd.read_csv, paths))
        df.to_csv(output[0], index = False)
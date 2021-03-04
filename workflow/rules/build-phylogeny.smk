### ======= Generate consensus sequences and make a phylogeny. ======= ###
#
# Author: Will Hannon 
# Email: wwh22@uw.edu
# Date: 01/29/2021
#

rule make_boat_consensus:
    input: bam=expand(join(config['align_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.sorted.marked.merged.bam"), accession=list(set(pd.read_csv(config['samples']['file'])['Run'])), aligner=['BWA']),
           bai=expand(join(config['align_dir'], "{aligner}", "{accession}", "{accession}.{aligner}.sorted.marked.merged.bam.bai"), accession=list(set(pd.read_csv(config['samples']['file'])['Run'])), aligner=['BWA']), 
           genome=join(config['index_dir']['samtools'], 'SARS2.fa')
    output: join(config['consensus_dir'], "merged.consenus.fa")
    conda: "../envs/pysam.yml"
    script: "../scripts/make-consensus-sequence.py"

rule iqtree_boat_phylogeny:
    input: join(config['consensus_dir'], "merged.consenus.fa")
    output: join(config['consensus_dir'], "merged.consenus.fa.iqtree")
    conda: "../envs/iqtree.yml"
    shell: "iqtree -s {input}"


rule align_GISAID_sequences:
    input: seqs=join(config['gisaid_dir'], "GISAID.fasta"),
           reference=get_genome
    output: join(config['gisaid_dir'], "GISAID.aligned.fasta")
    conda: "../envs/iqtree.yml"
    threads: 4
    shell: "mafft --6merpair --thread {threads} --keeplength --addfragments {input.seqs} {input.reference} > {output}"


rule iqtree_GISAID_phylogeny: 
    input: join(config['gisaid_dir'], "GISAID.aligned.fasta")
    output: join(config['gisaid_dir'], "GISAID.aliged.fasta.iqtree")
    conda: "../envs/iqtree.yml"
    shell: "iqtree -s {input} -m GTR+I+G"


def hamming_distance(chaine1, chaine2):
    """
    Calculate the edit distance between masked consensus genomes.
    """
    return sum(c1 != c2 for c1, c2 in zip(chaine1, chaine2))


rule calculate_edit_distance():
    input: join(config['consensus_dir'], "merged.consenus.fa")
    output: join(config['consensus_dir'], "edit_distance.csv")
    conda: "../envs/pysam.yml"
    run:
        # Import the fasta sequences
        fasta_sequences = SeqIO.parse(open(input),'fasta')
        # Save the names and the genomes in a list
        names = []
        genomes = []
        # Iterate over all of the fasta sequences
        for fasta in fasta_sequences:
            names.append(fasta.id)
            genomes.append(str(fasta.seq))
        # Save the hamming distances in a list
        distances = []
        for m, i in zip(names, genomes):
            for n, j in zip(names, genomes):
                distances.append((f"{m}", f"{n}",hamming_distance(i,j)))    

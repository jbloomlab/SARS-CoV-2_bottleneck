### ======= Generate consensus sequences and make a phylogeny. ======= ###
#
# Author: Will Hannon 
# Email: wwh22@uw.edu
# Date: 01/29/2021
#

rule get_blast_sequences:
    input: metadata = join(config['gisaid_dir'], "metadata_2021-03-04_10-31.tsv"),
           sequences = join(config['gisaid_dir'], "sequences_2021-03-04_08-34.fasta")
    output: blast_seqs = join(config['gisaid_dir'], "bastn_gisaid_ids.csv"), 
            boat_seqs = join(config['gisaid_dir'], "boat_gisaid_ids.csv")
    conda: "../envs/r.yml"
    script: "../scripts/blastn_gisaid_ids.R"


rule make_blast_database: 
    input: join(config['gisaid_dir'], "blastn_covid_database.fasta")
    output: join(config['gisaid_dir'], "blastn_covid_database.fasta.nhr") 
    conda: "../envs/iqtree.yml"
    shell: "makeblastdb -in {input} -dbtype nucl"


rule search_blast_database: 
    input: query = join(config['gisaid_dir'], "boat_genomes.fasta"),
            database = join(config['gisaid_dir'], "blastn_covid_database.fasta.nhr") 
    output: join(config['gisaid_dir'], "similar_to_boat_genomes.out")
    params: blastdb = join(config['gisaid_dir'], "blastn_covid_database.fasta")
    conda: "../envs/iqtree.yml"
    shell: "blastn -db {params.blastdb} -query {input.query} -out {output}"


rule parse_blast_results:
    input: blast_results = join(config['gisaid_dir'], "similar_to_boat_genomes.out"), 
           metadata = join(config['gisaid_dir'], "metadata_2021-03-04_10-31.tsv")
    output: join(config['gisaid_dir'], "smilar_to_boat_genomes.csv")
    params: max_included = 50
    run: 
        with open(str(input.blast_results), 'r') as infile:     
            data = False
            dataframelist = []
            counter = 0
            for line in infile:
                if counter == int(params.max_included):
                    break
                current_line = line.strip()
                if current_line.startswith("Sequences producing significant alignments:"):
                    data = True
                if data:
                    row = current_line.split()
                    if len(row) == 3:
                        counter += 1
                        dataframelist.append(row)
        df = pd.DataFrame(dataframelist, columns = ['strain', 'score', 'value'])
        metadata_df = pd.read_csv(input.metadata, sep='\t', low_memory=False) 
        metadata_df = metadata_df[metadata_df.strain.isin(df.strain)]
        metadata_df.to_csv(str(output), index=False )


rule align_GISAID_sequences:
    input: seqs=join(config['gisaid_dir'], "phylogeny.fasta"),
           reference=join(config['ref_dir'], "SARS2.fa")
    output: join(config['gisaid_dir'], "phylogeny.aligned.fasta")
    conda: "../envs/iqtree.yml"
    threads: 4
    shell: "mafft --6merpair --thread {threads} --keeplength --addfragments {input.seqs} {input.reference} > {output}"


rule iqtree_GISAID_phylogeny: 
    input: join(config['gisaid_dir'], "phylogeny.aligned.fasta")
    output: join(config['gisaid_dir'], "phylogeny.aligned.fasta.treefile")
    conda: "../envs/iqtree.yml"
    shell: "iqtree -s {input} -m GTR+I+G"



def hamming_distance(chaine1, chaine2):
    """
    Calculate the edit distance between masked consensus genomes.
    """
    return sum(c1 != c2 for c1, c2 in zip(chaine1, chaine2))


rule calculate_edit_distance:
    input: join(config['consensus_dir'], "merged.consenus.fa")
    output: join(config['consensus_dir'], "edit_distance.csv")
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

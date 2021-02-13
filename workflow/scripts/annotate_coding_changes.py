"""
Annotate coding changes in the SARS-CoV-2 genome. 
"""
__author__ = "Will Hannon"
__copyright__ = "Copyright 2020 Will Hannon"
__email__ = "wwh22@uw.edu"
__license__ = "MIT"

from Bio import SeqIO,  Entrez

def parse_genbank(email = "wwh22@uw.edu", ref_id = "NC_045512.2"): 
    """
    This function parses the genbank file assocated with the reference ID that gets passed to 
    it and builds a dictionary of coding sequence indicies that can easily be used to extract 
    the coding sequences from the reference genome. The reference ID should match that used
    to create the BAM file. 
    
    Parameters
    ----------
    email : str
        email of the user, necessary for accessing NCBI servers
        
    ref_id: str
        name of the reference genome

    Returns
    -------
    dict
        dictionary with key = gene names, value = Bio.SeqFeature.SeqFeature object
        
    Bio.Seq.Seq
        sequence of the reference in the genbank record
    
    """
    ## ============ Fetch genbank record ============ ##
    # Set email 
    Entrez.email = email
    # Make handel object 
    handle = Entrez.efetch(db="nuccore", id=ref_id, rettype="gb", retmode="text")
    # Save the record -- only extract first record (there should be only one)
    record = next(SeqIO.parse(handle, "gb"))
    
    ## ============ Parse genbank record ============ ##
    # Dictionary to hold the open reading frames
    ORFS = dict()
    for feature in record.features:
        # Only extract the coding sequences
        if feature.type == "CDS":  
            # Special considerations for overlapping ORF
            if feature.qualifiers.get("gene")[0] == "ORF1ab":
                # Get the open reading frame that contains the ribosomal slippage
                if "-1 ribosomal frameshift" in str(feature.qualifiers.get("note")): 
                    # Extract the non-overlapping and frameshifted indices
                    name = "ORF1ab"
                    ORFS[name] = feature
                # Get the open reading frame that just contains the 'a' portion
                else:
                    # Extract the non-overlapping and frameshifted indices
                    name = "ORF1a"
                    ORFS[name] = feature
            # Iterate ove the remaining trivial CDS   
            else:
                # Build the lookup dictionary with the normal sequences
                name = feature.qualifiers.get("gene")[0]
                ORFS[name] = feature
    # Return Lookup dictionary
    return ORFS, record.seq


def translate(codon):
    """
    Translate a three letter DNA string into 
    a one letter amino acid code. 

    Parameters
    ----------
    codon : str
        three letter DNA sequence

    Returns
    -------
    str
        one letter amino acid code

    Raises
    ------
    AssertionError
        error if codon sequence is invalid
        
    """
    
    table = { 
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                  
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
        'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*', 
        'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W', 
    } 
    
    assert codon in table.keys(), "Not a valid codon sequence."
    
    return table[codon]


def annotate_effect(cds_dict, genome, snp): 
    """
    Translate a three letter DNA string into 
    a one letter amino acid code. 

    Parameters
    ----------
    cds_dict : dict
        dictionary with gene name and Bio.SeqFeature.SeqFeature object
        
    genome : Bio.Seq.Seq
        object holding the sequence of the reference genome
        
    snp : tuple
        contains the position of the SNP (1-indexed) and alternative allele

    Returns
    -------
    tuple
        codon position (int), coding change (str), effect (str)
    """
    # List to save the coding effect
    coding_effect = []
    
    # Change the SNP position from 1-indexed to 0-indexed
    snp = (snp[0]-1, snp[1])
    
    # Determine which genes the SNP is located in
    genes = []
    for k,v in cds_dict.items():
        if snp[0] in range(v.location.start, v.location.end): 
            genes.append(k)
    # Check that SNP is in a gene
    if genes: 
        # Some SNPs will be in more than one gene, SARS has overlaping ORFs
        for gene in genes: 
            gene_tuple = list(zip(list(cds_dict[gene].location), cds_dict[gene].location.extract(genome)))
            # Get the indicies relative to the gene, add 1 to get 1-indexed values
            indicies = [x + 1 for x, y in enumerate(gene_tuple) if y[0] == snp[0]]
            # Determine codon position from gene index
            for i in indicies:
                # First position in codon
                if i % 3 == 1:
                    codonpos = 1
                    wtcodon = [gene_tuple[i-1], gene_tuple[i], gene_tuple[i+1]]
                # Second position in codon
                elif i % 3 == 2:
                    codonpos = 2
                    wtcodon = [gene_tuple[i-2], gene_tuple[i-1], gene_tuple[i]]
                # Third position in codon 
                elif i % 3 == 0:
                    codonpos = 3
                    wtcodon = [gene_tuple[i-3], gene_tuple[i-2], gene_tuple[i-1]]
                
                # From the wt codon sequence, determine the alterative codon, coding change, and effect
                altcodon = [snp if i == (codonpos-1) else b for i, b in enumerate(wtcodon)]
                wtaa = translate("".join(y for x,y in wtcodon))
                altaa = translate("".join(y for x,y in altcodon))
                if wtaa == altaa:
                    effect = "synonymous"
                elif wtaa != altaa and altaa == '*':
                    effect = "nonsense"
                elif wtaa != altaa and altaa != '*':
                    effect = "missense"
                # Save the codon effects and information
                coding_effect.append((codonpos, f"{wtaa}{-(i // -3)}{altaa}", effect, gene))
    # If the SNP isn't in a gene, it's intergeneic and has no coding effect
    else:
        coding_effect.append(("NA", "NA", "NA", "intergeneic"))
    
    
    # Deal with SNPs in multiple genes with multiple effects 
    if len(coding_effect) == 1:
        return list(coding_effect[0])
    else: 
        if len(set([(a,b,c) for a,b,c,d in coding_effect])) == 1: 
            return list(list(set(coding_effect))[0])
        # TODO: Deal with ambiguous sequences
        else:
            return ["NA", "NA", "NA", "ambiguous"] 
            

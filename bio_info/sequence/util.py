amino_acids = ["G", "A", "S", "P", "V", "T", "C", "I",
               "L", "N", "D", "K", "Q", "E", "M", "H",
               "F", "R", "Y", "W"]


rna_amino_dict = {
    "AAA": "K", "AAC": "N", "AAG": "K", "AAU": "N",
    "ACA": "T", "ACC": "T", "ACG": "T", "ACU": "T",
    "AGA": "R", "AGC": "S", "AGG": "R", "AGU": "S",
    "AUA": "I", "AUC": "I", "AUG": "M", "AUU": "I",
    "CAA": "Q", "CAC": "H", "CAG": "Q", "CAU": "H",
    "CCA": "P", "CCC": "P", "CCG": "P", "CCU": "P",
    "CGA": "R", "CGC": "R", "CGG": "R", "CGU": "R",
    "CUA": "L", "CUC": "L", "CUG": "L", "CUU": "L",
    "GAA": "E", "GAC": "D", "GAG": "E", "GAU": "D",
    "GCA": "A", "GCC": "A", "GCG": "A", "GCU": "A",
    "GGA": "G", "GGC": "G", "GGG": "G", "GGU": "G",
    "GUA": "V", "GUC": "V", "GUG": "V", "GUU": "V",
    "UAA": "*", "UAC": "Y", "UAG": "*", "UAU": "Y",
    "UCA": "S", "UCC": "S", "UCG": "S", "UCU": "S",
    "UGA": "*", "UGC": "C", "UGG": "W", "UGU": "C",
    "UUA": "L", "UUC": "F", "UUG": "L", "UUU": "F"
}


amino_weight_dict = {
    "G": 57, "A": 71, "S": 87, "P": 97,
    "V": 99, "T": 101, "C": 103, "I": 113,
    "L": 113, "N": 114, "D": 115, "K": 128,
    "Q": 128, "E": 129, "M": 131, "H": 137,
    "F": 147, "R": 146, "Y": 163, "W": 186
}


def kmer_composition(k, dna, sort=False):
    kmers = list()
    for i in range(len(dna) - k+1):
        kmers.append(dna[i:i+k])
    if sort:
        kmers.sort()
    return kmers


def n_50(contigs):
    # TODO write function
    total_length = sum([len(c) for c in contigs])
    halfway = int(total_length / 2)
    length_passed = 0
    for contig in sorted(contigs, reverse=True):
        length_passed += len(contig)
        if length_passed > halfway:
            return len(contig)


def rna_to_amino(rna):
    amino = list()
    i = 0
    while i < len(rna):
        amino.append(rna_amino_dict[rna[i:i + 3]])
        i += 3
    return ''.join(amino)


def get_peptide_mass(peptide):
    mass = 0
    for amino in peptide:
        mass += amino_weight_dict[amino]
    return mass

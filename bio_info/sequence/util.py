def kmer_composition(k, dna):
    kmers = list()
    for i in range(len(dna) - k+1):
        kmers.append(dna[i:i+k])
    kmers.sort()
    return kmers
def kmer_composition(k, dna, sort=False):
    kmers = list()
    for i in range(len(dna) - k+1):
        kmers.append(dna[i:i+k])
    if sort:
        kmers.sort()
    return kmers
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

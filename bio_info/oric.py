def get_kmer_count_dict(dna, k, precision='exact', d=None):
    """
    Creates a dictionary of counts for each k-mer in a string.
    
    :param dna: a string of DNA
    :param k: the length of each substring
    :param precision: how much precision is required in obtaining kmer_count
            exact: exact string matches
            mismatch: string matches with d or fewer differences
            reverse: string matches reverse complement as well
            loose: includes mismatch and reverse options
    :param d: optional param indicating maximum number of mismatches (Hamming distance)
    
    :return: a dictionary of each k-mer and the number of times it appears, the largest count
    """
    k_mers = dict()
    largest_count = 0

    # Run for all precision levels
    for i in range(0, len(dna) - k):
        seq = dna[i:i + k]
        if seq in k_mers:
            k_mers[seq] += 1
            if k_mers[seq] > largest_count:
                largest_count = k_mers[seq]
        else:
            k_mers[seq] = 1

    if precision == 'mismatch' or precision == 'loose':
        for k_mer in k_mers.keys():
            indexes = find_approximate_match_indexes(dna, k_mer, d)
            k_mers[k_mer] = len(indexes)
            if k_mers[k_mer] > largest_count:
                largest_count = k_mers[k_mer]

    if precision == 'reverse' or precision == 'loose':
        new_kmers = dict()
        for k_mer in k_mers.keys():
            rc_kmer = reverse_complement(k_mer)
            # Only add if reverse complement isn't already in new dictionary
            if rc_kmer not in new_kmers:
                # Determine if reverse complement is also present in the DNA
                if rc_kmer in k_mers:
                    # Add reverse complement count to current count
                    total_count = k_mers[k_mer] + k_mers[rc_kmer]
                else:
                    total_count = k_mers[k_mer]
                # Add to new list of k-mers
                new_kmers[k_mer] = total_count
                if new_kmers[k_mer] > largest_count:
                    largest_count = new_kmers[k_mer]
        k_mers = new_kmers

    # Return dictionary of existing k-mers and their counts
    return k_mers, largest_count
    
    
def reverse_complement(dna):
    """
    Find the reverse complement of a given DNA string.
    
    :param dna: the string of DNA
    :return: the reverse complement of :param dna
    """
    tmp = "U"

    # Swap A & T
    c_dna = dna.replace("A", tmp)
    c_dna = c_dna.replace("T", "A")
    c_dna = c_dna.replace(tmp, "T")

    # Swap C & G
    c_dna = c_dna.replace("C", tmp)
    c_dna = c_dna.replace("G", "C")
    c_dna = c_dna.replace(tmp, "G")

    # Reverse string
    rc_dna = ""
    for i in range(len(c_dna)):
        rc_dna += c_dna[len(c_dna) - i - 1]

    return rc_dna
    
    
def find_pattern_clumps(dna, k, L, t):
    """
    Find all patterns forming (L, t)-clumps in a given sequence of DNA.
    
    :param dna: a string of DNA
    :param k: the length of each k-mer
    :param L: the size of the window determining clumps
    :param t: the number of required occurrences
    
    :return: a list of patterns that meet the clumping requirements
    """
    # Create dictionary of k-mers
    k_mers = dict()
    for i in range(len(dna) - k):
        seq = dna[i:i + k]
        if seq in k_mers:
            k_mers[seq] += 1
        else:
            k_mers[seq] = 1

    # Get list of more common k_mers
    sequences = set()
    for sequence, frequency in k_mers.items():
        if frequency >= t:
            sequences.add(sequence)

    #  Create list of indexes for each common k_mer
    k_mer_indexes = dict()
    for pattern in sequences:
        pattern_lookup = '(?=' + pattern + ')'  # to find overlapping matches
        indexes = [s.start() for s in re.finditer(pattern_lookup, dna)]
        k_mer_indexes[pattern] = indexes

    # Look for clumps
    sequences.clear()
    for pattern in k_mer_indexes:
        indexes = k_mer_indexes.get(pattern)
        for i in range(0, len(indexes) - t + 1):
            if indexes[i + t - 1] + k - indexes[i] <= L:  # make sure it completely fits
                sequences.add(pattern)    
                
                
def find_minimum_skew(dna, save_skew=False, file_name=None):
    """
    Generates a list of indexes with the minimum G - C skew.
    
    :param dna: a string of DNA
    :param save_skew: whether to save the skew list to a file
    :param file_name: the file to save the skew to
    
    :return: a list of DNA indexes with the lowest skew
    """
    skew = list()
    current_skew = 0
    lowest_skew = 0
    skew.append(current_skew)

    # Find the skew for each position
    for i in range(len(dna)):
        if dna[i] == 'C':
            current_skew -= 1
        elif dna[i] == 'G':
            current_skew += 1
        skew.append(current_skew)
        if lowest_skew > current_skew:
            lowest_skew = current_skew

    if save_skew:
        # Save skew to file
        with open(file_name, 'w') as s:
            for value in skew:
                s.write("{},".format(value))

    # Find the indexes with the lowest skew
    minimum_indexes = list()
    for i in range(len(skew)):
        if skew[i] == lowest_skew:
            minimum_indexes.append(i)

    return minimum_indexes



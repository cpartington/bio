import math
import random


def get_nucleotide_count(dna):
    """
    Count the number of each nucleotide in the given DNA string.
    
    :param dna: a string of DNA
    :return: the number of A's, C's, G's, and T's
    """
    a, c, g, t = 0, 0, 0, 0
    for i in range(len(dna)):
        val = dna[i].upper()
        if val == 'A':
            a += 1
        elif val == 'C':
            c += 1
        elif val == 'G':
            g += 1
        elif val == 'T':
            t += 1
    return a, c, g, t
    
    
def hamming_distance(dna1, dna2):
    """
    Finds the Hamming Distance between two strings of DNA.
    
    :param dna1: a string of DNA
    :param dna2: a string of DNA

    :return: the computed Hamming Distance
    """
    distance = 0
    for i in range(len(dna1)):
        if dna1[i] != dna2[i]:
            distance += 1

    return distance    
    
    
def pattern_to_number(pattern):
    """
    Converts a DNA string to a number.
    
    :param pattern:
    :return:
    """
    if len(pattern) == 0:
        return 0
    last = pattern[len(pattern) - 1]
    prefix = pattern[:len(pattern) - 1]
    return 4 * pattern_to_number(prefix) + _symbol_to_number_(last)


def _symbol_to_number_(symbol):
    if symbol == 'A':
        return 0
    elif symbol == 'C':
        return 1
    elif symbol == 'G':
        return 2
    elif symbol == 'T':
        return 3


def number_to_pattern(num, k):
    """
    Converts a number back to a DNA string.
    
    :param num:
    :param k:
    :return:
    """
    if k == 1:
        return _number_to_symbol_(num)
    prefix_num = int(num / 4)
    remainder = num % 4
    symbol = _number_to_symbol_(remainder)
    prefix_pattern = number_to_pattern(prefix_num, k - 1)
    return prefix_pattern + symbol


def _number_to_symbol_(num):
    if num == 0:
        return 'A'
    elif num == 1:
        return 'C'
    elif num == 2:
        return 'G'
    elif num == 3:
        return 'T'


def read_fasta(file_name, join=True):
    """
    Reads a fasta file and returns a single DNA sequence.

    :param file_name: the name of the fasta file
    :param join: whether to join the lines of DNA

    :return: a single string containing all DNA in the file
    """
    dna = list()
    with open(file_name) as f:
        for line in f:
            if not line.startswith('>'):
                dna.append(line.strip())
    if join:
        return ''.join(dna)
    else:
        return dna


def get_entropy(a_freq, c_freq, g_freq, t_freq):
    """ 
    Get entropy for a set of nucleotide frequencies. 
    """
    if a_freq == 0:
        a_ent = 0
    else:
        a_ent = a_freq * math.log(a_freq, 2)

    if c_freq == 0:
        c_ent = 0
    else:
        c_ent = c_freq * math.log(c_freq, 2)

    if g_freq == 0:
        g_ent = 0
    else:
        g_ent = g_freq * math.log(g_freq, 2)

    if t_freq == 0:
        t_ent = 0
    else:
        t_ent = t_freq * math.log(t_freq, 2)

    return -1 * (a_ent + c_ent + g_ent + t_ent)    


def get_d_neighborhood(dna, d):
    """

    
    :param dna:
    :param d:
    :return:
    """
    if d == 0:
        return dna
    if len(dna) == 1:
        return {'A', 'C', 'G', 'T'}
    neighborhood = set()
    suffix_neighbors = get_d_neighborhood(dna[1:], d)
    for pattern in suffix_neighbors:
        if hamming_distance(dna[1:], pattern) < d:
            for nucleotide in {'A', 'C', 'G', 'T'}:
                neighborhood.add('{}{}'.format(nucleotide, pattern))
        else:
            neighborhood.add('{}{}'.format(dna[0], pattern))
    return neighborhood


def find_approximate_match_indexes(dna, pattern, d):
    """
    Finds all indexes of substrings that match a pattern with at most d mismatches.
    
    :param dna: a string of DNA
    :param pattern: the pattern to be matched
    :param d: maximum number of mismatches (Hamming Distance)
    
    :return: a list of start indexes for matched strings
    """
    indexes = list()
    k = len(pattern)
    for i in range(len(dna) - k+1):
        substr = dna[i:i + k]
        if hamming_distance(pattern, substr) <= d:
            indexes.append(i)
    return indexes


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


# TODO add random DNA / peptide generation class
def generate(length, count=1, type="dna"):
    if type == "peptide":
        data = "ILVFMCAGPTSYWQNHEDKR"
    else:
        data = "ACGT"
    if count == 1:
        return ''.join(random.choice(data) for i in range(length))
    else:
        dna = list()
        for i in range(count):
            dna.append(''.join(random.choice(data) for i in range(length)))
        return dna

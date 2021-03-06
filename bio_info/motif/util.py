from .profile import Profile
from ..util import *


def entropy_score(motif_list):
    """
    Get entropy for a list of motifs.

    :param motif_list: a list of DNA strings of equal length
    :return: the sum of the entropy scores for each index / column of the DNA strings
    """
    ent_score = 0
    profile = Profile(motif_list)
    for i in range(len(profile.get(0))):
        ent_score += get_entropy(profile.get(0, i), profile.get(1, i),
                                 profile.get(2, i), profile.get(3, i))
    return ent_score


def mismatch_score(motif_list):
    """
    Get mismatch count for a list of motifs.

    :param motif_list: a list of DNA strings of equal length
    :return: the sum of the non-matching nucleotides in each position
    """
    total_mismatch = 0
    for i in range(len(motif_list[0])):
        col = [motif[i] for motif in motif_list]
        a, c, g, t = get_nucleotide_count(col)
        largest = max(a, c, g, t)
        total_mismatch += len(col) - largest
    return total_mismatch


def distance_pattern_strings(pattern, dna_list):
    """
    Finds the sum of the minimum Hamming distances between a DNA
    pattern and a list of longer DNA sequences

    :param pattern: a DNA sequence
    :param dna_list: a list of DNA strings of equal length

    :return: the sum of the minimum Hamming distance for each string
    """
    k = len(pattern)
    total_d = 0
    for dna in dna_list:
        hamming_d = float("inf")  # infinity
        for i in range(len(dna) - k+1):
            dist = hamming_distance(pattern, dna[i:i+k])
            if dist < hamming_d:
                hamming_d = dist
        total_d += hamming_d
    return total_d


def form_consensus_string(motif_list):
    """
    Generates a consensus string for each position in a motif matrix using the
    nucleotide frequencies.

    :param motif_list: a list of DNA strings of equal length
    :return: a single DNA sequence representing the most common nucleotide for
             each position in the given list of motifs
    """
    consensus = list()
    for i in range(len(motif_list[0])):
        a, c, g, t = get_nucleotide_count(''.join(motif[i] for motif in motif_list))
        high_count = max(a, c, g, t)
        if a == high_count:
            consensus.append("A")
        elif c == high_count:
            consensus.append("C")
        elif g == high_count:
            consensus.append("G")
        elif t == high_count:
            consensus.append("T")
    return ''.join(consensus)

import random
import numpy as np

from progressbar import progressbar

from .general import *


def find_profile_probable_kmer(profile, k, dna):
    """
    Finds the most likely k-mer given a profile.

    :param profile: a profile matrix stored as a list of lists
    :param k: the length of a k-mer
    :param dna: the DNA sequence to search
    :return: the k-mer in :param dna: with the highest probability
    """
    index = 0
    top_prob = 0
    for i in range(len(dna) - k+1):
        pattern = dna[i:i+k]
        prob = profile[pattern_to_number(pattern[0])][0]
        for j in range(1, len(pattern)):
            prob *= profile[pattern_to_number(pattern[j])][j]
        if prob > top_prob:
            top_prob = prob
            index = i
    return dna[index:index+k]


def entropy_score(motif_list):
    """ 
    Get entropy for a list of motifs. 
    
    :param motif_list: a list of DNA strings of equal length
    :return: the sum of the entropy scores for each index / column of the DNA strings
    """
    ent_score = 0
    profile = profile_matrix(motif_list)
    for i in range(len(profile[0])):
        ent_score += get_entropy(profile[0][i], profile[1][i],
                                 profile[2][i], profile[3][i])
    return ent_score


def mismatch_score(motif_list):
    """ 
    Get mismatch count for a list of motifs. 
    """
    total_mismatch = 0
    for i in range(len(motif_list[0])):
        col = [motif[i] for motif in motif_list]
        a, c, g, t = get_nucleotide_count(col)
        largest = max(a, c, g, t)
        total_mismatch += len(col) - largest
    return total_mismatch


def profile_matrix(motif_list):
    """
    Get profile matrix for a list of motifs.
    """
    k = len(motif_list[0])
    t = len(motif_list) + 4
    profile = list()
    profile.append(list())
    profile.append(list())
    profile.append(list())
    profile.append(list())

    for i in range(k):
        col = [motif[i] for motif in motif_list]
        a_num, c_num, g_num, t_num = get_nucleotide_count(col)
        profile[0].append((a_num + 1) / t)
        profile[1].append((c_num + 1) / t)
        profile[2].append((g_num + 1) / t)
        profile[3].append((t_num + 1) / t)

    return profile    


def greedy_motif_search(k, dna_list, scoring='entropy'):
    t = len(dna_list)
    # Generate initial motif matrix
    best_motifs = list()
    for dna in dna_list:
        pattern = dna[0:k]
        best_motifs.append(pattern)
    if scoring == 'mismatch':
        b_score = mismatch_score(best_motifs)
    else:
        b_score = entropy_score(best_motifs)

    # Iteratively find best motif matrix
    motif_list = ['A'] * t
    for i in range(len(dna_list[0]) - k+1):  # for each k-mer in first DNA string
        motif_list[0] = dna_list[0][i:i + k]
        for j in range(1, t):
            profile = profile_matrix(motif_list[:j])
            motif_list[j] = (find_profile_probable_kmer(
                              profile, k, dna_list[j]))
        # Get motif list score
        if scoring == 'mismatch':
            m_score = mismatch_score(motif_list)
        else:
            m_score = entropy_score(motif_list)
        if m_score < b_score:
            # Update best motifs
            best_motifs = motif_list.copy()
            b_score = m_score

    return best_motifs    


def randomized_motif_search(dna_list, k, n=1000, score='entropy'):
    """
    A wrapper for the randomized motif search algorithm.

    :param dna_list: a list of DNA strings
    :param k: length of k-mer
    :param n: number of iterations
    :param score: entropy or mismatch count

    :return: the set of motifs with the lowest score
    """
    best_motifs = list()
    best_score = float("inf")
    for i in progressbar(range(n)):
        motifs, low_score = _randomized_motif_search_(dna_list, k, score)
        if low_score < best_score:
            best_motifs = motifs
            best_score = low_score
    return best_motifs


def _randomized_motif_search_(dna_list, k, score):
    t = len(dna_list)
    # Get random k-mers from each string
    motifs = list()
    for i in range(t):
        index = random.randint(0, len(dna_list[0]) - k)
        motifs.append(dna_list[i][index:index+k])
    best_motifs = motifs.copy()
    if score == 'entropy':
        b_score = entropy_score(best_motifs)
    else:
        b_score = mismatch_score(best_motifs)
    while True:
        profile = profile_matrix(motifs)
        # Get profile-most-probable k-mer of each DNA sequence
        motifs.clear()
        for dna in dna_list:
            motifs.append(find_profile_probable_kmer(profile, k, dna))
        # Compare results
        if score == 'entropy':
            m_score = entropy_score(motifs)
        else:
            m_score = mismatch_score(motifs)
        if m_score < b_score:
            best_motifs = motifs.copy()
            b_score = m_score
        else:
            return best_motifs, b_score


def gibbs_sampler(dna_list, k, it=20, n=1000, score='entropy'):
    """
    Runs Gibbs sampler on a given list of DNA sequences to find the best
    set of motifs.

    :param dna_list: a list of DNA strings
    :param k: the desired k-mer length
    :param it: the number of sampler iterations
    :param n: to be deprecated
    :param score: entropy or mismatch count

    :return: the set of motifs with the lowest score
    """
    best_motifs = list()
    best_score = len(dna_list) * k  # impossibly large score
    for i in progressbar(range(it)):
        motifs, low_score = _gibbs_sampler_(dna_list, k, score, n)
        if low_score < best_score:
            best_motifs = motifs
            best_score = low_score
    return best_motifs


def _gibbs_sampler_(dna_list, k, score, n):
    t = len(dna_list)
    # Randomly select k-mer from each DNA string
    motifs = list()
    for i in range(t):
        index = random.randint(0, len(dna_list[0]) - k)
        motifs.append(dna_list[i][index:index+k])
    best_motifs = motifs.copy()
    # Get initial best motifs score
    if score == 'entropy':
        b_score = entropy_score(best_motifs)
    else:
        b_score = mismatch_score(best_motifs)

    # Replace motifs randomly
    for x in range(n):
        # Randomly remove one string
        i = random.randint(0, t-1)  # choose random string to remove
        # Generate profile
        profile = profile_matrix(motifs[:i] + motifs[i + 1:])
        # Generate probabilities for removed string
        removed = dna_list[i]
        probabilities = list()
        for m in range(len(removed) - k+1):
            kmer = removed[m:m+k]
            probabilities.append(kmer_profile_probability(profile, kmer))
        # Choose profile-random k-mer from removed DNA sequence
        j = _biased_random_(len(removed) - k, probabilities)
        motifs[i] = removed[j:j +k]
        # Compare scores
        if score == 'entropy':
            m_score = entropy_score(motifs)
        else:
            m_score = mismatch_score(motifs)
        if m_score < b_score:
            best_motifs = motifs.copy()
            b_score = m_score

    return best_motifs, b_score


def kmer_profile_probability(profile, kmer):
    """
    Find the probability of a k-mer given a profile.

    :param profile: a matrix of probabilities in the form of a list of lists
    :param kmer: the DNA string
    :return: the probability of :param kmer
    """
    prob = profile[pattern_to_number(kmer[0])][0]
    for i in range(1, len(kmer)):
        prob *= profile[pattern_to_number(kmer[i])][i]
    return prob


def _biased_random_(max_index, probabilities):
    total = sum(probabilities)
    probs = [p / total for p in probabilities]
    return np.random.choice(max_index+1, p=probs)


def distance_pattern_strings(pattern, dna_list):
    """
    Finds the sum of the minimum Hamming distances between a DNA
    pattern and a list of longer DNA sequences

    :param pattern: a DNA sequence
    :param dna_list: a list of DNA strings

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


def median_string_search(k, dna_list):
    """
    Finds the median string of a set of DNA strings.

    :param k: the desired k-mer length
    :param dna_list: a list of DNA strings

    :return: the median string
    """
    best_d = float("inf")  # infinity
    for i in range(1, 4 ** k):
        pattern = number_to_pattern(i, k)
        new_d = distance_pattern_strings(pattern, dna_list)
        if new_d < best_d:
            best_d = new_d
            median = pattern
    return median


def form_consensus_string(motif_list):
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

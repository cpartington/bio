from .util import *


def overlap_graph(pattern_list):
    adjacencies = dict()

    for i in range(len(pattern_list)):
        # Setup
        pattern = pattern_list[i]
        adjacencies[pattern] = list()
        suffix = pattern[1:]
        other_patterns = pattern_list[:i] + pattern_list[i+1:]
        # Check other patterns for matching prefixes
        for p in other_patterns:
            prefix = p[:-1]
            if suffix == prefix:
                # Add to adjacency list
                adjacencies.get(pattern).append(p)

    return adjacencies


def de_brujin_string(k, dna):
    adj = dict()
    edges = kmer_composition(k, dna)
    # Add first k-mer
    prefix = edges[0][:-1]
    suffix = edges[0][1:]
    adj[prefix] = list()
    adj.get(prefix).append(suffix)
    for i in range(1, len(edges)):
        prefix = suffix
        if prefix not in adj:
            adj[prefix] = list()
        suffix = edges[i][1:]
        adj.get(prefix).append(suffix)
    return adj

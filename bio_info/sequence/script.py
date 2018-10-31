from .graph import de_bruijn, eulerian_path


def reconstruct(dna_list):
    """
    Given a series of reads, reconstruct a genome using
    a de Bruijn graph and a Eulerian path algorithm.

    :param dna_list: a list of DNA string reads

    :return: the combined DNA sequence
    """
    print("Building graph...")
    g = de_bruijn(dna_list)

    print("Finding path...")
    path = eulerian_path(g)

    print("Building sequence")
    seq_list = list()
    for node in path:
        seq_list.append(node.label[0])
    seq_list.append(path[-1].label[1:])
    return "".join(seq_list)
from tabulate import tabulate

from ..general import *


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


def profile_to_string(profile):
    return tabulate(profile, tablefmt="plain")

from ..general import *


class Profile:

    def __init__(self, motif_list):
        self.profile = self.make(motif_list)

    def make(self, motif_list):
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

    def get(self, row=None, col=None):
        if row is None and col is None:
            return
        elif row is None:
            return [r[col] for r in self.profile]
        elif col is None:
            return self.profile[row]
        else:
            return self.profile[row][col]

    def to_string(self):
        string = list()
        for row in self.profile:
            for col in row:
                string.append("{:.3f}  ".format(col))
            string.append('\n')
        return ''.join(string).strip()

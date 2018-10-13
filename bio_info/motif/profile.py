from ..util import *


class Profile:

    def __init__(self, motif_list):
        val = motif_list[0][0]
        if type(val) == str:
            self.profile = self.from_dna_list(motif_list)
        elif type(val) == int:
            self.profile = self.from_count(motif_list)
        elif type(val) == float:
            self.profile = self.from_profile(motif_list)

    def from_dna_list(self, motif_list):
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

    def from_count(self, count_matrix):
        k = len(count_matrix[0])
        profile = list()
        profile.append(list())
        profile.append(list())
        profile.append(list())
        profile.append(list())

        for i in range(k):
            col = [count[i] for count in count_matrix]
            t = sum([col[0], col[1], col[2], col[3]]) + 4
            profile[0].append((col[0] + 1) / t)
            profile[1].append((col[1] + 1) / t)
            profile[2].append((col[2] + 1) / t)
            profile[3].append((col[3] + 1) / t)

        return profile


    def from_profile(self, profile):
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

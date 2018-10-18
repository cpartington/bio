from ..util import *


class Profile:
    # TODO option to build profile from peptides
    def __init__(self, data=None, k = None):
        if data is None and k is None:
            return
        if data is None:
            self.profile = self.from_count([[1] * k] * 4)
            return
        val = data[0][0]
        if type(val) == str:
            self.profile = self.from_dna_list(data)
        elif type(val) == int:
            self.profile = self.from_count(data)
        elif type(val) == float:
            self.profile = data
        self.m = 4
        self.n = len(data[0])

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

    def get(self, row=None, col=None):
        """
        Gets the data at the specified index(es).

        :param row: the row index
        :param col: the column index
        :return: a single value if row and column are specified
                 a list if only one is specified
        """
        if row is None and col is None:
            return
        elif row is None:
            return [r[col] for r in self.profile]
        elif col is None:
            return self.profile[row]
        else:
            return self.profile[row][col]

    def entropy_score(self):
        """
        Get entropy for a list of motifs.

        :return: the sum of the entropy scores for each column in the profile
        """
        ent_score = 0
        for i in range(len(self.get(0))):
            ent_score += get_entropy(self.get(0, i), self.get(1, i),
                                     self.get(2, i), self.get(3, i))
        return ent_score

    def kmer_probability(self, kmer):
        """
        Find the probability of a k-mer given a profile.

        :param profile: a matrix of probabilities in the form of a list of lists
        :param kmer: the DNA string
        :return: the probability of :param kmer
        """
        if len(kmer) != self.n:
            return
        prob = self.get(row=pattern_to_number(kmer[0]), col=0)
        for i in range(1, len(kmer)):
            prob *= self.get(row=pattern_to_number(kmer[i]), col=i)
        return prob

    def consensus_string(self):
        consensus = list()
        for i in range(self.n):
            col = self.get(col=i)
            a, c, g, t = col[0], col[1], col[2], col[3]
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

    def to_string(self):
        string = list()
        for row in self.profile:
            for col in row:
                string.append("{:.3f}  ".format(col))
            string.append('\n')
        return ''.join(string).strip()

    def max_entropy(self):
        """
        Get the maximum possible entropy based on the size of the profile.

        :return: maximum entropy
        """
        uniform_profile = Profile(self.from_count([[1] * self.n] * 4))
        return uniform_profile.entropy_score()

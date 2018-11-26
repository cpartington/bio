import math
from .util import amino_acids


amino_mass_dict = {
    "": 0, "G": 57, "A": 71, "S": 87, "P": 97,
    "V": 99, "T": 101, "C": 103, "I": 113,
    "L": 113, "N": 114, "D": 115, "K": 128,
    "Q": 128, "E": 129, "M": 131, "H": 137,
    "F": 147, "R": 156, "Y": 163, "W": 186
}


def get_peptide_mass(peptide):
    """
    Gets the mass of a given peptide.

    :param peptide: a peptide in string or list form; can
           be a combination of integer masses and peptide
           character representations

    :return: the total mass of the peptide
    """
    mass = 0
    if isinstance(peptide, int):
        return peptide
    for amino in peptide:
        if isinstance(amino, int):
            mass += amino
        else:
            mass += amino_mass_dict[amino]
    return mass


def expected_spectrum_size(peptide):
    return len(peptide) * (len(peptide) - 1) + 2


def expected_peptide_length(spectrum):
    a = 1
    b = -1
    c = -len(spectrum) + 2
    delta = (b ** 2) - (4 * a * c)
    return int(-b + math.sqrt(delta) / (2 * a))


class Spectrum:
    """
    A class to do mass spectrum operations.
    """
    def __init__(self, spectrum, from_peptide=False, cyclic=False):
        """
        Initializes spectrum class with already-existing spectrum or creates a
        new spectrum using a provided peptide.

        :param spectrum: the data, either a spectrum or a peptide
        :param from_peptide: if True, treat :param spectrum as a peptide
        :param cyclic: only used if :param from_peptide is True; determines
               how to build the spectrum from the provided peptide (cyclic or
               linear)
        """
        if from_peptide:
            self.spectrum = self.build_spectrum(spectrum, cyclic)
        else:
            self.spectrum = spectrum

    def build_spectrum(self, peptide, cyclic=False):
        """
        Builds a spectrum given a peptide.

        :param peptide: the peptide to build a spectrum for
        :param cyclic: if True, treat the string as circular

        :return: the created spectrum as a list of integers
        """
        masses = [0]
        if isinstance(peptide, int):
            masses += [peptide]
            return masses
        for i in range(len(peptide)):
            masses += [masses[i] + get_peptide_mass(peptide[i])]
        spec = [0]
        if cyclic:
            peptide_mass = masses[len(peptide)]
        for i in range(len(peptide) - 1):
            for j in range(i + 1, len(peptide)):
                spec += [masses[j] - masses[i]]
                if cyclic:
                    spec += [peptide_mass - (masses[j] - masses[i])]
        spec += [masses[len(peptide)]]
        spec.sort()
        return spec

    def score(self, peptide, cyclic=False):
        """
        Gives a score of a peptide against the spectrum, where the
        score is calculated as the number of common masses between the
        peptide and the spectrum.

        :param peptide: the peptide string to generate a score for
        :param cyclic: whether to consider the peptide cyclic when
               generating its spectrum

        :return: the determined score
        """
        p_spectrum = self.build_spectrum(peptide, cyclic)
        s_spectrum = self.spectrum[:]
        common = 0
        for p in p_spectrum:
            if p in s_spectrum:
                common += 1
                s_spectrum.remove(p)
        return common

    def expected_peptide_length(self):
        """
        Uses the quadratic formula to find an approximation of the
        expected peptide length given the size of the spectrum.

        :return: expected peptide length
        """
        a = 1
        b = -1
        c = -len(self.spectrum) + 2
        delta = (b ** 2) - (4 * a * c)
        return int(-b + math.sqrt(delta) / (2 * a))

    # TODO combine sequencing functions
    def cyclopeptide_sequencing(self, return_as='peptides'):
        peptides = list()
        peptides.extend([a for a in amino_acids if get_peptide_mass(a) in self.spectrum])
        p_size = expected_peptide_length(self.spectrum)
        for _ in range(1, p_size):
            peptides = self._c_extend_(peptides)
        # Compare spectrums
        final_peptides = list()
        for peptide in peptides:
            if self.build_spectrum(peptide, cyclic=True) == self.spectrum:
                final_peptides.append(peptide)
        if return_as == 'peptides':
            return final_peptides
        elif return_as == 'masses':
            masses = set(['-'.join([str(get_peptide_mass(a)) for a in p]) for p in final_peptides])
            return masses

    def _c_extend_(self, peptides):
        new_peptides = list()
        for peptide in peptides:
            for a in amino_acids:
                if get_peptide_mass(peptide + a) in self.spectrum:
                    new_peptides.append(peptide + a)
        return new_peptides

    def leaderboard_sequencing(self, n):
        total_mass = max(self.spectrum)
        leaderboard = [""]
        leader_peptide = ""
        leader_score = self.score(leader_peptide, cyclic=True)

        # Look for peptides until no more peptides of correct mass can be found
        while len(leaderboard) > 0:
            leaderboard = [p + a for a in amino_acids for p in leaderboard]
            new_board = []
            new_score = []
            for p in leaderboard:
                p_mass = get_peptide_mass(p)
                if p_mass == total_mass:
                    p_score = self.score(p, cyclic=True)
                    if p_score > leader_score:
                        leader_peptide = p
                        leader_score = p_score
                elif p_mass < total_mass:
                    new_board += [p]
                    new_score += [self.score(p, self.spectrum)]
            leaderboard = self._l_cut_(new_board, new_score, n)
        return leader_peptide

    def _l_cut_(self, leaderboard, scores, n):
        best = [peptide for _, peptide in sorted(zip(scores, leaderboard), reverse=True)]
        cutoff = n
        while len(best) > cutoff and self.score(best[cutoff - 1]) == self.score(best[cutoff]):
            cutoff += 1
        print(n, len(best[:cutoff]), len(best))
        return best[:cutoff]

    def convolution(self):
        diff_counts = {}
        # diff_matrix = []
        for i in range(len(self.spectrum)):
            # diff_matrix += [[]]
            for j in range(len(self.spectrum)):
                if i > j:
                    diff = self.spectrum[i] - self.spectrum[j]
                    # diff_matrix[i] += [diff]
                    if diff in diff_counts:
                        diff_counts[diff] += 1
                    elif diff != 0:
                        diff_counts[diff] = 1
        # for l in diff_matrix:
        #     print(*l)
        return diff_counts

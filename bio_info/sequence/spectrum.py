import math
from .util import amino_acids


amino_mass_dict = {
    "G": 57, "A": 71, "S": 87, "P": 97,
    "V": 99, "T": 101, "C": 103, "I": 113,
    "L": 113, "N": 114, "D": 115, "K": 128,
    "Q": 128, "E": 129, "M": 131, "H": 137,
    "F": 147, "R": 156, "Y": 163, "W": 186
}


def get_peptide_mass(peptide):
    mass = 0
    for amino in peptide:
        mass += amino_mass_dict[amino]
    return mass


def spectrum(peptide, cyclic=False):
    masses = [0]
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


def expected_spectrum_size(peptide):
    return len(peptide) * (len(peptide) - 1) + 2


def expected_peptide_length(spectrum):
    a = 1
    b = -1
    c = -len(spectrum) + 2
    delta = (b ** 2) - (4 * a * c)
    return int(-b + math.sqrt(delta) / (2 * a))


def cyclopeptide_sequencing(spec, return_as='peptides'):
    peptides = list()
    peptides.extend([a for a in amino_acids if get_peptide_mass(a) in spec])
    p_size = expected_peptide_length(spec)
    for _ in range(1, p_size):
        peptides = _extend_peptides_(peptides, spec)
    # Compare spectrums
    final_peptides = list()
    for peptide in peptides:
        if spectrum(peptide, cyclic=True) == spec:
            final_peptides.append(peptide)
    if return_as == 'peptides':
        return final_peptides
    elif return_as == 'masses':
        masses = set(['-'.join([str(get_peptide_mass(a)) for a in p]) for p in final_peptides])
        return masses


def _extend_peptides_(peptides, spectrum):
    new_peptides = list()
    for peptide in peptides:
        for a in amino_acids:
            if get_peptide_mass(peptide + a) in spectrum:
                new_peptides.append(peptide + a)
    return new_peptides


def score(peptide, m_spec, cyclic=False):
    p_spectrum = spectrum(peptide, cyclic=False)
    return len([s for s in p_spectrum if s in m_spec])

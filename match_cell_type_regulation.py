import numpy as np
from consts import PARTIAL_SUMS_GH37, PARTIAL_SUMS_GH38

PARTIAL_SUMS = {
    37: PARTIAL_SUMS_GH37,
    38: PARTIAL_SUMS_GH38
}


def annotate_cell_types_vector(nucleotides_db, annotations_db, reference_genome, chromosome, start_position, end_position):
    partial_sums = PARTIAL_SUMS[reference_genome]
    index = partial_sums[chromosome - 1]
    first_nuc = annotations_db[nucleotides_db[index+start_position-1]]
    last_nuc = annotations_db[nucleotides_db[index+end_position]]
    return [max(first_nuc[0], last_nuc[0]), max(first_nuc[1], last_nuc[1]), np.amax(np.array([first_nuc[2], last_nuc[2]]), axis=0)]


def annotate_cell_types_matrix(nucleotides_db, annotations_db, reference_genome, chromosome, start_position, end_position):
    partial_sums = PARTIAL_SUMS[reference_genome]
    index = partial_sums[chromosome - 1]
    result = annotations_db[nucleotides_db[index+start_position-1:index+end_position]]
    return result


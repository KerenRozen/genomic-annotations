import numpy as np
from consts import PARTIAL_SUMS_GH37, PARTIAL_SUMS_GH38

PARTIAL_SUMS = {
    37: PARTIAL_SUMS_GH37,
    38: PARTIAL_SUMS_GH38
}


def annotate_cell_types_vector(nucleotides_db, annotations_db, reference_genome, chromosome, start_position, end_position):
    """
    Receives nucleotides_db, annotations_db, reference_genome and chromosome, start_position and end_position of a sample.
    Calculates the output based on the first and last nucleotides of the sample.
    :return: a 3-long list:
            list[0] = sum_score in format numpy.float64
            list[1] = mean_score in format numpy.float64
            list[2] = 164 values from 0 to 9 corresponding to the label for each cell type in format numpy.ndarray with shape (164,) and dtype=numpy.uint8.
                The labels encoding: 0 - no label, 9 - Promoter, 8 -  Enhancer, 7 - Transcribed, 6 - Bivalent, 5 - RegPermissive,
                4 - ConstitutiveHet, 3 - FacultativeHet, 2 - LowConfidence, 1 -  Quiescent
    """
    partial_sums = PARTIAL_SUMS[reference_genome]
    index = partial_sums[chromosome - 1]
    first_nuc = annotations_db[nucleotides_db[index+start_position-1]]
    last_nuc = annotations_db[nucleotides_db[index+end_position]]
    return [sum([first_nuc[0], last_nuc[0]])/2, sum([first_nuc[1], last_nuc[1]])/2, np.amax(np.array([first_nuc[2], last_nuc[2]]), axis=0)]


def annotate_cell_types_matrix(nucleotides_db, annotations_db, reference_genome, chromosome, start_position, end_position):
    """
    Receives nucleotides_db, annotations_db, reference_genome and chromosome, start_position and end_position of a sample.
    :return: a matrix represented as a numpy ndarray with shape(sample_length, 3) and dtype=numpy.object.
            Each row is a single nucleotide annotation in format numpy.ndarray with shape(3,) and dtype=numpy.object:
                row[0] = sum_score in format numpy.float16
                row[1] = mean_score in format numpy.float16
                row[2] = 164 values from 0 to 9 corresponding to the label for each cell type in format numpy.ndarray with shape (164,) and dtype=numpy.uint8
                    The labels encoding: 0 - no label, 9 - Promoter, 8 -  Enhancer, 7 - Transcribed, 6 - Bivalent, 5 - RegPermissive,
                    4 - ConstitutiveHet, 3 - FacultativeHet, 2 - LowConfidence, 1 -  Quiescent
    """
    partial_sums = PARTIAL_SUMS[reference_genome]
    index = partial_sums[chromosome - 1]
    result = annotations_db[nucleotides_db[index+start_position-1:index+end_position]]
    return result


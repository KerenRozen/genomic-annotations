import numpy as np
from consts import PARTIAL_SUMS_GH37, PARTIAL_SUMS_GH38

PARTIAL_SUMS = {
    37: PARTIAL_SUMS_GH37,
    38: PARTIAL_SUMS_GH38
}


def annotate_methylation_vector(nucleotides_db, annotations_db, reference_genome, chromosome, start_position, end_position):
    """
    Receives nucleotides_db, annotations_db, reference_genome and chromosome, start_position and end_position of a sample.
    Calculates the output based on the first and last nucleotides of the sample (calculating the mean score for each cell type).
    :return: a numpy.ndarray with shape (39,0) and elements dtype=numpy.float16.
    """
    partial_sums = PARTIAL_SUMS[reference_genome]
    index = partial_sums[chromosome - 1]
    first_nuc = annotations_db[nucleotides_db[index + start_position - 1]]
    last_nuc = annotations_db[nucleotides_db[index + end_position - 1]]
    return np.mean(np.array([first_nuc, last_nuc]), axis=0)


def annotate_methylation_matrix(nucleotides_db, annotations_db, reference_genome, chromosome, start_position, end_position):
    """
    Receives nucleotides_db, annotations_db, reference_genome and chromosome, start_position and end_position of a sample.
    :return: a matrix represented as a numpy ndarray with shape(sample_length, 39) and dtype=numpy.object.
            Each row is a single nucleotide annotation in format numpy.ndarray with shape(39,) and elements dtype=numpy.float16.
    """
    partial_sums = PARTIAL_SUMS[reference_genome]
    index = partial_sums[chromosome - 1]
    return annotations_db[nucleotides_db[index+start_position-1:index+end_position]]




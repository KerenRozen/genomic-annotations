import numpy as np

from consts import PARTIAL_SUMS_GH38, PARTIAL_SUMS_GH37

PARTIAL_SUMS = {
    37: PARTIAL_SUMS_GH37,
    38: PARTIAL_SUMS_GH38
}
NEGATIVE_STRAND = {
    37: 83,
    38:193
}
NUM_OF_REGULATORY_REGIONS= {
    37: 6,
    38: 5
}


def decode_strand(flag: int) -> str:
    """
    Receives sample's BAM FLAG.
    :return: sample's strand: '+' OR '-' based on:
            BAM flag: 0x10	SEQ being reverse complemented
    """
    strand_mask = 0b00010000
    reverse = int(bin(flag), 2) & strand_mask
    return '+' if reverse == 0 else '-'


def match_regulatory_regions_vector(db, reference_genome, chromosome: int, start_pos: int, end_pos: int, flag: int) -> list:
    """
    Receives a sample.
    :return: a list where each regulation is set to 1 if the sample within it and 0 if not.
    len(output) = 6 for hg37 and 5 for hg38.
    The regulatory regions classifications order:
    hg 37 ['enhancer', 'TF_binding_site', 'promoter_flanking_region', 'promoter', 'open_chromatin_region', 'CTCF_binding_site']
    hg 38 ['CTCF_binding_site', 'promoter', 'open_chromatin_region', 'TF_binding_site', 'enhancer']
    """
    strand = decode_strand(flag)
    partial_sums = PARTIAL_SUMS[reference_genome]
    neg_strand = NEGATIVE_STRAND[reference_genome]
    num_of_regulatory_regions = NUM_OF_REGULATORY_REGIONS[reference_genome]
    index = partial_sums[chromosome-1] if strand == '+' else partial_sums[chromosome+neg_strand]
    res = np.bitwise_or.reduce(db[index+start_pos-1 :index+end_pos])

    return np.unpackbits(np.array([res], dtype='>i4').view(np.uint8))[-num_of_regulatory_regions:].tolist()


def match_regulatory_regions_matrix(db, reference_genome, chromosome: int, start_pos: int, end_pos: int, flag: int) -> list:
    """
    Receives a sample.
    :return: a matrix with shape (sample_length, number_of_different_classifications) represented as a list of numpy arrays.
    The number_of_different_classifications for hg37 is 6 and for hg38 is 5.
    Each row is a single nucleotide annotation: a numpy ndarray in which each regulation is set to 1 if the sample within it and 0 if not.
    The regulatory regions classifications order:
    hg 37 ['enhancer', 'TF_binding_site', 'promoter_flanking_region', 'promoter', 'open_chromatin_region', 'CTCF_binding_site']
    hg 38 ['CTCF_binding_site', 'promoter', 'open_chromatin_region', 'TF_binding_site', 'enhancer']
    """
    strand = decode_strand(flag)
    partial_sums = PARTIAL_SUMS[reference_genome]
    neg_strand = NEGATIVE_STRAND[reference_genome]
    num_of_regulatory_regions = NUM_OF_REGULATORY_REGIONS[reference_genome]
    index = partial_sums[chromosome-1] if strand == '+' else partial_sums[chromosome+neg_strand]

    return [np.unpackbits(np.array([nuc], dtype='>i4').view(np.uint8))[-num_of_regulatory_regions:] for nuc in db[index+start_pos-1 :index+end_pos]]
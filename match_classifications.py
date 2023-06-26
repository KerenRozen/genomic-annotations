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

def decode_strand(flag: int) -> str:
    """
    Receives sample's BAM FLAG.
    :return: sample's strand: '+' OR '-' based on:
            BAM flag: 0x10	SEQ being reverse complemented
    """
    strand_mask = 0b00010000
    reverse = int(bin(flag), 2) & strand_mask
    return '+' if reverse == 0 else '-'


def match_classifications_vector(db, reference_genome, chromosome: int, start_pos: int, end_pos: int, flag: int) -> list:
    """
    Receives a sample.
    :return: a 25-long list where each classification is set to 1 if the sample within it and 0 if not.
    The classifications order: ['scaffold', 'pseudogene', 'lnc_RNA', 'ncRNA', 'unconfirmed_transcript',
                                'V_gene_segment', 'biological_region', 'snRNA', 'D_gene_segment', 'five_prime_UTR',
                                'pseudogenic_transcript', 'gene', 'mRNA', 'scRNA', 'snoRNA', 'tRNA', 'J_gene_segment',
                                'ncRNA_gene', 'exon', 'rRNA', 'miRNA', 'three_prime_UTR', 'transcript',
                                'C_gene_segment', 'CDS']

    """
    partial_sums = PARTIAL_SUMS[reference_genome]
    neg_strand = NEGATIVE_STRAND[reference_genome]
    strand = decode_strand(flag)
    index = partial_sums[chromosome-1] if strand == '+' else partial_sums[chromosome+neg_strand]
    res = np.bitwise_or.reduce(db[index+start_pos-1 :index+end_pos])

    return np.unpackbits(np.array([res], dtype='>i4').view(np.uint8))[-25:].tolist()


def match_classifications_matrix(db, reference_genome, chromosome: int, start_pos: int, end_pos: int, flag: int) -> list:
    """
    Receives a sample.
    :return: a matrix with shape (sample_length, 25). represented as a nested lists.
    Each row is a single nucleotide annotation: 25-long list where each classification is set to 1 if the sample within it and 0 if not.
    The classifications order: ['scaffold', 'pseudogene', 'lnc_RNA', 'ncRNA', 'unconfirmed_transcript',
                                'V_gene_segment', 'biological_region', 'snRNA', 'D_gene_segment', 'five_prime_UTR',
                                'pseudogenic_transcript', 'gene', 'mRNA', 'scRNA', 'snoRNA', 'tRNA', 'J_gene_segment',
                                'ncRNA_gene', 'exon', 'rRNA', 'miRNA', 'three_prime_UTR', 'transcript',
                                'C_gene_segment', 'CDS']

    """
    partial_sums = PARTIAL_SUMS[reference_genome]
    neg_strand = NEGATIVE_STRAND[reference_genome]
    strand = decode_strand(flag)
    index = partial_sums[chromosome-1] if strand == '+' else partial_sums[chromosome+neg_strand]

    return [np.unpackbits(np.array([nuc], dtype='>i4').view(np.uint8))[-25:].tolist() for nuc in db[index+start_pos-1 :index+end_pos]]




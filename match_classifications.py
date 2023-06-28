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
NUM_OF_CLASSIFICATIONS = {
    37: 32,
    38: 25
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
    :return: a list where each classification is set to 1 if the sample within it and 0 if not.
    The number_of_different_classifications for hg37 is 32 and for hg38 is 25.
    The classifications order:
    hg37 ['RNA', 'snoRNA_gene', 'lincRNA', 'VD_gene_segment', 'NMD_transcript_variant', 'exon',
         'miRNA', 'biological_region', 'mRNA', 'snRNA', 'rRNA_gene', 'miRNA_gene', 'gene',
         'rRNA', 'nc_primary_transcript', 'CDS', 'snoRNA', 'processed_pseudogene',
         'V_gene_segment', 'three_prime_UTR', 'mt_gene', 'processed_transcript', 'pseudogene',
         'snRNA_gene', 'lincRNA_gene', 'pseudogenic_transcript', 'J_gene_segment', 'supercontig',
         'C_gene_segment', 'aberrant_processed_transcript', 'transcript', 'five_prime_UTR']
    hg38 ['scaffold', 'pseudogene', 'lnc_RNA', 'ncRNA', 'unconfirmed_transcript',
        'V_gene_segment', 'biological_region', 'snRNA', 'D_gene_segment', 'five_prime_UTR',
        'pseudogenic_transcript', 'gene', 'mRNA', 'scRNA', 'snoRNA', 'tRNA', 'J_gene_segment',
        'ncRNA_gene', 'exon', 'rRNA', 'miRNA', 'three_prime_UTR', 'transcript',
        'C_gene_segment', 'CDS']

    """
    partial_sums = PARTIAL_SUMS[reference_genome]
    neg_strand = NEGATIVE_STRAND[reference_genome]
    num_of_classifications = NUM_OF_CLASSIFICATIONS[reference_genome]
    strand = decode_strand(flag)
    index = partial_sums[chromosome-1] if strand == '+' else partial_sums[chromosome+neg_strand]
    res = np.bitwise_or.reduce(db[index+start_pos-1 :index+end_pos])

    return np.unpackbits(np.array([res], dtype='>i4').view(np.uint8))[-num_of_classifications:].tolist()


def match_classifications_matrix(db, reference_genome, chromosome: int, start_pos: int, end_pos: int, flag: int):
    """
    Receives a sample.
    :return: a matrix with shape (sample_length, number_of_different_classifications) represented as a list of numpy arrays.
    The number_of_different_classifications for hg37 is 32 and for hg38 is 25.
    Each row is a single nucleotide annotation: a numpy ndarray where each classification is set to 1 if the sample within it and 0 if not.
    The classifications order:
    hg37 ['RNA', 'snoRNA_gene', 'lincRNA', 'VD_gene_segment', 'NMD_transcript_variant', 'exon',
         'miRNA', 'biological_region', 'mRNA', 'snRNA', 'rRNA_gene', 'miRNA_gene', 'gene',
         'rRNA', 'nc_primary_transcript', 'CDS', 'snoRNA', 'processed_pseudogene',
         'V_gene_segment', 'three_prime_UTR', 'mt_gene', 'processed_transcript', 'pseudogene',
         'snRNA_gene', 'lincRNA_gene', 'pseudogenic_transcript', 'J_gene_segment', 'supercontig',
         'C_gene_segment', 'aberrant_processed_transcript', 'transcript', 'five_prime_UTR']
    hg38 ['scaffold', 'pseudogene', 'lnc_RNA', 'ncRNA', 'unconfirmed_transcript',
        'V_gene_segment', 'biological_region', 'snRNA', 'D_gene_segment', 'five_prime_UTR',
        'pseudogenic_transcript', 'gene', 'mRNA', 'scRNA', 'snoRNA', 'tRNA', 'J_gene_segment',
        'ncRNA_gene', 'exon', 'rRNA', 'miRNA', 'three_prime_UTR', 'transcript',
        'C_gene_segment', 'CDS']
    """
    partial_sums = PARTIAL_SUMS[reference_genome]
    neg_strand = NEGATIVE_STRAND[reference_genome]
    num_of_classifications = NUM_OF_CLASSIFICATIONS[reference_genome]
    strand = decode_strand(flag)
    index = partial_sums[chromosome-1] if strand == '+' else partial_sums[chromosome+neg_strand]

    return [np.unpackbits(np.array([nuc], dtype='>i4').view(np.uint8))[-num_of_classifications:] for nuc in db[index+start_pos-1 :index+end_pos]]




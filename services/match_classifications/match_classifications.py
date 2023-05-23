import numpy as np

from services.consts import PARTIAL_SUMS_GH38


def decode_strand(flag: int) -> str:
    """
    Receives sample's BAM FLAG.
    :return: sample's strand: '+' OR '-' based on:
            BAM flag: 0x10	SEQ being reverse complemented
    """
    strand_mask = 0b00010000
    reverse = int(bin(flag), 2) & strand_mask
    return '+' if reverse == 0 else '-'


def match_classifications(db, chromosome: int, start_pos: int, end_pos: int, flag: int) -> list:
    """
    Receives a sample.
    :return: a 26-long list where each classification is set to 1 if the sample within it and 0 if not.
    The classifications order: ['scaffold', 'pseudogene', 'lnc_RNA', 'chromosome', 'ncRNA', 'unconfirmed_transcript',
                                'V_gene_segment', 'biological_region', 'snRNA', 'D_gene_segment', 'five_prime_UTR',
                                'pseudogenic_transcript', 'gene', 'mRNA', 'scRNA', 'snoRNA', 'tRNA', 'J_gene_segment',
                                'ncRNA_gene', 'exon', 'rRNA', 'miRNA', 'three_prime_UTR', 'transcript',
                                'C_gene_segment', 'CDS']

    """
    strand = decode_strand(flag)
    index = PARTIAL_SUMS_GH38[chromosome-1] if strand == '+' else PARTIAL_SUMS_GH38[chromosome+83]
    res = np.bitwise_or.reduce(db[index+start_pos-1 :index+end_pos])

    return np.unpackbits(np.array([res], dtype='>i4').view(np.uint8))[-26:].tolist()


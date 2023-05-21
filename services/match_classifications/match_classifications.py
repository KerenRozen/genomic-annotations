import numpy as np

from services.consts import PARTIAL_SUMS


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
    :return: a 32-long list where each classification is set to 1 if the sample within it and 0 if not.
    The classifications order:  ['RNA', 'snoRNA_gene', 'lincRNA', 'VD_gene_segment', 'NMD_transcript_variant', 'exon',
                                 'miRNA', 'biological_region', 'mRNA', 'snRNA', 'rRNA_gene', 'miRNA_gene', 'gene',
                                 'rRNA', 'nc_primary_transcript', 'CDS', 'snoRNA', 'processed_pseudogene',
                                 'V_gene_segment', 'three_prime_UTR', 'mt_gene', 'processed_transcript', 'pseudogene',
                                 'snRNA_gene', 'lincRNA_gene', 'pseudogenic_transcript', 'J_gene_segment', 'supercontig',
                                 'C_gene_segment', 'aberrant_processed_transcript', 'transcript', 'five_prime_UTR']
    """
    strand = decode_strand(flag)
    index = PARTIAL_SUMS[chromosome-1] if strand == '+' else PARTIAL_SUMS[chromosome+83]
    res = np.bitwise_or.reduce(db[index+start_pos :index+end_pos+1])

    return [(res >> i) & 1 for i in range(31, -1, -1)]
    #return res


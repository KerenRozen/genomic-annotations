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


def match_regulations(db, chromosome: int, start_pos: int, end_pos: int, flag: int) -> list:
    """
    Receives a sample.
    :return: a 5-long list where each regulation is set to 1 if the sample within it and 0 if not.
    The regulations order: ['CTCF_binding_site', 'promoter', 'open_chromatin_region', 'TF_binding_site', 'enhancer']
    """
    strand = decode_strand(flag)
    index = PARTIAL_SUMS_GH38[chromosome-1] if strand == '+' else PARTIAL_SUMS_GH38[chromosome+193]
    res = np.bitwise_or.reduce(db[index+start_pos-1 :index+end_pos])

    return np.unpackbits(np.array([res], dtype='>i4').view(np.uint8))[-5:].tolist()


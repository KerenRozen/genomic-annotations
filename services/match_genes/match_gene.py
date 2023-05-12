import logging
from base64 import b64decode
from functools import lru_cache

from bitarray import bitarray

from services.consts import CHROMOSOMES_INDEX_PATH
from services.common.compressed_json import load_json


@lru_cache(None)
def read_index():
    """
    Read the dictionary from the disk and prepare it to be used by decoding the base64 encoded bitarrays
    """
    logging.info('Loading index')
    index = load_json(CHROMOSOMES_INDEX_PATH)
    for key, value in index.items():
        index[key] = bitarray()
        index[key].frombytes(b64decode(value.encode('ascii')))
    logging.info('Finished loading')
    return index


def decode_strand(flag: int) -> str:
    """
    Receives sample's BAM FLAG.
    :return: sample's strand: '+' OR '-' based on:
            BAM flag: 0x10	SEQ being reverse complemented
    """
    strand_mask = 0b00010000
    reverse = int(bin(flag), 2) & strand_mask
    return '+' if reverse == 0 else '-'


def sample_matches_any_gene(chromosome: int, start_pos: int, end_pos: int, flag: int) -> bool:
    """
    Receives parameters of a sample.
    :return: Whether the sample intersects with any of the genes in the index.
    """
    strand = decode_strand(flag)
    index = read_index()
    return index.get(chromosome, bitarray())[start_pos:end_pos+1].any() if strand == '+'\
        else index.get(-chromosome, bitarray())[start_pos:end_pos+1].any()

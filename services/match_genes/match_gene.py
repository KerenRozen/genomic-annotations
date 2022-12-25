import logging
from base64 import b64decode
from functools import lru_cache

from bitarray import bitarray

from services.consts import CHROMOSOMES_INDEX_PATH
from services.match_genes.compressed_json import load_json


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


def sample_matches_any_gene(chromosome: int, start_pos: int, end_pos: int) -> bool:
    """
    Receives parameters of a sample.
    :return: Whether the sample intersects with any of the genes in the index.
    """
    index = read_index()
    return index.get(chromosome, bitarray())[start_pos:end_pos+1].any()

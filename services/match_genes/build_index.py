"""
Creating index to match samples with chromosome genes.
Index is of type dict[int, bitarray].
The keys of the index are chromosomes (e.g. 1, 43, 74...)
For each key the value is a bitarray with the length of the chromosome.
In the bitarray, all bits are 0 except for bits inside genes of this chromosome (which are 1).
"""

from bitarray import bitarray
from base64 import b64encode

from services.consts import CHROMOSOMES_INDEX_PATH, RAW_LENGTHS_PATH, RAW_GENES_PATH
from services.match_genes.compressed_json import dump_json


def _init_index():
    genome = {}
    with open(RAW_LENGTHS_PATH, 'r') as f1:
        for row in f1:
            chromosome, length = row.split()
            chromosome, length = int(chromosome), int(length)
            genome[chromosome] = bitarray(length)
            genome[chromosome].setall(0)

    return genome


def _fill_index(index):
    with open(RAW_GENES_PATH, 'r') as f2:
        next(f2)  # Skip header
        for row in f2:
            _, chromosome, start, end = row.split()
            chromosome, start, end = int(chromosome), int(start), int(end)
            index[chromosome][start:end+1] = 1


def _write_index(index):
    """
    Prepare the index dictionary to be written on disk by decoding the bitarrays to base64 text
    Dump the index using compressed_json
    """
    for key, value in index.items():
        index[key] = b64encode(value.tobytes()).decode('ascii')
    dump_json(CHROMOSOMES_INDEX_PATH, index)


def build_index():
    index = _init_index()
    _fill_index(index)
    _write_index(index)

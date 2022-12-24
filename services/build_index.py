import json
import zlib

from bitarray import bitarray
from base64 import b64encode


RAW_GENES_PATH = 'C:/Users/keren/PycharmProjects/project/genes.csv'
RAW_LENGTHS_PATH = 'C:/Users/keren/PycharmProjects/project/GRch37_chrom_lengths.tsv'
CHROMOSOMES_BITARRAY_PATH = 'C:/Users/keren/PycharmProjects/project/chromosomes_index.json'


def _parse_row(row):
    _, chromosome, start, end = row.split()
    chromosome, start, end = int(chromosome), int(start), int(end)
    return chromosome, start, end


def _init_index(chrom_lengths_path):
    genome = {}
    with open(chrom_lengths_path, 'r') as f1:
        for row in f1:
            chromosome, length = row.split()
            chromosome, length = int(chromosome), int(length)
            genome[chromosome] = bitarray(length)
            genome[chromosome].setall(0)

    return genome


def _fill_index(index, genes_csv_path):
    with open(genes_csv_path, 'r') as f2:
        next(f2)  # Skip header
        for row in f2:
            chromosome, start, end = _parse_row(row)
            index[chromosome][start:end+1] = 1


def _dump_index(index, json_path):
    for key, value in index.items():
        index[key] = b64encode(value.tobytes()).decode('ascii')

    with open(json_path, 'wb') as f:
        f.write(
            zlib.compress(
                json.dumps(index).encode('ascii')
            )
        )


def build_chrom_bitarray_index(chromosomes_lengths_path, raw_genes_path, index_path):
    # create key-val dict where keys are chromosomes and values are bitarray in chromosome's length
    index = _init_index(chromosomes_lengths_path)
    _fill_index(index, raw_genes_path)
    _dump_index(index, index_path)


if __name__ == '__main__':
    build_chrom_bitarray_index(RAW_LENGTHS_PATH, RAW_GENES_PATH, CHROMOSOMES_BITARRAY_PATH)

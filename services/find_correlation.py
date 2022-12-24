import json
import zlib
from base64 import b64decode
import time
from functools import lru_cache

from bitarray import bitarray
from tqdm import tqdm

CHROMOSOMES_BITARRAY_PATH = 'C:/Users/keren/PycharmProjects/project/chromosomes_index.json'
TEST_SAMPLES = 'C:/Users/keren/PycharmProjects/project/example_read_file.tsv'


def _parse_row(row):
    _, chromosome, start, end = row.split()
    chromosome, start, end = int(chromosome), int(start), int(end)
    return chromosome, start, end


@lru_cache(None)
def _load_index(json_path):
    print('Loading index')
    with open(json_path, 'rb') as f:
        index = json.loads(
            zlib.decompress(
                f.read()
            ).decode('ascii'),
            object_hook=lambda d: {int(k): v for k, v in d.items()})
    for key, value in index.items():
        index[key] = bitarray()
        index[key].frombytes(b64decode(value.encode('ascii')))
    print('Finished loading')
    return index


def does_gene_correlate(chromosome, start_pos, end_pos):
    index = _load_index(CHROMOSOMES_BITARRAY_PATH)
    return index.get(chromosome, bitarray())[start_pos:end_pos].any()  # Testing for any correlation between index and sample


if __name__ == '__main__':
    _load_index(CHROMOSOMES_BITARRAY_PATH)

    with open(TEST_SAMPLES, 'r') as f:
        next(f)
        while f:
            try:
                start_time = time.time()
                success = 0
                for i in tqdm(range(100000)):
                    row = next(f)
                    chromosome, start, end = _parse_row(row)
                    if chromosome != 100:
                        success += int(does_gene_correlate(chromosome, start, end))
                print(f'Time took {(time.time() - start_time)}', )
                print(success, 'Samples were success')
            except StopIteration:
                break



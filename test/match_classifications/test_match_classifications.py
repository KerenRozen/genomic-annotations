import time

import numpy as np
import pytest

from services.consts import TEST_SAMPLES, CLASSIFICATIONS_DB
from services.match_classifications.match_classifications import match_classifications


def test__runtime():
    db = np.load(CLASSIFICATIONS_DB)

    with open(TEST_SAMPLES, 'r') as f:
        next(f)  # Skip header
        row = next(f)
        rows = 1
        start_time = time.time()
        while row:
            _, chromosome, start, end, flag = row.split()
            chromosome, start, end, flag = int(chromosome), int(start), int(end), int(flag)
            res = match_classifications(db, chromosome, start, end, flag)
            rows += 1
            row = next(f, None)
        diff = time.time() - start_time
        assert(100000 * diff / rows <= 1)



@pytest.mark.parametrize(
    ['chromosome', 'start', 'end', 'strand', 'expected_result'],
    [
        (1, 235448349, 235448440, 15, [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1]),
        (10, 112445437, 112445650, 16, [0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0]),
        (11, 407315, 414923, 16, [0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1]),
        (12, 45050420, 45051603, 15, [0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0]),
        (14, 21269298, 21269410, 16, [0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0]),
        (3, 100, 300, 15, [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
        (18, 641327, 641541, 15, [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1]),
        (23, 51340284, 51340402, 16, [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0]),
        (68, 300, 800, 15, [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    ], ids=[
        "1.a 'lnc_RNA', 'mRNA', 'gene', 'exon', 'three_prime_UTR', 'CDS'",
        "1.b 'lnc_RNA', 'five_prime_UTR', 'mRNA', 'gene', 'exon'",
        "1.c 'transcript', 'lnc_RNA', 'five_prime_UTR', 'mRNA', 'biological_region', 'gene', 'exon', 'CDS'",
        "1.d 'exon', 'biological_region', 'lnc_RNA', 'ncRNA_gene'",
        "1.e 'lnc_RNA', 'five_prime_UTR', 'mRNA', 'biological_region', 'gene', 'exon'",
        "1.f nothing",
        "1.g 'exon', 'mRNA', 'CDS', 'gene'",
        "2.a 'exon', 'lnc_RNA', 'ncRNA_gene'",
        "2.b 'scaffold'"
    ])
def test__assert_classifications(chromosome, start, end, strand, expected_result):
    db = np.load(CLASSIFICATIONS_DB)
    actual_result = match_classifications(db, chromosome, start, end, strand)
    assert actual_result == expected_result, f"expected {expected_result}, but got {actual_result}"


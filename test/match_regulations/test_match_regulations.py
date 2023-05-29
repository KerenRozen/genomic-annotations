import time

import numpy as np
import pytest

from services.consts import TEST_SAMPLES, REGULATIONS_DB
from services.match_regulations.match_regulations import match_regulations


def test__runtime():
    db = np.load(REGULATIONS_DB)

    with open(TEST_SAMPLES, 'r') as f:
        next(f)  # Skip header
        row = next(f)
        rows = 1
        start_time = time.time()
        while row:
            _, chromosome, start, end, flag = row.split()
            chromosome, start, end, flag = int(chromosome), int(start), int(end), int(flag)
            res = match_regulations(db, chromosome, start, end, flag)
            rows += 1
            row = next(f, None)
        diff = time.time() - start_time
        assert(100000 * diff / rows <= 1)



@pytest.mark.parametrize(
    ['chromosome', 'start', 'end', 'strand', 'expected_result'],
    [
        (1, 235448349, 235448440, 15, [0, 0, 0, 0, 0]),
        (10, 112445437, 112445650, 16, [0, 1, 0, 0, 0]),
        (11, 407315, 414923, 16, [1, 1, 0, 0, 0]),
        (12, 45050420, 45051603, 15, [0, 1, 0, 0, 0]),
        (14, 21269298, 21269410, 16, [0, 1, 0, 0, 0]),
        (3, 100, 300, 15, [0, 0, 0, 0, 0]),
        (18, 641327, 641541, 15, [0, 0, 0, 0, 0]),
        (23, 51340284, 51340402, 16, [0, 0, 0, 0, 0]),
        (68, 300, 800, 15, [0, 0, 0, 0, 0]),
        (18, 35116801, 35120999, 1, [0, 0, 0, 0, 1]),
        (8, 37967115, 37967453, 83, [1, 0, 0, 1, 0]),
        (23, 90438125, 90438335, 196, [0, 0, 1, 0, 0])
    ], ids=[
        "1.a nothing",
        "1.b 'promoter'",
        "1.c 'CTCF_binding_site', 'promoter'",
        "1.d 'promoter'",
        "1.e 'promoter'",
        "1.f nothing",
        "1.g nothing",
        "2.a nothing",
        "2.b nothing",
        "2.c 'enhancer'",
        "2.d 'CTCF_binding_site', 'TF_binding_site'",
        "2.e 'open_chromatin_region'"
    ])
def test__assert_classifications(chromosome, start, end, strand, expected_result):
    db = np.load(REGULATIONS_DB)
    actual_result = match_regulations(db, chromosome, start, end, strand)
    assert actual_result == expected_result, f"expected {expected_result}, but got {actual_result}"


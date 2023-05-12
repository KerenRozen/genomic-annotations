import time

import pytest

from services.common.compressed_json import load_json
from services.consts import REGULATORY_REGIONS_DB, TEST_SAMPLES
from services.regulatory_regions.match_regulatory_regions import match_regulatory_regions


def test__runtime():
    db = load_json(REGULATORY_REGIONS_DB)

    with open(TEST_SAMPLES, 'r') as f:
        next(f)  # Skip header
        row = next(f)
        rows = 1
        start_time = time.time()
        while row:
            _, chromosome, start, end, flag = row.split()
            chromosome, start, end, flag = int(chromosome), int(start), int(end), int(flag)
            res = match_regulatory_regions(db, chromosome, start, end, flag)
            rows += 1
            row = next(f, None)
        diff = time.time() - start_time
        assert (100000 * diff / rows <= 1)


@pytest.mark.parametrize(
    ['chromosome', 'start', 'end', 'strand', 'expected_result'],
    [
        (51, 101295, 101509, 15, [0, 0, 0, 1, 0, 0]),
        (51, 18915, 19219, 1, [0, 0, 0, 0, 1, 0]),
        (24 ,100023401, 100024400, 1, [1, 0, 0, 0, 0, 0]),
        (24, 10003801, 10004400, 8, [0, 0, 0, 0, 0, 1]),
        (24,  100074201, 100076800, 9, [0, 1, 0, 0, 0, 0]),
        (24, 100000602, 100003999, 15, [0, 0, 1, 0, 0, 0]),
        (24, 100, 250, 16, [0, 0, 0, 0, 0, 0]),
        (24, 69674800, 69881201, 1, [1, 1, 0, 1, 0, 0])
    ], ids=[
        "1.a open_chromatin_region",
        "1.b TF_binding_site",
        "1.c CTCF_binding_site",
        "1.d enhancer",
        "1.e promoter",
        "1.f promoter_flanking_region",
        "1.g nothing",
        "1.h CTCF_binding_site, open_chromatin_region, promoter"
    ])
def test__assert_regulatory_regions(chromosome, start, end, strand, expected_result):
    db = load_json(REGULATORY_REGIONS_DB)
    actual_result = match_regulatory_regions(db, chromosome, start, end, strand)
    assert actual_result == expected_result, f"expected {expected_result}, but got {actual_result}"
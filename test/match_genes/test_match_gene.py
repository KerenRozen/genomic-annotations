"""
1. Test that match gene works for all expected inputs
    1.a. Complete overlap
    1.a0. (Mukhal)
    1.b Overlap from right
    1.c Overlap from left
    1.d No overlap
2. Test that for 100k samples takes less than 1s
"""
import time

import pytest

from services.consts import TEST_SAMPLES
from services.match_genes.match_gene import read_index, sample_matches_any_gene


def test__assert_mass_run_time_and_success_rate():
    read_index()  # For cache
    with open(TEST_SAMPLES, 'r') as f:
        next(f)  # Skip header
        start_time = time.time()
        row = next(f)
        rows = 1
        success = 0
        while row:
            _, chromosome, start, end, flag = row.split()
            chromosome, start, end, flag = int(chromosome), int(start), int(end), int(flag)
            #if chromosome != 100:
            success += int(sample_matches_any_gene(chromosome, start, end, flag))
            rows += 1
            row = next(f, None)
        diff = time.time() - start_time
        assert 100000 * diff / rows <= 1
        #assert success == 445435


@pytest.mark.parametrize(
    ['chromosome', 'start', 'end', 'strand', 'expected_result'],
    [
        (1, 69091, 70008, 15, True),
        (1, 69091, 70008, 16, False),
        (1, 69900, 70000, 15, True),
        (1, 69900, 70000, 16, False),
        (1, 69500, 70100, 15, True),
        (1, 69500, 70100, 16, False),
        (1, 69000, 69500, 15, True),
        (1, 69000, 69500, 16, False),
        (1, 71000, 71200, 15, False),
        (1, 71000, 71200, 16, False),
        (1, 249250600, 250250680, 16, False),
        (1, -10, 100, 15, False),
        (1, 134751, 134901, 16, True),
        (1,	367490, 367640,	15, True),
        (1, 140339, 140489, 16, True),
        (1, 70008, 70158, 15, True)


    ], ids=[
        "1.a Complete overlap",
        "1.a0. Wrong strand complete overlap",
        "1.b (Mukhal)",
        "1.b0. (Mukhal) wrong strand",
        "1.c Overlap from right",
        "1.c0. Overlap from right wrong strand",
        "1.d Overlap from left",
        "1.d0. Overlap from left wrong strand",
        "1.e No overlap on pos strand",
        "1.e0. No overlap on neg strand",
        "1.f End_pos is longer than chromosome's length",
        "1.g Start_pos < 0",
        "1.h Last base of interval match neg strand",
        "1.h0. Last base of interval match pos strand",
        "1.i First base of interval match neg strand",
        "1.i0. First base of interval match pos strand",
    ])
def test__assert_match_gene(chromosome, start, end, strand, expected_result):
    assert sample_matches_any_gene(chromosome, start, end, strand) is expected_result

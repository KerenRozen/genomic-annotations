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
            _, chromosome, start, end = row.split()
            chromosome, start, end = int(chromosome), int(start), int(end)
            if chromosome != 100:
                success += int(sample_matches_any_gene(chromosome, start, end))
            rows += 1
            row = next(f, None)
        diff = time.time() - start_time
        assert 100000 * diff / rows <= 1
        assert success == 445429


@pytest.mark.parametrize(
    ['chromosome', 'start', 'end', 'expected_result'],
    [
        (1, 69091, 70008, True),
        (1, 69900, 70000, True),
        (1, 69500, 70100, True),
        (1, 69000, 69500, True),
        (1, 71000, 71200, False)
    ], ids=[
        "1.a. Complete overlap",
        "1.a0. (Mukhal)",
        "1.b Overlap from right",
        "1.c Overlap from left",
        "1.d No overlap",
    ])
def test__assert_match_gene(chromosome, start, end, expected_result):
    assert sample_matches_any_gene(chromosome, start, end) is expected_result

import time

import pytest

from services.common.compressed_json import load_json
from services.consts import TEST_SAMPLES, CLASSIFICATIONS_DB
from services.match_classifications.match_classifications import match_classifications


def test__runtime():
    db = load_json(CLASSIFICATIONS_DB)

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
        assert(100000 * diff / rows <= 1.2)



@pytest.mark.parametrize(
    ['chromosome', 'start', 'end', 'strand', 'expected_result'],
    [
        (1, 69091, 70008, 15, [0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
        (1, 69091, 70008, 16, [0 for i in range(32)]),
        (1, 11870, 12030, 15, [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0]),
        (1, 134751, 134901, 16, [0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
        (1, 141474, 150000, 16, [0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0]),
        (1, 157800, 160500, 16, [0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0]),
        (1, 367650, 367660, 15, [0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1])
    ], ids=[
        "1.a gene",
        "1.b nothing",
        "1.c pseudogene, processed_transcript, exon",
        "1.d exon, mRNA, gene, three_prime_UTR",
        "1.e processed_transcript, transcript, exon, lincRNA",
        "1.f processed_transcript, exon, lincRNA, snRNA, snRNA_gene",
        "1.g lincRNA_gene, lincRNA, gene, mRNA, five_prime_UTR, exon, CDS"
    ])
def test__assert_classifications(chromosome, start, end, strand, expected_result):
    db = load_json(CLASSIFICATIONS_DB)
    actual_result = match_classifications(db, chromosome, start, end, strand)
    assert actual_result == expected_result, f"expected {expected_result}, but got {actual_result}"


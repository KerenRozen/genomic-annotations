import argparse
import gzip
import numpy as np

from pathlib import Path
from consts import CHROMOSOME_TO_INT_GH37, CHROMOSOME_TO_INT_GH38, PARTIAL_SUMS_GH37, \
    PARTIAL_SUMS_GH38

VALUES_LENGTH = {
    37: 7104162,
    38: 8011007
}
GENOME_LENGTH = {
    37: 3101804739,
    38: 3099411205,
}
PARTIAL_SUMS = {
    37: PARTIAL_SUMS_GH37,
    38: PARTIAL_SUMS_GH38
}
CHROMOSOME_TO_INT = {
    37: CHROMOSOME_TO_INT_GH37,
    38: CHROMOSOME_TO_INT_GH38
}


def init_db(input_file, reference_genome):
    """
    Using '2 lists' method.
    Input:  input_file- path to 'atlas.summary'
            reference_genome- 37 for gh37 or 38 for gh38
    """
    chromosome_int = CHROMOSOME_TO_INT[reference_genome]
    partial_sums = PARTIAL_SUMS[reference_genome]
    genome_length = GENOME_LENGTH[reference_genome]
    nucleotides_db = np.zeros(genome_length, dtype=np.uintc)
    values_db = np.empty(VALUES_LENGTH[reference_genome]+1, dtype=object)
    values_db[0] = np.zeros(39, dtype=np.float16)
    with gzip.open(Path(input_file), 'rt') as f:
        next(f) # Skip header
        for line_number, row in enumerate(f, start=1):
            fields = row.strip().split('\t')
            fields = list(map(lambda x: x.replace('NA', '0'), fields))    # convert NA to 0
            if fields[0] == 'chrM':   # convert M to MT to be able using the chromosome-int dict.
                fields[0] = 'chrMT'
            chromosome = chromosome_int[fields[0][3:]]
            start = int(fields[1])
            end = int(fields[2])
            values_db[line_number] = np.array(fields[5:], dtype=np.float16)
            index = partial_sums[chromosome-1]
            nucleotides_db[index+start-1:index+end] = line_number

    f.close()

    return nucleotides_db, values_db


def main():
    """
    Receives from user: path to atlas.summary file, path to save the DB, genome_reference (37 or 38)
    :return: Saves .npz file.
    """

    parser = argparse.ArgumentParser()
    parser.add_argument("methylation_file", help="path to atlas.summary file")
    parser.add_argument("db_file", help="path to save DB file")
    parser.add_argument("genome_reference", type=int, help="37 for gh37 or 38 for gh38")
    args = parser.parse_args()
    nucleotides_db, values_db = init_db(args.methylation_file, args.genome_reference)
    np.savez_compressed(args.db_file, nucleotides_db, values_db)


if __name__ == '__main__':
    main()
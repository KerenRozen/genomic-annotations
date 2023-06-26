import argparse
import gzip
import numpy as np
from pathlib import Path
import time
import sys

from consts import CHROMOSOME_TO_INT_GH38, PARTIAL_SUMS_GH38, CLASSIFICATIONS_ORDER_GH38, PARTIAL_SUMS_GH37, \
    CHROMOSOME_TO_INT_GH37, CLASSIFICATIONS_ORDER_GH37

GENOME_LENGTH = {
    37: 6203609478,
    38: 6198822410,
}
PARTIAL_SUMS = {
    37: PARTIAL_SUMS_GH37,
    38: PARTIAL_SUMS_GH38
}
CHROMOSOME_TO_INT = {
    37: CHROMOSOME_TO_INT_GH37,
    38: CHROMOSOME_TO_INT_GH38
}
CLASSIFICATIONS_ORDER = {
    37: CLASSIFICATIONS_ORDER_GH37,
    38: CLASSIFICATIONS_ORDER_GH38
}
NEGATIVE_STRAND = {
    37: 83,
    38:193
}
def init_classifications_db(input_file, reference_genome):
    chromosome_int = CHROMOSOME_TO_INT[reference_genome]
    partial_sums = PARTIAL_SUMS[reference_genome]
    genome_length = GENOME_LENGTH[reference_genome]
    classifications_order = CLASSIFICATIONS_ORDER[reference_genome]
    neg_strand = NEGATIVE_STRAND[reference_genome]
    db = np.zeros(genome_length, dtype = np.uintc)
    with gzip.open(Path(input_file), 'rt') as f:
        for row in f:
            if not row.startswith('#'):
                fields = row.strip().split('\t')
                seg_type = fields[2]
                if fields[2] != 'chromosome':
                    chromosome = chromosome_int[fields[0]]
                    start_pos = int(fields[3])
                    end_pos = int(fields[4])
                    strand = fields[6]
                    if strand == '+':
                        index = partial_sums[chromosome-1]
                        db[index+start_pos-1: index+end_pos] |= 2 ** classifications_order[seg_type]
                    elif strand == '-':
                        index = partial_sums[chromosome+neg_strand]
                        db[index+start_pos-1: index+end_pos] |= 2 ** classifications_order[seg_type]
                    else:
                        index1 = partial_sums[chromosome-1]
                        index2 = partial_sums[chromosome+neg_strand]
                        db[index1+start_pos-1: index1+end_pos] |= 2**classifications_order[seg_type]
                        db[index2+start_pos-1 :index2+end_pos] |= 2 ** classifications_order[seg_type]

    f.close()
    return db


def main():
    """
    Create classifications DB in .npy file from the gff file.
    The GFF path , the path where the DB will be stored and the genome_reference (37 or 38) are received as input from the user.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("genome_file", help="path to the genome file")
    parser.add_argument("db_file", help="path to save the DB file")
    parser.add_argument("genome_reference", type=int, help="37 for gh37 or 38 for gh38")
    args = parser.parse_args()
    start_time = time.time()
    db = init_classifications_db(args.genome_file, args.genome_reference)
    np.save(args.db_file, db)
    print(time.time()-start_time)
    print(sys.getsizeof(db))


if __name__ == '__main__':
    main()


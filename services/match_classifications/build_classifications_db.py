import argparse
import gzip
from pathlib import Path

from services.common.compressed_json import dump_json
from services.consts import CHROMOSOME_TO_INT, CLASSIFICATIONS_ORDER


def init_classifications_db(path, db_size=None):
    db = {}
    with gzip.open(Path(path), 'rt') as f:
        for row in f:
            if db_size and len(db) >= db_size:
                break

            if not row.startswith('#'):
                fields = row.strip().split('\t')
                seg_type = fields[2]
                if fields[2] != 'chromosome':
                    chromosome = CHROMOSOME_TO_INT[fields[0]]
                    start_pos = int(fields[3])
                    end_pos = int(fields[4])
                    strand = fields[6]
                    for i in range(start_pos, end_pos+1):
                        key = f'{chromosome}-{strand}-{i}'
                        if (db.get(key, 0) & 2**CLASSIFICATIONS_ORDER[seg_type]) == 0: # checking if this is the first time this nucleotide comes across the current classification
                            db[key] = db.get(key, 0)
                            db[key] += (2**CLASSIFICATIONS_ORDER[seg_type])
    f.close()
    return db


def main():
    """
    Create classifications DB in .json file from the gff file.
    The GFF path and the path where the DB will be stored are received as input from the user.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("genome_file", help="path to the genome file")
    parser.add_argument("classifications_db_file", help="path to the classifications DB will be saved (should end with .json)")
    args = parser.parse_args()
    db = init_classifications_db(args.genome_file, 10**7)
    dump_json(Path(args.classifications_db_file), db)


if __name__ == '__main__':
    main()


import gzip

from services.common.compressed_json import dump_json
from services.consts import CHROMOSOME_TO_INT, REGULATORY_REGIONS_ORDER, \
    REGULATORY_REGIONS, REGULATORY_REGIONS_DB

"""
599404 lines.
6 types of regulatory regions: ['CTCF_binding_site', 'promoter', 'promoter_flanking_region', 'open_chromatin_region', 'TF_binding_site', 'enhancer']
No strand info.  
{'CTCF_binding_site': 82274800, 'promoter': 71073936, 'promoter_flanking_region': 231927635, 'open_chromatin_region': 34109794, 'TF_binding_site': 11866041, 'enhancer': 100783808}
#Active nucleotides = 498248256
"""


def init_regulatory_regions(path, db_size=None):
    db = {}
    with gzip.open(REGULATORY_REGIONS, 'rt') as f:
        for row in f:
            if db_size and len(db) >= db_size:
                break

            fields = row.strip().split('\t')
            regulatory_type = fields[2]
            chromosome = CHROMOSOME_TO_INT[fields[0]]
            start_pos = int(fields[3])
            end_pos = int(fields[4])
            for i in range(start_pos, end_pos+1):
                key = f'{chromosome}-{i}'
                if db.get(key, 0) & 2 ** REGULATORY_REGIONS_ORDER[regulatory_type] == 0: # checking if this is the first time this nucleotide comes across the current classification
                    db[key] = db.get(key, 0)
                    db[key] += (2 ** REGULATORY_REGIONS_ORDER[regulatory_type])

    f.close()
    return db


def main():
    """
    Prepare the index dictionary to be written on disk by decoding the bitarrays to base64 text
    Dump the index using compressed_json
    """
    db = init_regulatory_regions(REGULATORY_REGIONS)
    dump_json(REGULATORY_REGIONS_DB, db)


if __name__ == '__main__':
    main()






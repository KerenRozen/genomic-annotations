import gzip

import numpy as np
#from bitarray import bitarray

from services.consts import CHROMOSOMES_LENGTHS_GH38, CH38_GENOME, CHROMOSOME_TO_INT_GH38, CLASSIFICATIONS_DB

#
# def active():
#     genome = {}
#     for k, v in CHROMOSOMES_LENGTHS_GH38.items():
#             genome[k] = bitarray(v)
#             genome[k].setall(0)
#             genome[-k] = bitarray(v)
#             genome[-k].setall(0)
#
#     with gzip.open(CH38_GENOME, 'rt') as f:
#         for row in f:
#             if not row.startswith('#'):
#                 fields = row.strip().split('\t')
#                 seg_type = fields[2]
#                 if fields[2] != 'chromosome':
#                     chromosome = CHROMOSOME_TO_INT_GH38[fields[0]]
#                     start= int(fields[3])
#                     end= int(fields[4])
#                     strand = fields[6]
#                     if strand == '+':
#                         genome[chromosome][start:end + 1] = 1
#                     elif strand == '-':
#                         genome[-chromosome][start:end + 1] = 1
#     f.close()
#     cnt = 0
#     for k, v in genome.items():
#         cnt += v.count(1)
#     return cnt, sum(CHROMOSOMES_LENGTHS_GH38.values())

#
# def verify_len():
#     d = {}
#     l = []
#     with open (RAW_LENGTHS_PATH, 'r') as f:
#         for row in f:
#             chromosome, length = row.split()
#             chromosome, length = int(chromosome), int(length)
#             d[chromosome] = length
#             l.append(length)
#     print(sum(d.values()))
#     print(sum(l))
#     print(d)
#     print(l)


if __name__ == '__main__':
    db = np.load(CLASSIFICATIONS_DB)
    active = np.count_nonzero(db)
    print("Number of active nucleotides=", active)
    print("Total nucleotides=", len(db))
    #
    # c, l = active()
    # print(c)
    # print(l)
    # #
    # # verify_len()
    # print(sum(CHROMOSOMES_LENGTHS.values()))
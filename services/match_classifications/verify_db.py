import gzip
import sys
import numpy as np

from services.consts import CHROMOSOMES_LENGTHS_GH38, GH38_GENOME, CHROMOSOME_TO_INT_GH38, CLASSIFICATIONS_DB, \
    CLASSIFICATIONS_ORDER_GH38, PARTIAL_SUMS_GH38


#
# def local_active():
#     genome = {}
#     #genome = {chr :np.zeros([CHROMOSOMES_LENGTHS_GH38[chr]]), -chr: np.zeros([CHROMOSOMES_LENGTHS_GH38[chr]])}
#     for k, v in CHROMOSOMES_LENGTHS_GH38.items():
#             genome[k] = bitarray(v)
#             genome[k].setall(0)
#             genome[-k] = bitarray(v)
#             genome[-k].setall(0)
#     with gzip.open(GH38_GENOME, 'rt') as f:
#         for row in f:
#             if not row.startswith('#'):
#                 fields = row.strip().split('\t')
#                 seg_type = fields[2]
#                 if fields[2] != 'chromosome':
#                     chromosome = CHROMOSOME_TO_INT_GH38[fields[0]]
#                     # if chromosome == chr:
#                     start= int(fields[3])
#                     end= int(fields[4])
#                     strand = fields[6]
#                     if strand == '+':
#                         genome[chromosome][start-1:end] = 1
#                     elif strand == '-':
#                         genome[-chromosome][start-1:end] = 1
#                     else:
#                         genome[chromosome][start-1:end] = 1
#                         genome[-chromosome][start-1:end] = 1
#
#     f.close()
#     #res = np.array([])
#     cnt = 0
#     for k, v in genome.items():
#         #res = np.concatenate((res, v))
#         cnt += v.count(1)
#     return cnt
#
#
# def numpy_active(chr):
#     partial = [0, CHROMOSOMES_LENGTHS_GH38[chr]]
#     db = np.zeros(CHROMOSOMES_LENGTHS_GH38[chr]*2, dtype = np.uintc) # Gh38
#     with gzip.open(GH38_GENOME, 'rt') as f:
#         for row in f:
#             if not row.startswith('#'):
#                 fields = row.strip().split('\t')
#                 seg_type = fields[2]
#                 if fields[2] != 'chromosome':
#                     chromosome = CHROMOSOME_TO_INT_GH38[fields[0]]
#                     if chromosome == chr:
#                         start_pos = int(fields[3])
#                         end_pos = int(fields[4])
#                         strand = fields[6]
#                         if strand == '+':
#                             index = partial[0]
#                             db[index+start_pos-1: index+end_pos] |= 2 ** CLASSIFICATIONS_ORDER_GH38[seg_type]
#                         elif strand == '-':
#                             index = partial[1]
#                             db[index+start_pos-1: index+end_pos] |= 2 ** CLASSIFICATIONS_ORDER_GH38[seg_type]
#                         else:
#                             index1 = partial[0]
#                             index2 =partial[1]
#                             db[index1+start_pos-1: index1+end_pos] |= 2**CLASSIFICATIONS_ORDER_GH38[seg_type]
#                             db[index2+start_pos-1 :index2+end_pos] |= 2 ** CLASSIFICATIONS_ORDER_GH38[seg_type]
#     f.close()
#     return db


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
    print('size:', sys.getsizeof(db))
    active = np.count_nonzero(db)
    print("Number of active nucleotides=", active)
    print("Total nucleotides=", len(db))
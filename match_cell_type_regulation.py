import numpy as np
import time
from services.consts import PARTIAL_SUMS_GH37, PARTIAL_SUMS_GH38
import random

PARTIAL_SUMS = {
    37: PARTIAL_SUMS_GH37,
    38: PARTIAL_SUMS_GH38
}


def annotate_cell_types_vector(nucleotides_db, annotations_db, reference_genome, chromosome, start_position, end_position):
    partial_sums = PARTIAL_SUMS[reference_genome]
    index = partial_sums[chromosome - 1]
    first_nuc = annotations_db[nucleotides_db[index+start_position-1]]
    last_nuc = annotations_db[nucleotides_db[index+end_position]]
    return [max(first_nuc[0], last_nuc[0]), max(first_nuc[1], last_nuc[1]), np.amax(np.array([first_nuc[2], last_nuc[2]]), axis=0)]


def annotate_cell_types_matrix(nucleotides_db, annotations_db, reference_genome, chromosome, start_position, end_position):
    partial_sums = PARTIAL_SUMS[reference_genome]
    index = partial_sums[chromosome - 1]
    result = annotations_db[nucleotides_db[index+start_position-1:index+end_position]]
    return result


if __name__ == '__main__':
    db = np.load('db38.npz', allow_pickle=True)
    nucleotides = db["arr_0"]
    annotations = db["arr_1"]
    start = time.time()
    for i in range(100000):
        item = [random.randint(1, 24), random.randint(1, 900000)]
        res = annotate_cell_types_vector(nucleotides, annotations,38, item[0], item[1], item[1]+150)
        #res = annotate_cell_types_matrix(nucleotides, annotations,37, 24, 70100, 70250)
        #res = annotate_cell_types_matrix(nucleotides, annotations, 38, item[0], item[1], item[1]+150)
    print(time.time()-start)
    print(type(res))
    print(res)

    #annotate_cell_types_matrix(nucleotides, annotations, 37, 24, 70100, 70120)


    # for n in [100,1000, 100000, 1000000]:
    #     cnt = 0
    #     for i in range(n):
    #         item = [random.randint(1, 24), random.randint(1, 900000), random.randint(1, 100000)]
    #         index = PARTIAL_SUMS_GH37[item[0]-1]
    #         start = index+item[1]-1
    #         end = start+150
    #         for j in range(start, end):
    #             nuc_j = nucleotides[j]
    #             first_val = annotations[nuc_j]
    #             nuc_j_plus1 = nucleotides[j+1]
    #             sec_val = annotations[nuc_j_plus1]
    #             if  first_val[0] != sec_val[0] or first_val[1] != sec_val[1] or not np.array_equal(first_val[2], sec_val[2]):
    #                 cnt+=1
    #                 break
    #     print(n , cnt, (cnt/n)*100, '%')

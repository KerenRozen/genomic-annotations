import gzip, bitarray
import sys

from services.consts import SEGWAY_ENCYCLOPEDIA_37, CHROMOSOMES_LENGTHS_GH37, CHROMOSOME_TO_INT_GH37, \
    SEGWAY_ENCYCLOPEDIA_38, CHROMOSOMES_LENGTHS_GH38, CHROMOSOME_TO_INT_GH38


def unzip_file(output_file, input_file):
    with gzip.open(input_file, 'rt') as f, open(output_file, 'w') as f2:
        for row in f:
            f2.write(row)

    f.close()
    f2.close()


def overlaps(input_file):
    ones = {}
    d = {}
    chromosomes = set()
    for k in CHROMOSOMES_LENGTHS_GH38:
        d[k] = bitarray.bitarray(CHROMOSOMES_LENGTHS_GH38[k])
        d[k].setall(0)

    with gzip.open(input_file, 'rt') as f:
        next(f)
        rows = 0
        minimum = sys.maxsize
        for row in f:
            rows+=1
            fields = row.strip().split('\t')
            #end_index =  min(fields[0].find(' ', 3), fields[0].find('_', 3))
            chromosome = CHROMOSOME_TO_INT_GH38[fields[0][3:]] if '_' not in fields[0] else CHROMOSOME_TO_INT_GH38[fields[0][3:fields[0].index('_')]]
            chromosomes.add(chromosome)
            start = int(fields[1])
            end = int(fields[2])
            #minimum = min(minimum, end-start+1)
            #
            ones.setdefault(chromosome, 0)
            ones[chromosome] += (end-start+1)
            # if d[chromosome][start-1:end].any():
            #     print(d[chromosome][start-1:end])
            #     print(chromosome, start, end)
            d[chromosome][start-1:end] = 1


    f.close()



    cntr = {}
    for k in d.keys():
        c = d[k].count(1)
        if c != 0:
            cntr[k] = c

    print(sum(cntr.values()))
    print(sum(ones.values()))
    #
    # for i in range(1,23):
    #     if ones[i] != cntr[i]:
    #         print(i, " in file:", ones[i], " in bitarray:", cntr[i])
    #         print('difference:', ones[i]-cntr[i])
    print(chromosomes)
    print(rows)


if __name__ == '__main__':
    #unzip_file('segway_encyclopedia_38', 'segway_encyclopedia_with_header_hg38.bed.gz')
    overlaps('segway_encyclopedia_with_header_hg38.bed.gz')

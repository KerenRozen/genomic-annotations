import math
import time
import pybedtools
import pysam
import matplotlib.pyplot as plt
from services.consts import READS, GENOME, TEST_SAMPLES, BEDTOOLS_TIMES, MY_ALGORITHM_TIMES
from services.match_genes.match_gene import sample_matches_any_gene, read_index


def check_bedtools_times():
    """
    Measures runtime of bedtools intersect method (between samples and genes)
    """
    genome = pybedtools.BedTool(GENOME)
    genes = genome.filter(lambda x: x[2] == 'gene')

    reads = pysam.AlignmentFile(READS, 'rb').head(1000)

    start_time = time.time()
    for read in reads:
        if not read.is_unmapped:
            chrom = read.reference_name
            start = read.reference_start
            end = read.reference_end
            name = read.query_name
            score = read.mapping_quality
            strand = '-' if read.is_reverse else '+'
            read= pybedtools.BedTool(f'{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}', from_string=True)
            read.intersect(genes, u=True, s=True)

    diff = time.time() - start_time
    print((1000, diff))


#### TEST MY CODE'S PERFORMANCE
def check_my_algorithm_times():
    """
    Measures runtime of my algorithm times (for intersection between samples and genes)
    """
    my_code_times = []
    cnt = [10,100,1000,10000,75000]

    read_index()
    # For cache
    for num in cnt:
        with open(TEST_SAMPLES, 'r') as f:
            next(f)  # Skip header
            row = next(f)
            start_time = time.time()
            for i in range(num):
                _, chromosome, start, end, flag = row.split()
                chromosome, start, end, flag = int(chromosome), int(start), int(end), int(flag)
                sample_matches_any_gene(chromosome, start, end, flag)
        diff = time.time() - start_time
        f.close()
        my_code_times.append((num, diff))
    return my_code_times


def plot_time_performances():
    samples = [10,100,1000,10000,75000]
    my_algorithm_times = list(map(lambda t: math.log2(t), MY_ALGORITHM_TIMES))
    bedtools_times = list(map(lambda t: math.log2(t), BEDTOOLS_TIMES))


    fig, ax = plt.subplots()

    #ax.plot(samples, my_algorithm_times, c='green', label='My algorithm')
    ax.plot(samples, bedtools_times, c='blue', label='bedtools')

    #plt.scatter(samples, my_algorithm_times, marker='o', color='green')
    plt.scatter(samples, bedtools_times, marker='o', color='blue')

    plt.xlabel('Number of samples')
    plt.ylabel('log2(runtime)')
    #plt.title('My algorithm vs bedtools runtimes')
    #plt.title('My algorithm runtimes')
    plt.title('bedtools runtimes')

    plt.legend()

    #plt.savefig()

    #plt.show()


if __name__ == '__main__':
    plot_time_performances()

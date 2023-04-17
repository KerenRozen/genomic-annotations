import csv
import gzip
import pybedtools
from services.match_genes.match_gene import sample_matches_any_gene

human_genome = 'Homo_sapiens.GRCh37.87.gff3.gz'
samples = 'SLX-11873.D707_D502.HFNWFBBXX.s_5.GRCh37.bwa2.bamdownsample001.sorted.bam'
pybedtools_intersect = 'bedtools_intersect_method_match_gene_results.csv'
my_code_intersect = 'my_algorithm_match_gene_results.csv'
test_samples = 'test_samples_with_flag.tsv'
bam_reads = 'SLX-11873.D707_D502.HFNWFBBXX.s_5.GRCh37.bwa2.bamdownsample001.sorted.bam'

CHROMOSOME_TO_INT = {'10': 10, '11': 11, '12': 12, '13': 13, '14': 14, '15': 15, '16': 16, '17': 17, '18': 18,
                     '19': 19, '1': 1, '20': 20, '21': 21, '22': 22, '2': 2, '3': 3, '4': 4, '5': 5, '6': 6,
                     '7': 7, '8': 8, '9': 9, 'MT': 23, 'X': 24, 'Y': 25, "GL000192.1": 26, "GL000225.1": 27,
                     "GL000194.1": 28, "GL000193.1": 29, "GL000200.1": 30, "GL000222.1": 31, "GL000212.1": 32,
                     "GL000195.1": 33, "GL000223.1": 34, "GL000224.1": 35, "GL000219.1": 36, "GL000205.1": 37,
                     "GL000215.1": 38, "GL000216.1": 39, "GL000217.1": 40, "GL000199.1": 41, "GL000211.1": 42,
                     "GL000213.1": 43, "GL000220.1": 44, "GL000218.1": 45, "GL000209.1": 46, "GL000221.1": 47,
                     "GL000214.1": 48, "GL000228.1": 49, "GL000227.1": 50, "GL000191.1": 51, "GL000208.1": 52,
                     "GL000198.1": 53, "GL000204.1": 54, "GL000233.1": 55, "GL000237.1": 56, "GL000230.1": 57,
                     "GL000242.1": 58, "GL000243.1": 59, "GL000241.1": 60, "GL000236.1": 61, "GL000240.1": 62,
                     "GL000206.1": 63, "GL000232.1": 64, "GL000234.1": 65, "GL000202.1": 66, "GL000238.1": 67,
                     "GL000244.1": 68, "GL000248.1": 69, "GL000196.1": 70, "GL000249.1": 71, "GL000246.1": 72,
                     "GL000203.1": 73, "GL000197.1": 74, "GL000245.1": 75, "GL000247.1": 76, "GL000201.1": 77,
                     "GL000235.1": 78, "GL000239.1": 79, "GL000210.1": 80, "GL000231.1": 81, "GL000229.1": 82,
                     "GL000226.1": 83, "GL000207.1": 84}

# a function to extract gene features from the genome file
def filter_genes_from_genome(genome):
    gene_list = []
    with gzip.open(genome, 'rt') as f:
        for line in f:
            if not line.startswith('#'):
                fields = line.strip().split('\t')
                if fields[2] == 'gene':
                    gene_list.append([CHROMOSOME_TO_INT[fields[0]], fields[3], fields[4], fields[6]])
    return gene_list


# a function to create test samples file in format: [chromosome, start_pos, end_pos, flag] from the BAM reads.
def create_samples_with_flag():
    mapped_mask = 0b00000100
    with open(test_samples, 'w', newline='') as file:
        writer = csv.writer(file, delimiter='\t')
        writer.writerow(['', 'chromosome', 'start_pos', 'end_pos', 'strand'])
        reads = pybedtools.BedTool(samples)
        for i, read in enumerate(reads):
            bam_flag = int(read.strand)
            mapped = int(bin(bam_flag), 2) & mapped_mask
            if mapped == 0:
                chromosome = CHROMOSOME_TO_INT[read.chrom]
                start = read.start
                end = read.end
                writer.writerow([i, chromosome, start, end, bam_flag])

    file.close()


# a function to create a csv file with the intersection reads calculated by bedtool.intersect
def find_genes_and_reads_intersection_with_pybedtools():
    genome = pybedtools.BedTool(human_genome)
    genes = genome.filter(lambda x: x[2] == 'gene')
    reads = pybedtools.BedTool(samples)
    intersecting_reads = reads.intersect(genes, u=True, s=True)

    mask = 0b00010000

    with open(pybedtools_intersect, 'w', newline='') as file:
        writer = csv.writer(file)
        for read in intersecting_reads:
            chromosome = CHROMOSOME_TO_INT[read.chrom]
            start = read.start
            end = read.end
            strand = int(read.strand)
            if ((int(bin(strand), 2)) & mask) != 0:
                strand = '-'
            else:
                strand = '+'
            writer.writerow([chromosome, start, end, strand])
    file.close()


# a function to create a csv file with the intersection reads calculated by my code.
def find_genes_and_reads_intersection_with_mycode():
    """
    The way to extract mapped, strand info based on bitwise AND on the read's BAM flag:
    BAM flag: 0x4	segment unmapped
    BAM flag: 0x10	SEQ being reverse complemented
    """
    reads = pybedtools.BedTool(samples)
    mapped_mask = 0b00000100
    strand_mask = 0b00010000

    with open(my_code_intersect, 'w', newline='') as file:
        writer = csv.writer(file)
        for read in reads:
            bam_flag = int(read.strand)
            mapped = int(bin(bam_flag), 2) & mapped_mask
            reverse = int(bin(bam_flag), 2) & strand_mask
            if mapped == 0:
                chromosome = CHROMOSOME_TO_INT[read.chrom]
                start = read.start
                end = read.end
                if reverse != 0:
                    strand = '-'
                else:
                    strand = '+'
                if sample_matches_any_gene(chromosome, start, end, strand):
                    writer.writerow([chromosome, start, end, strand])
        file.close()


if __name__ == '__main__':
    create_samples_with_flag()

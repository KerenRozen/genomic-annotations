from pathlib import Path


DATA_PATH = Path("../../data")
#RAW_GENES_PATH = DATA_PATH / 'genes.csv'
RAW_GENES_PATH = DATA_PATH / 'genes_with_strand.csv'
RAW_LENGTHS_PATH = DATA_PATH / 'GRch37_chrom_lengths.tsv'
CHROMOSOMES_INDEX_PATH = DATA_PATH / 'chromosomes_index.json'
#TEST_SAMPLES = DATA_PATH / 'test_samples.tsv'
TEST_SAMPLES = DATA_PATH / 'test_samples_with_flag.tsv'
READS = DATA_PATH / 'SLX-11873.D707_D502.HFNWFBBXX.s_5.GRCh37.bwa2.bamdownsample001.sorted.bam'
GENOME = DATA_PATH / 'Homo_sapiens.GRCh37.87.gff3.gz'
INTERSECT = DATA_PATH / 'bedtools_intersect_method_match_gene_results.csv'
MY_CODE = DATA_PATH / 'my_algorithm_match_gene_results.csv'

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

### Match gene Runtimes
BEDTOOLS_TIMES = [29.757471561431885, 30.437360048294067, 43.943015813827515, 106.90546178817749, 1528.7368199825287]
MY_ALGORITHM_TIMES = [0.000926971435546875, 0.0010945796966552734, 0.002707958221435547, 0.02183985710144043,
                 0.16087698936462402]
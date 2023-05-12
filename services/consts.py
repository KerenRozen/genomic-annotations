from pathlib import Path

DATA_PATH = Path("../../data")


RAW_LENGTHS_PATH = DATA_PATH / 'GRch37_chrom_lengths.tsv'
CHROMOSOMES_INDEX_PATH = DATA_PATH / 'chromosomes_index.json'
CLASSIFICATIONS_PATH = DATA_PATH / 'classifications_index.json'
REGULATORY_REGIONS = DATA_PATH / 'homo_sapiens.GRCh37.Regulatory_Build.regulatory_features.20201218.gff.gz'
REGULATORY_REGIONS_DB = DATA_PATH / 'regulatory_regions_db.json'
TEST_SAMPLES = DATA_PATH / 'test_samples_with_flag.tsv'
READS = DATA_PATH / 'SLX-11873.D707_D502.HFNWFBBXX.s_5.GRCh37.bwa2.bamdownsample001.sorted.bam'
GENOME = DATA_PATH / 'Homo_sapiens.GRCh37.87.gff3.gz'
CLASSIFICATIONS_DB = DATA_PATH / 'classifications_db.json'
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


CLASSIFICATIONS_ORDER = {'RNA': 31, 'snoRNA_gene': 30, 'lincRNA': 29, 'VD_gene_segment': 28, 'NMD_transcript_variant': 27,
                         'exon': 26, 'miRNA': 25, 'biological_region': 24, 'mRNA': 23, 'snRNA': 22, 'rRNA_gene': 21,
                         'miRNA_gene': 20, 'gene': 19, 'rRNA': 18, 'nc_primary_transcript': 17, 'CDS': 16, 'snoRNA': 15,
                         'processed_pseudogene': 14, 'V_gene_segment': 13, 'three_prime_UTR': 12, 'mt_gene': 11,
                         'processed_transcript': 10, 'pseudogene': 9, 'snRNA_gene': 8, 'lincRNA_gene': 7,
                         'pseudogenic_transcript': 6, 'J_gene_segment': 5, 'supercontig': 4, 'C_gene_segment': 3,
                         'aberrant_processed_transcript': 2, 'transcript': 1, 'five_prime_UTR': 0}


REGULATORY_REGIONS_ORDER = {'CTCF_binding_site': 5, 'promoter': 4, 'promoter_flanking_region':3,
                            'open_chromatin_region':2, 'TF_binding_site':1, 'enhancer':0}


CHROMOSOMES_LENGTHS ={1: 249250621, 10: 35534747, 11: 135006516, 12: 133851895, 13: 115169878, 14: 107349540, 15: 102531392,
                      16: 90354753, 17: 81195210, 18: 78077248, 19: 59128983, 2: 243199373, 20:	63025520, 21: 48129895,
                      22: 51304566, 3: 198022430, 4: 191154276, 5: 180915260, 6: 171115067,  7:	159138663,
                      8: 146364022, 9: 141213431, 23: 16569, 24: 155270560, 25:	59373566, 26: 547496, 27: 211173, 28: 191469,
                      29: 189789, 30: 187035, 31: 186861, 32: 186858, 33: 182896, 34: 180455, 35: 179693, 36: 179198,
                      37: 174588, 38: 172545, 39: 172294, 40: 172149, 41: 169874, 42: 166566, 43: 164239, 44: 161802,
                      45: 161147, 46: 159169, 47: 155397, 48: 137718, 49: 129120, 50: 128374, 51: 106433, 52: 92689,
                      53: 90085, 54: 81310, 55: 45941, 56: 45867, 57: 43691, 58: 43523, 59:	43341, 60: 42152, 61: 41934,
                      62: 41933, 63: 41001, 64:	40652, 65: 40531, 66: 40103, 67: 39939, 68:	39929, 69: 39786, 70: 38914,
                      71: 38502, 72: 38154, 73:	37498, 74: 37175, 75: 36651, 76: 36422, 77:	36148, 78: 34474, 79: 33824,
                      80: 27682, 81: 27386, 82:	19913, 83: 15008, 84: 4262}

TYPES_MASKS = {'RNA': '00000000000000000000000000000001', 'snoRNA_gene': '00000000000000000000000000000010',
               'lincRNA': '00000000000000000000000000000100', 'VD_gene_segment': '00000000000000000000000000001000',
               'NMD_transcript_variant': '00000000000000000000000000010000', 'exon': '00000000000000000000000000100000',
               'miRNA': '00000000000000000000000001000000', 'biological_region': '00000000000000000000000010000000',
               'mRNA': '00000000000000000000000100000000', 'snRNA': '00000000000000000000001000000000',
               'rRNA_gene': '00000000000000000000010000000000', 'miRNA_gene': '00000000000000000000100000000000',
               'gene': '00000000000000000001000000000000', 'rRNA': '00000000000000000010000000000000',
               'nc_primary_transcript': '00000000000000000100000000000000', 'CDS': '00000000000000001000000000000000',
               'snoRNA': '00000000000000010000000000000000', 'processed_pseudogene': '00000000000000100000000000000000',
               'V_gene_segment': '00000000000001000000000000000000', 'three_prime_UTR': '00000000000010000000000000000000',
               'mt_gene': '00000000000100000000000000000000', 'processed_transcript': '00000000001000000000000000000000',
               'pseudogene': '00000000010000000000000000000000', 'snRNA_gene': '00000000100000000000000000000000',
               'lincRNA_gene': '00000001000000000000000000000000', 'pseudogenic_transcript': '00000010000000000000000000000000',
               'J_gene_segment': '00000100000000000000000000000000', 'supercontig': '00001000000000000000000000000000',
               'C_gene_segment': '00010000000000000000000000000000', 'aberrant_processed_transcript': '00100000000000000000000000000000',
               'transcript': '01000000000000000000000000000000', 'five_prime_UTR': '10000000000000000000000000000000'}


### Match gene Runtimes
BEDTOOLS_TIMES = [29.757471561431885, 30.437360048294067, 43.943015813827515, 106.90546178817749, 1528.7368199825287]
MY_ALGORITHM_TIMES = [0.000926971435546875, 0.0010945796966552734, 0.002707958221435547, 0.02183985710144043,
                 0.16087698936462402]
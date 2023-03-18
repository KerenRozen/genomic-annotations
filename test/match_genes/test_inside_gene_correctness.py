import pysam
from pybedtools import BedTool
import pandas as pd
from services.consts import READS, CHROMOSOME_TO_INT, GENOME, INTERSECT, MY_CODE
from services.match_genes.match_gene import sample_matches_any_gene


def test__my_code_using_pysam():
    """
    Checking my code using pysam, by iterating through the reads.
    """
    bamfile = pysam.AlignmentFile(READS, 'rb')
    #print(bamfile.count()) ---> returns 77678

    success = 0
    unmapped = 0

    for read in bamfile:
        if read.is_unmapped:
            #print(f'Read {read.query_name} is unmapped')
            unmapped += 1
        else:
            chromosome = CHROMOSOME_TO_INT[bamfile.get_reference_name(read.reference_id)]
            start_pos = read.reference_start
            end_pos = read.reference_end
            success += int(sample_matches_any_gene(chromosome, start_pos, end_pos))

    bamfile.close()
    print('#Success: ', success)
    print('#Unmapped: ', unmapped)
    print('The BAM file contains 77678 reads')
    print(f'{success/77678}% of samples are inside genes')


def test__my_code_with_bedtools():
    """
    Checking my code using pybedtools intersect method, by iterating through the reads.
    """
    reads = BedTool(READS)
    success =  0

    for read in reads:
        chromosome = CHROMOSOME_TO_INT.get(read.chrom, -1)
        if chromosome == -1:
            continue
        else:
            start = read.start
            end = read.end
            if sample_matches_any_gene(chromosome, start, end):
                success += 1

    return success


def test__intersection_using_bedtools():
    """
    1. Finding the intersection between the reads & genes by using pybedtools.intersect method..
    2. Checking my code on every read among the intersection reads.
    """
    # Load the GFF3 file as a BedTool object
    genome = BedTool(GENOME)
    genes = genome.filter(lambda x: x[2] == 'gene')
    # Load the BAM file as a BedTool object
    reads = BedTool(READS)

    # Find intersecting reads with gene regions
    intersecting_reads = reads.intersect(genes, u=True)

    # Iterate over intersecting reads
    for read in intersecting_reads:
        chromosome = CHROMOSOME_TO_INT[read.chrom]
        start = read.start
        end = read.end

        assert sample_matches_any_gene(chromosome, start, end) == True


def test__find_diff():
    """
    This function print the samples that received different result from my code and from pybedtools.intersect method.
    (There's a difference only in a single read: (1, 171086959 ,171087109)
    """
    bedtools_samples = pd.read_csv(INTERSECT)

    my_code_samples = pd.read_csv(MY_CODE)
    bedtools_samples.columns = my_code_samples.columns = ['chromosome', 'start', 'end']
    merged = bedtools_samples.merge(my_code_samples, on=['chromosome', 'start', 'end'], how='outer', indicator=True)

    # Filter the merged dataframe to show only rows with differences
    diffs = merged[merged['_merge'] != 'both']

    # Add a column to show which dataframe each row came from
    diffs['From'] = diffs['_merge'].apply(lambda x: 'bedtools' if x == 'left_only' else 'my_code')

    # Drop the '_merge' column
    diffs.drop('_merge', axis=1, inplace=True)

    # Print the resulting dataframe showing the differences and which dataframe each row came from
    print(diffs)

import gzip

from services.consts import GH38_GENOME, CHROMOSOMES_LENGTHS_GH38, CHROMOSOME_TO_INT_GH38


def create_partial_sums():
    sorted_lengths = [value for _, value in sorted(CHROMOSOMES_LENGTHS_GH38.items())] * 2
    partial = [0]
    for i in range(len(sorted_lengths)):
        partial.append(partial[i] + sorted_lengths[i])
    print(partial)

def build_lengths_dict():
    chromosomes_to_int = {}
    chromosomes_legths = {}
    with open('../../data/Gh38_lengths.tsv', 'rt') as f:
        cnt = 1
        for row in f:
            fields = row.strip().split()
            chromosome = fields[1]
            chromosomes_to_int[chromosome] = cnt
            length = fields[3]
            chromosomes_legths[cnt] = int(length)
            cnt += 1
    print(chromosomes_to_int)
    print(chromosomes_legths)
    print(sum(chromosomes_legths.values()))


if __name__ == '__main__':
    #build_lengths_dict()
    #create_partial_sums()
    types = set()
    chrs = set()
    strands = set()
    with gzip.open(GH38_GENOME, 'rt') as f:
        for row in f:
            if not row.startswith('#'):
                fields = row.strip().split('\t')
                seg_type = fields[2]
                types.add(seg_type)
                if fields[2] != 'chromosome':
                    chromosome = CHROMOSOME_TO_INT_GH38[fields[0]]
                    chrs.add(chromosome)
                    strand = fields[6]
                    strands.add(strand)
    f.close()
    print(len(types))
    print(len(chrs))
    print(types)
    print(chrs)
    print(strands)

    # types = {'scaffold', 'pseudogene', 'lnc_RNA', 'chromosome', 'ncRNA', 'unconfirmed_transcript', 'V_gene_segment', 'biological_region', 'snRNA', 'D_gene_segment', 'five_prime_UTR', 'pseudogenic_transcript', 'gene', 'mRNA', 'scRNA', 'snoRNA', 'tRNA', 'J_gene_segment', 'ncRNA_gene', 'exon', 'rRNA', 'miRNA', 'three_prime_UTR', 'transcript', 'C_gene_segment', 'CDS'}

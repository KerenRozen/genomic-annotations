import gzip

import redis
from bitarray import bitarray

from services.consts import GENOME, CHROMOSOME_TO_INT, CHROMOSOMES_LENGTHS, TYPES_MASKS


def init_db():
    redis_db = redis.Redis()
    segments = {} #maps chromosomes to segments within them
    with gzip.open(GENOME, 'rt') as f:
        for row in f:
            if (not row.startswith('#')):
                fields = row.strip().split('\t')
                if fields[2] != 'chromosome':
                    chromosome = CHROMOSOME_TO_INT[fields[0]]
                    if chromosome not in segments:
                        segments[chromosome] = []
                    segments[chromosome].append((fields[2], fields[3], fields[4], fields[6]))
    f.close()
    for chrom in segments:
        segments[chrom].sort(key=lambda x: int(x[1])) # sorting segments by start position

    chromosomes = [i for i in range(1,85)]  # list of chromosomes in the genome
    for chrom in chromosomes:
        # Create a bitarray for each strand of the chromosome
        segments_forward = bitarray(32 * CHROMOSOMES_LENGTHS[chrom])
        segments_reverse = bitarray(32 * CHROMOSOMES_LENGTHS[chrom])
        segments_forward.setall(0)
        segments_reverse.setall(0)

        for segment in segments[chrom]:
            seg_type = segment[0]
            start = int(segment[1]) - 1  # GFF is 1-based, but Python is 0-based
            end = int(segment[2])
            strand = segment[3]

            # Set the relevant bits in the bitarray
            if strand == '+':
                for i in range(start, end):
                    idx = i * 32
                    a = segments_forward[idx:idx + 32]
                    b = bitarray(TYPES_MASKS[seg_type])
                    segments_forward[idx:idx + 32] = a | b
            elif strand == '-':
                for i in range(start, end):
                    idx = i * 32
                    a = segments_reverse[idx:idx + 32]
                    b = bitarray(TYPES_MASKS[seg_type])
                    segments_reverse[idx:idx + 32] = a | b

        # Upload the bitarrays to Redis

        for i in range(CHROMOSOMES_LENGTHS[chrom]):
            pos = f'{chrom}:+:{i+1}'
            val = segments_forward[i * 32:(i + 1) * 32].tobytes()
            redis_db.set(pos, val)
            pos = f'{chrom}:-:{i+1}'
            val = segments_reverse[i * 32:(i + 1) * 32].tobytes()
            redis_db.set(pos, val)
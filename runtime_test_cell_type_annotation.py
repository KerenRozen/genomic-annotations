import argparse
import numpy as np
import random
import time

from match_cell_type_regulation import annotate_cell_types_vector, annotate_cell_types_matrix


def generate_random_items(n):
    items = []
    for _ in range(n):
        item = [
            random.randint(1, 24),
            random.randint(1, 100000),
            random.randint(1, 90)
        ]
        items.append(item)
    return items

def main():
    # Create the argument parser
    parser = argparse.ArgumentParser()
    parser.add_argument("db_file", help="Path to DB file")
    parser.add_argument("genome_reference", choices=[37, 38], type=int, help="37 for gh37 or 38 for gh38")
    parser.add_argument("number_of_samples_to_annotate", type=int, nargs='?', default=1, help="The number of samples you want to annotate")
    parser.add_argument("output_format", choices=['flat', 'matrix'], nargs='?', default='flat', help="Specify the output format: flat or matrix")
    parser.add_argument("-s", "--input_sample", nargs='+', type=int,
                        help="Specify the input sample (only for 'sample' operation)")

    # Parse the command-line arguments
    args = parser.parse_args()

    # Access the parsed arguments
    db_file = args.db_file
    genome_reference = args.genome_reference
    number_of_samples_to_annotate = args.number_of_samples_to_annotate
    output_format = args.output_format
    input_sample = args.input_sample

    db = np.load(db_file, allow_pickle=True)
    nucleotides = db["arr_0"]
    values = db["arr_1"]

    if input_sample is not None:
        return annotate_cell_types_vector(nucleotides, values, genome_reference, input_sample[0], input_sample[1], input_sample[2])\
            if output_format == 'flat' else\
            annotate_cell_types_matrix(nucleotides, values, genome_reference, input_sample[0], input_sample[1], input_sample[2])

    start_time_init = time.time()
    # Generate random items
    random_items = generate_random_items(number_of_samples_to_annotate)
    end_time_init = time.time()
    temp = []
    if output_format == 'flat':
        start_time = time.time()
        for item in random_items:
            temp.append(
                annotate_cell_types_vector(nucleotides, values, genome_reference, chromosome=item[0], start_position=item[1], end_position=item[1] + 200))
        end_time = time.time()
    else:
        start_time = time.time()
        for item in random_items:
            temp.append(
                annotate_cell_types_matrix(nucleotides, values, genome_reference, chromosome=item[0], start_position=item[1], end_position=item[1] + 200))
        end_time = time.time()

    print("Elapsed init time: {:.4f} seconds".format(end_time_init - start_time_init))
    print("Elapsed time: {:.4f} seconds".format(end_time - start_time))


if __name__ == '__main__':
    main()
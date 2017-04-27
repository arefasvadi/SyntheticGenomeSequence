from ReferenceSequence import ReferenceSequence
from SyntheticSequence import SyntheticSequence
from datetime import datetime


def main():
    start_time = datetime.now()
    SYNTHETIC_DATASET_SIZE = 500
    print "\nStart time is: " + str(start_time.time()) + "\n"
    reference_file_path = "references/reference3500/"
    reference_file_name = "reference3500.txt"
    variation_file_path = "variation-types/"
    variation_file_name = "variations.txt"
    output_file_path = "references/reference3500/500-sample1/"

    # print "Reading the reference file"
    ref_seq = ReferenceSequence(reference_file_path + reference_file_name, variation_file_path + variation_file_name)
    ref_seq.accept_synthetic_sequences(create_empty_synthetic_sequence(SYNTHETIC_DATASET_SIZE, ref_seq))
    ref_seq.generate_synthetic_sequence(output_file_path)
    # print "Reading the reference file was successful!"
    # print "Processing the variation file"
    finish_time = datetime.now()
    elapsed_time = finish_time - start_time

    print "\nFinish time is: " + str(finish_time.time())
    print "It took " + str(elapsed_time)


def create_empty_synthetic_sequence(size, ref_seq):
    print "*********Started building empty synthetic sequences*********"
    synthetic_sequences = list()
    for i in range(1, size + 1):
        synthetic_sequences.append(SyntheticSequence(i, ref_seq))

    print "Synthetic list's size is: " + str(len(synthetic_sequences))
    print "*********Finished building empty synthetic sequences*********\n"
    return synthetic_sequences


if __name__ == '__main__':
    main()

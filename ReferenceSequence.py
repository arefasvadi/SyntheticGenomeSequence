from Utility import Utility
from VariationController import VariationController, DeletionController, InsertionController, SubstitutionController
from SyntheticSequence import SyntheticSequence
from sortedcontainers import SortedSet, SortedList, SortedDict
from random import Random
from datetime import datetime


class ReferenceSequence(object):
    def __init__(self, reference_file_path_fname, variation_file_path_name):

        print "*********Started reading sequence file*********"
        self.sequence = Utility.readfile(reference_file_path_fname).split("\n")
        self.sequence_id = self.sequence[0][:-1]
        self.sequence = self.sequence[1].upper()
        self.sequence_size = len(self.sequence)
        # print "Reference sequence is: " + self.sequence
        print "Reference id is: " + str(self.sequence_id)
        print "Reference length is: " + str(self.sequence_size)
        print "*********Finished reading sequence file*********\n"

        self.variation_controller = VariationController(variation_file_path_name, self)
        self.__choose_all_hot_spots()
        self.__choose_high_variance_hot_spots()
        self.__assign_variation_type_to_hot_spots()

    def generate_synthetic_sequence(self):

        print "*********Started generating synthetic sequences*********"
        # cnt = 0
        for synthetic_seq in self.synthetic_sequences:
            # if(cnt > 2):
            #     break
            # Choose random hot spots for this sequence
            synthetic_seq.hot_locuses = self.__synthetic_sequence_hot_spots()
            # Identify the type of change
            start_pos = 0
            for locus in synthetic_seq.hot_locuses:
                synthetic_seq.fill_before_locus(start_pos, locus - 1)
                if (self.hot_spots_variation_dict[locus] == "sub"):
                    high_variance = False
                    if (self.high_variance_hot_spots.__contains__(locus)):
                        high_variance = True

                    start_pos = locus + \
                                self.__substitution_handler(synthetic_seq, high_variance, locus)
                elif (self.hot_spots_variation_dict[locus] == "ins"):
                    self.__insertion_handler(synthetic_seq)
                elif (self.hot_spots_variation_dict[locus] == "del"):
                    self.__deletion_handler(synthetic_seq)
            if (start_pos < self.sequence_size):
                synthetic_seq.fill_before_locus(start_pos, synthetic_seq.get_sequence_size() - 1)
                    # cnt = cnt + 1
            print synthetic_seq.print_vcf()
        print "*********Finished generating synthetic sequences*********\n"

    def __substitution_handler(self, synthetic_seq, high_variance, locus):
        variations_list, variations_set = self.variation_controller.sub_controller.get_locuses(locus)
        sum = 0.0
        chosen_size = 0
        distribution_table_list = list()
        for key in self.variation_controller.sub_controller.possible_sizes_prob.keys():
            distribution_table_list.append(
                [sum, sum + self.variation_controller.sub_controller.possible_sizes_prob[key], key])
            sum = sum + self.variation_controller.sub_controller.possible_sizes_prob[key]

        if ((variations_set is None) or len(variations_set) < self.variation_controller.MAX_VARIATION_CAP):
            rvalue = self.variation_controller.sub_controller.r.uniform(0.0, 1.0)
            for i in range(len(distribution_table_list)):
                if (rvalue >= distribution_table_list[i][0] and rvalue <= distribution_table_list[i][1]):
                    chosen_size = distribution_table_list[i][2]
                    break
            if (chosen_size == 0):
                raise ValueError("The chosen size is zero!")

            ref_string = self.sequence[locus:locus + chosen_size - 1]
            seq_string = self.variation_controller.sub_controller.substitution_mutation_generator(ref_string)
            # self.variation_controller.sub_controller.add_value_to_locus(locus, ref_string)
            synthetic_seq.sequence = synthetic_seq.sequence + seq_string
            self.variation_controller.sub_controller.add_value_to_locus(locus, seq_string)
            for i in range(len(seq_string)):
                synthetic_seq.vcf.add([locus + i, "sub", seq_string[i]])

        else:
            most_common_value, max_occurrence, total_length = self.variation_controller.sub_controller.most_common_occurrence(
                locus)
            lower_threshold = 0.0
            upper_threshold = 1.0
            if (high_variance):
                lower_threshold = self.variation_controller.min_percent_high_variance_share
                upper_threshold = self.variation_controller.max_percent_high_variance_share
            elif (not high_variance):
                lower_threshold = self.variation_controller.min_percent_low_variance_share
                upper_threshold = self.variation_controller.max_percent_low_variance_share

            current_occurrence = float(max_occurrence / total_length)
            if (current_occurrence > lower_threshold and current_occurrence < upper_threshold):
                rvalue = 0
                while True:
                    rvalue = self.variation_controller.sub_controller.r.randint(0, len(variations_set) - 1)
                    if (variations_set[rvalue] != most_common_value):
                        break
                seq_string = variations_set[rvalue]
                synthetic_seq.sequence = synthetic_seq.sequence + seq_string
                self.variation_controller.sub_controller.add_value_to_locus(locus, seq_string)
                chosen_size = len(seq_string)
                for i in range(len(seq_string)):
                    synthetic_seq.vcf.add([locus + i, "sub", seq_string[i]])
            else:
                synthetic_seq.sequence = synthetic_seq.sequence + most_common_value
                self.variation_controller.sub_controller.add_value_to_locus(locus, most_common_value)
                for i in range(len(most_common_value)):
                    synthetic_seq.vcf.add([locus + i, "sub", most_common_value[i]])
                chosen_size = len(most_common_value)
        return chosen_size

    def __synthetic_sequence_hot_spots(self):
        synthetic_seq_hot_spots = SortedSet()
        r = Random()
        r.seed(datetime.now())
        while (len(synthetic_seq_hot_spots) <
                   int(self.variation_controller.expected_pairwise_variation_to_reference * self.sequence_size)):
            rvalue = r.randint(0, len(self.all_hot_spots) - 1)
            synthetic_seq_hot_spots.add(self.all_hot_spots[rvalue])
        return synthetic_seq_hot_spots

    def accept_synthetic_sequences(self, synthetic_seqs):
        self.synthetic_sequences = synthetic_seqs

    def __assign_variation_type_to_hot_spots(self):

        print "*********Started assigning variation types to each hot spot*********"
        self.hot_spots_variation_dict = SortedDict()
        rand = Random()
        rand.seed(datetime.now())

        for hot_spot in self.all_hot_spots:
            rvalue = rand.uniform(0.0, 1.0)
            if (rvalue >= 0.0 and rvalue < SubstitutionController.OCCURRENCE_RATIO):
                self.hot_spots_variation_dict[hot_spot] = "sub"
            elif (rvalue >= SubstitutionController.OCCURRENCE_RATIO and
                          rvalue < SubstitutionController.OCCURRENCE_RATIO + InsertionController.OCCURRENCE_RATIO):
                self.hot_spots_variation_dict[hot_spot] = "ins"
            else:
                self.hot_spots_variation_dict[hot_spot] = "del"
        print self.hot_spots_variation_dict
        print "*********Finished assigning variation types to each hot spot*********\n"

    def __choose_high_variance_hot_spots(self):

        print "*********Started picking high variance hot locations*********"
        self.high_variance_hot_spots = SortedSet()
        r = Random()
        r.seed(datetime.now())
        while (len(self.high_variance_hot_spots) < int(
                    self.variation_controller.high_variance_hot_spot_percentage * self.sequence_size)):
            rvalue = r.randint(0, len(self.all_hot_spots) - 1)
            self.high_variance_hot_spots.add(self.all_hot_spots[rvalue])
        print "Size of high variance hot spot locations is: " + str(len(self.high_variance_hot_spots))
        print "*********Finished picking high variance hot locations*********\n"

    def __choose_all_hot_spots(self):

        print "*********Started generating random locations*********"
        self.all_hot_spots = SortedSet()
        r = Random()
        r.seed(datetime.now())
        while (len(self.all_hot_spots) < int(self.variation_controller.all_hot_spot_percentage * self.sequence_size)):
            rvalue = r.randint(0, self.sequence_size - 1)
            if (self.all_hot_spots.__contains__(rvalue)):
                continue
            elif (not self.__is_hot_spot_valid_in_distance(rvalue)):
                continue
            self.all_hot_spots.add(rvalue)
            if (len(self.all_hot_spots) % 100 == 0):
                print len(self.all_hot_spots)

        print "Size of all random location hot spots is: " + str(len(self.all_hot_spots))
        # print self.all_hot_spots
        print "*********Finished generating random locations*********\n"

    def __is_hot_spot_valid_in_distance(self, rvalue):

        valid = False
        self.all_hot_spots.add(rvalue)
        index = self.all_hot_spots.index(rvalue)
        if (len(self.all_hot_spots) == 1):
            valid = True
        else:
            if (index - 1 >= 0):
                if (abs(self.all_hot_spots[index - 1] - rvalue) <= int(
                            self.variation_controller.hot_spot_max_interval_to_size_ratio * self.sequence_size)
                    and abs(self.all_hot_spots[index - 1] - rvalue) >= int(
                            self.variation_controller.hot_spot_min_interval_to_size_ratio * self.sequence_size)):
                    valid = True
            if (index + 1 <= len(self.all_hot_spots) - 1):
                if (valid == True):
                    if (abs(self.all_hot_spots[index + 1] - rvalue) <= int(
                                self.variation_controller.hot_spot_max_interval_to_size_ratio * self.sequence_size)
                        and abs(self.all_hot_spots[index + 1] - rvalue) >= int(
                                self.variation_controller.hot_spot_min_interval_to_size_ratio * self.sequence_size)):
                        valid = True
                    else:
                        valid = False
                else:
                    if (abs(self.all_hot_spots[index + 1] - rvalue) <= int(
                                self.variation_controller.hot_spot_max_interval_to_size_ratio * self.sequence_size)
                        and abs(self.all_hot_spots[index + 1] - rvalue) >= int(
                                self.variation_controller.hot_spot_min_interval_to_size_ratio * self.sequence_size)):
                        valid = True

        self.all_hot_spots.remove(rvalue)
        return valid

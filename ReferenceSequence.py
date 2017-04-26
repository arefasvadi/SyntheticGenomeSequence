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
            for locus in synthetic_seq.hot_locuses:
                if (self.hot_spots_variation_dict[locus] == "sub"):
                    high_variance = False
                    if (self.high_variance_hot_spots.__contains__(locus)):
                        high_variance = True
                    self.__substitution_handler(synthetic_seq, high_variance, locus)
                elif (self.hot_spots_variation_dict[locus] == "ins"):
                    self.__insertion_handler(synthetic_seq)
                elif (self.hot_spots_variation_dict[locus] == "del"):
                    self.__deletion_handler(synthetic_seq)

                    # cnt = cnt + 1
        print "*********Finished generating synthetic sequences*********\n"

    def __substitution_handler(self, synthetic_seq, high_variance, locus):
        self.variation_controller.sub_controller.get_locuses()[locus]

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
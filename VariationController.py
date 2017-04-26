from Utility import Utility
from sortedcontainers import SortedDict, SortedList, SortedSet
from random import Random
from datetime import datetime

class VariationController(object):
    MAX_VARIATION_CAP = 10

    def __init__(self, pathandname, ref_seq_instance):

        configs = Utility.readfile(pathandname)
        self.sub_controller = SubstitutionController()
        self.ins_controller = InsertionController()
        self.del_controller = DeletionController()
        self.__parse_configs(configs)
        self.sub_controller.read_variation_type_specific_configs("variation-types/substitution.txt")
        self.ins_controller.read_variation_type_specific_configs("variation-types/insertion.txt")
        self.del_controller.read_variation_type_specific_configs("variation-types/deletion.txt")

    def __parse_configs(self, configs):
        print "*********Started reading configs*********"
        configs = configs.split("\n")
        for line in configs:
            if (str(line).startswith("#") or str(line) == ""):
                continue

            line_portions = line.split("=>")
            for portion in line_portions:
                portion = portion.strip()
            # if ''.join(sorted(line_portions[0])).strip() == ''.join(sorted("reference-all-hot-spots-locus")).strip():
            if (''.join(sorted(line_portions[0])).strip() == ''.join(sorted("reference-all-hot-spots-locus")).strip()):
                print line_portions[0] + " => " + line_portions[1]
                self.all_hot_spot_percentage = float(line_portions[1])

            elif (''.join(sorted(line_portions[0])).strip() == ''.join(
                    sorted("reference-hot-spots-variance-locus")).strip()):
                print line_portions[0] + " => " + line_portions[1]
                self.high_variance_hot_spot_percentage = float(line_portions[1])

            elif (''.join(sorted(line_portions[0])).strip() == ''.join(
                    sorted("min-percent-high-variance-share")).strip()):
                print line_portions[0] + " => " + line_portions[1]
                self.min_percent_high_variance_share = float(line_portions[1])

            elif (''.join(sorted(line_portions[0])).strip() == ''.join(
                    sorted("max-percent-high-variance-share")).strip()):
                print line_portions[0] + " => " + line_portions[1]
                self.max_percent_high_variance_share = float(line_portions[1])

            elif (
                        ''.join(sorted(line_portions[0])).strip() == ''.join(
                        sorted("min-percent-low-variance-share")).strip()):
                print line_portions[0] + " => " + line_portions[1]
                self.min_percent_low_variance_share = float(line_portions[1])

            elif (
                        ''.join(sorted(line_portions[0])).strip() == ''.join(
                        sorted("max-percent-low-variance-share")).strip()):
                print line_portions[0] + " => " + line_portions[1]
                self.max_percent_low_variance_share = float(line_portions[1])

            elif (''.join(sorted(line_portions[0])).strip() == ''.join(
                    sorted("max-cap-hot-spot-type-variation")).strip()):
                print line_portions[0] + " => " + line_portions[1]
                VariationController.MAX_VARIATION_CAP = float(line_portions[1])

            elif (''.join(sorted(line_portions[0])).strip() == ''.join(
                    sorted("expected-pairwise-variation-to-reference")).strip()):
                print line_portions[0] + " => " + line_portions[1]
                self.expected_pairwise_variation_to_reference = float(line_portions[1])

            elif (''.join(sorted(line_portions[0])).strip() == ''.join(
                    sorted("hot-spot-min-interval-to-size-ratio")).strip()):
                print line_portions[0] + " => " + line_portions[1]
                self.hot_spot_min_interval_to_size_ratio = float(line_portions[1])

            elif (''.join(sorted(line_portions[0])).strip() == ''.join(
                    sorted("hot-spot-max-interval-to-size-ratio")).strip()):
                print line_portions[0] + " => " + line_portions[1]
                self.hot_spot_max_interval_to_size_ratio = float(line_portions[1])

            elif (''.join(sorted(line_portions[0])).strip() == ''.join(
                        sorted("substitution_occurrence_ratio")).strip()):
                print line_portions[0] + " => " + line_portions[1]
                SubstitutionController.OCCURRENCE_RATIO = float(line_portions[1])

            elif (''.join(sorted(line_portions[0])).strip() == ''.join(sorted("insertion_occurrence_ratio")).strip()):
                print line_portions[0] + " => " + line_portions[1]
                InsertionController.OCCURRENCE_RATIO = float(line_portions[1])
            elif (''.join(sorted(line_portions[0])).strip() == ''.join(sorted("deletion_occurrence_ratio")).strip()):
                print line_portions[0] + " => " + line_portions[1]
                DeletionController.OCCURRENCE_RATIO = float(line_portions[1])
            else:
                print "Unknown Parameter/\\/\\/\\/\\/\\"
        print "*********Finished reading configs*********\n"


class VariationTypeController(object):

    def __init__(self):
        self.locuses = SortedDict()
        self.locuses_set = SortedDict()
        self.r = Random()
        self.r.seed(datetime.now())

    def get_locuses(self):
        return self.locuses, self.locuses_set

    def add_value_to_locus(self, key, value):
        if (not self.locuses.has_key(key)):
            values_list = SortedList(value)
            self.locuses[key] = values_list
        else:
            values_list = self.locuses[key]
            values_list.add(value)
            self.locuses[key] = values_list

        if (not self.locuses_set.has_key(key)):
            values_set = SortedSet(value)
            self.locuses_set[key] = values_set
        else:
            values_set = self.locuses_set[key]
            values_set.add(value)
            self.locuses_set[key] = values_set

    # def set_occurance_probability(self,prob):
    #     self.occur_prob = prob

    def set_possible_sizes(self, min, max):
        self.possibles_sizes = range(min, max + 1, 1)
        self.possible_sizes_prob = SortedDict()
        for key in self.possibles_sizes:
            self.possible_sizes_prob[key] = None

    def set_possible_size_probability(self, key, value):

        if (self.possible_sizes_prob.has_key(key) == False):
            raise RuntimeError("The key associated with size does not exist!")
        else:
            self.possible_sizes_prob[key] = value

    def most_common_occurrence(self, key):
        if (not self.locuses.has_key(key)):
            raise RuntimeError("The key associated with size does not exist!")
        else:
            values_list = SortedList(self.locuses[key])
            max_occurrence = 0;
            most_common_value = ""

            for i in values_list:
                if str(most_common_value) == str(i):
                    continue
                if (values_list.count(i) > max_occurrence):
                    max_occurrence = values_list.count(i)
                    most_common_value = i
        return most_common_value, max_occurrence

    def read_variation_type_specific_configs(self, pathandname):
        configs = Utility.readfile(pathandname)
        configs = configs.split("\n")
        min_is_checked = False
        max_is_checked = False
        min = 1
        max = 1

        for line in configs:
            if (str(line).startswith("#") or str(line) == ""):
                continue
            line_portions = line.split("=>")
            for portion in line_portions:
                portion = portion.strip()
            # if ''.join(sorted(line_portions[0])).strip() == ''.join(sorted("reference-all-hot-spots-locus")).strip():
            if (''.join(sorted(line_portions[0])).strip() == ''.join(sorted("min")).strip()
                and not min_is_checked):
                print line_portions[0] + " => " + line_portions[1]
                min = int(line_portions[1])
                min_is_checked = True
                if (min_is_checked and max_is_checked):
                    self.set_possible_sizes(min, max)

            elif (''.join(sorted(line_portions[0])).strip() == ''.join(sorted("max")).strip()
                  and not max_is_checked):
                print line_portions[0] + " => " + line_portions[1]
                max = int(line_portions[1])
                max_is_checked = True
                if (min_is_checked and max_is_checked):
                    self.set_possible_sizes(min, max)

            elif (max_is_checked and min_is_checked):
                print line_portions[0] + " => " + line_portions[1]
                index = [i for i, x in enumerate(self.possibles_sizes) if x == int(line_portions[0])]
                self.set_possible_size_probability(self.possibles_sizes[index[0]], float(line_portions[1]))
            else:
                print "Unknown Parameter/\\/\\/\\/\\/\\"


class SubstitutionController(VariationTypeController):
    OCCURRENCE_RATIO = 0.8

    def __init__(self):
        super(SubstitutionController, self).__init__()

    def read_variation_type_specific_configs(self, pathandname):
        print "********Started Reading Substitution Configs********"
        super(SubstitutionController, self).read_variation_type_specific_configs(pathandname)
        print "********Finished Reading Substitution Configs********\n"


class InsertionController(VariationTypeController):
    OCCURRENCE_RATIO = 0.1

    def __init__(self):
        super(InsertionController, self).__init__()

    def read_variation_type_specific_configs(self, pathandname):
        print "********Started Reading Insertion Configs********"
        super(InsertionController, self).read_variation_type_specific_configs(pathandname)
        print "********Finished Reading Insertion Configs********\n"


class DeletionController(VariationTypeController):
    OCCURRENCE_RATIO = 0.1

    def __init__(self):
        super(DeletionController, self).__init__()

    def read_variation_type_specific_configs(self, pathandname):
        print "********Started Reading Deletion Configs********"
        super(DeletionController, self).read_variation_type_specific_configs(pathandname)
        print "********Finished Reading Deletion Configs********\n"

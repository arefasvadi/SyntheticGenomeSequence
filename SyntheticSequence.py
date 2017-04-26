from sortedcontainers import SortedList, SortedSet


class SyntheticSequence(object):
    def __init__(self, id, ref_seq):
        self.refrence_seq = ref_seq
        self.sequence_id = id
        self.sequence = str("")
        self.vcf = SortedList()
        # self.sequence_size = 0
        self.hot_locuses = SortedSet()

    def sequence_size(self):
        return len(self.sequence)

    def fill_before_locus(self, start, end):
        self.hot_locuses = self.hot_locuses + self.refrence_seq.sequence[start:end]

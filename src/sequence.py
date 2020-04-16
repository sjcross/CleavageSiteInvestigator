import re

class Sequence():
    def __init__(self, seq):
        self._seq = seq

    def __str__(self):
        return self.get_sequence_string()

    def get_sequence_string(self):
        return self._seq

    def get_complement_sequence(self):
        complement = self._seq
        complement = complement.replace("A", "1")
        complement = complement.replace("a", "2")
        complement = complement.replace("C", "3")
        complement = complement.replace("c", "4")
        complement = complement.replace("T", "A")
        complement = complement.replace("t", "a")
        complement = complement.replace("G", "C")
        complement = complement.replace("g", "c")
        complement = complement.replace("1", "T")
        complement = complement.replace("2", "t")
        complement = complement.replace("3", "G")
        complement = complement.replace("4", "g")
        
        return Sequence(complement)

    def get_reverse_sequence(self):
        return Sequence(self._seq[::-1])

    def get_reverse_complement_sequence(self):
        return self.get_complement_sequence().get_reverse_sequence()

    def get_start_sequence(self, num_bases, offset=0):
        return Sequence(self._seq[offset:(num_bases+offset)])
        
    def get_end_sequence(self, num_bases, offset=0):
        if offset == 0:
            return Sequence(self._seq[(-num_bases-offset)::])
        else:
            return Sequence(self._seq[(-num_bases-offset):-offset])

    def find_sequence_match(self, test_seq):
        pattern = re.compile(test_seq.get_sequence_string(), re.IGNORECASE)
        indices = [m.start(0) for m in pattern.finditer(self._seq)]

        return len(indices)

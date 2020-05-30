from Bio import Align
from ends import Ends


class SequenceSearcher():
    def __init__(self, aligner, min_quality=1.0, num_bases=20, verbose=False):
        self._aligner = aligner
        self._min_quality = min_quality
        self._num_bases = num_bases
        self._verbose = verbose

    def get_aligner(self):
        return self._aligner

    def set_aligner(self, aligner):
        self._aligner = aligner

    def get_min_quality(self):
        return self._min_quality

    def set_min_quality(self, min_quality):
        self._min_quality = min_quality

    def get_num_bases(self):
        return self._num_bases

    def set_num_bases(self, num_bases):
        self._num_bases = num_bases

    def get_verbose(self):
        return self._verbose

    def set_verbose(self, verbose):
        self._verbose = verbose

    def process_sanger(self, ref, cass, test1, test2):
        # Finding break in test sequence 1
        if self._verbose:
            print("    Finding break in test sequence 1")

        (pos1, score) = self._find_cassette_end(cass, test1)
        (clevage_site1, isRC1) = self._find_clevage_site(ref, test1, pos1)
        if clevage_site1 is None:
            return(None,None)

        # Finding break in test sequence 2
        if self._verbose:
            print("    Finding break in test sequence 2")

        (pos2, score) = self._find_cassette_end(cass, test2)
        (clevage_site2, isRC2) = self._find_clevage_site(ref, test2, pos2)
        if clevage_site2 is None:
            return(None,None)

        # Only one should be RC
        if isRC1 and isRC2:
            print("    ERROR: Both identified as reverse complement")
            return(None,None)
        elif not isRC1 and not isRC2:
            print("    ERROR: Neither identified as reverse complement")
            return(None,None)

        # If clevage_site1 is RC, switch clevage sites
        if isRC2:
            (clevage_site1, clevage_site2) = (clevage_site2, clevage_site1)

        if self._verbose:
            print("\r")

        return (clevage_site1, clevage_site2)

    def process_other(self, ref, cass, test):
        # Finding start of test sequence
        if self._verbose:
            print("    Finding start of test sequence")

        alignment1 = self._search_sequence(cass, test, Ends.CASS_START)
        if alignment1.score == 0:
            if self._verbose:
                print("    Inverting cassette")

            cass = cass.reverse_complement()      
            alignment1 = self._search_sequence(cass, test, Ends.CASS_START)
        
        pos1 = alignment1.path[0][0]
        (alignment1, isRC1) = self._find_target_in_ref(ref, test, pos1-self._num_bases, self._num_bases)
        clevage_site1 = alignment1.path[1][0]

        # Finding end of test sequence
        if self._verbose:
            print("    Finding end of test sequence")

        pos2 = self._search_sequence(cass, test, Ends.CASS_END).path[1][0]
        (alignment2, isRC2) = self._find_target_in_ref(ref, test, pos2, self._num_bases)
        clevage_site2 = alignment2.path[0][0]    

        if self._verbose:
            print("\r")

        return (clevage_site1, clevage_site2)

    def _find_cassette_end(self, cass, test):
        # Finding cassette ends in test sequence
        if self._verbose:
            print("        Finding cassette end in test sequence")
        
        max_alignment = Align.PairwiseAlignment(
            target="", query="", path=((0, 0), (0, 0)), score=0.0)
        end = 0

        alignment = self._search_sequence(cass, test, Ends.CASS_START_RC)
        if alignment.score > 0:
            if alignment.score > max_alignment.score:
                max_alignment = alignment
                end = Ends.CASS_START_RC

        # Testing against cassette end
        alignment = self._search_sequence(cass, test, Ends.CASS_END)
        if alignment.score > max_alignment.score:
            max_alignment = alignment
            end = Ends.CASS_END

        if self._verbose:
            if end == Ends.CASS_START_RC:
                pos_string = "start RC"
            elif end == Ends.CASS_END:
                pos_string = "end"

            print("            Best match for cassette %s (%s)" %
                (pos_string, max_alignment.query))
            print("            Match score = %0.2f (quality %0.2f)" %
                (max_alignment.score, get_quality(max_alignment)))

        return (max_alignment.path[-1][0], max_alignment.score)

    def _find_clevage_site(self, ref, test, pos):
        # Getting region of test sequence to match to reference
        if self._verbose:
            print("        Finding cassette-adjacent sequence in reference")
            
        (alignment, isRC) = self._find_target_in_ref(ref, test, pos, self._num_bases)
        
        if alignment is None:
            return(alignment, False)

        if self._verbose:
            print("            Match score = %0.2f (quality %0.2f)" %
                (alignment.score, get_quality(alignment)))
        if isRC:
            clevage_position = alignment.path[-1][0]
        else:
            clevage_position = alignment.path[0][0]

        ref_break_seq_left = ref[clevage_position - 5:clevage_position]
        ref_break_seq_right = ref[clevage_position: clevage_position + 5]

        if self._verbose:
            print("            Reference break at position %i" % clevage_position)
        
        if self._verbose:
            print("            Reference break at sequence %s | %s" %
                (ref_break_seq_left, ref_break_seq_right))

        return (clevage_position, isRC)
    
    def _search_sequence(self, cass, test, end):
        max_alignment = Align.PairwiseAlignment(
            target="", query="", path=((0, 0), (0, 0)), score=0.0)

        if end is Ends.CASS_START:
            # Testing against cassette start
            cass_start = cass[0:self._num_bases]
            alignments = self._aligner.align(test, cass_start)
            for alignment in alignments:
                if alignment.score > max_alignment.score and get_quality(alignment) >= self._min_quality:
                    max_alignment = alignment

        elif end is Ends.CASS_START_RC:
            # Testing against cassette start (RC)
            cass_start_rc = cass[0:self._num_bases].reverse_complement()
            alignments = self._aligner.align(test, cass_start_rc)
            for alignment in alignments:
                if alignment.score > max_alignment.score and get_quality(alignment) >= self._min_quality:
                    max_alignment = alignment
        
        elif end is Ends.CASS_END:
            # Testing against cassette end
            cass_end = cass[-self._num_bases::]
            alignments = self._aligner.align(test, cass_end)
            for alignment in alignments:
                if alignment.score > max_alignment.score and get_quality(alignment) >= self._min_quality:
                    max_alignment = alignment

        elif end is Ends.CASS_END_RC:
            # Testing against cassette end
            cass_end_rc = cass[-self._num_bases::].reverse_complement()
            alignments = self._aligner.align(test, cass_end_rc)
            for alignment in alignments:
                if alignment.score > max_alignment.score and get_quality(alignment) >= self._min_quality:
                    max_alignment = alignment

        return max_alignment

    def _find_target_in_ref(self, ref, test, pos, search_length):
        (alignments, isRC) = self._find_all_targets_in_ref(ref, test, pos, search_length)

        # If more than one match was found, increasing the number of bases used in the recognition
        while len(alignments) > 1:
            search_length = search_length + 1
            (alignments, isRC) = self._find_all_targets_in_ref(ref, test, pos, search_length)

        if len(alignments) == 0:
            print("\nNo match found\n")
            return (None, False)

        return (alignments[0], isRC)

    def _find_all_targets_in_ref(self, ref, test, pos, search_length):
        test_target = test[pos: pos + search_length]

        # Finding test target in reference sequence
        alignments = self._aligner.align(ref, test_target)
        valid_alignments = self._get_valid_alignments(alignments)

        # If no matches were found, try the reverse complement of the test target
        isRC = False
        if len(valid_alignments) == 0:
            if self._verbose:
                print("        Reversing test target sequence")

            isRC = True
            test_target = test_target.reverse_complement()
            alignments = self._aligner.align(ref, test_target)
            valid_alignments = self._get_valid_alignments(alignments)

        if self._verbose:
            print("        Testing reference (%s) with %i bases (%i matches found)" %
                  (test_target, search_length, len(valid_alignments)))

        return (valid_alignments, isRC)

    def _get_valid_alignments(self, alignments):
        valid_alignments = []

        for alignment in alignments:
            quality = get_quality(alignment)
            if quality >= self._min_quality:
                valid_alignments.append(alignment)

        return valid_alignments        

def get_quality(alignment):
    return alignment.score / len(alignment.query)

def print_clevage_site(ref, clevage_site1, clevage_site2):
    if clevage_site1 is None or clevage_site2 is None:
        return

    # Identifying break type
    if clevage_site1 < clevage_site2:
        # 3' overhang
        left_seq1 = ref[clevage_site1 - 1 :clevage_site1]
        mid_seq1 = ref[clevage_site1:clevage_site2]
        right_seq1 = ref[clevage_site2: clevage_site2 + 1]

        left_seq2 = left_seq1.complement()
        mid_seq2 = mid_seq1.complement()
        right_seq2 = right_seq1.complement()

        print("    Restriction site (3' overhang):")
        print("        5'...%s %s↓%s...3'\r\n        3'...%s↑%s %s...5'\r\n" % (left_seq1, mid_seq1, right_seq1, left_seq2, mid_seq2, right_seq2))
        
    elif (clevage_site1 == clevage_site2):
        # Blunt end
        left_seq1 = ref[clevage_site1 - 3 :clevage_site1]
        right_seq1 = ref[clevage_site2: clevage_site2 + 3]

        left_seq2 = left_seq1.complement()
        right_seq2 = right_seq1.complement()

        print("    Restriction site (blunt end):")
        print("        5'...%s↓%s...3'\r\n        3'...%s↑%s...5'\r\n" % (left_seq1, right_seq1, left_seq2, right_seq2))

    elif clevage_site1 > clevage_site2:
        # 5' overhang
        left_seq1 = ref[clevage_site2 - 1 :clevage_site2]
        mid_seq1 = ref[clevage_site2:clevage_site1]
        right_seq1 = ref[clevage_site1: clevage_site1 + 1]

        left_seq2 = left_seq1.complement()
        mid_seq2 = mid_seq1.complement()
        right_seq2 = right_seq1.complement()

        print("    Restriction site (5' overhang):")
        print("        5'...%s↓%s %s...3'\r\n        3'...%s %s↑%s...5'\r\n" % (left_seq1, mid_seq1, right_seq1, left_seq2, mid_seq2, right_seq2))


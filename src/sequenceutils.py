from Bio import Align
from ends import Ends


class SequenceSearcher():
    def __init__(self, aligner, max_gap=10, min_quality=1.0, num_bases=20, verbose=False):
        self._aligner = aligner
        self._max_gap = max_gap
        self._min_quality = min_quality
        self._num_bases = num_bases
        self._verbose = verbose

    def get_aligner(self):
        return self._aligner

    def set_aligner(self, aligner):
        self._aligner = aligner

    def get_max_gap(self):
        return self._max_gap

    def set_max_gap(self, max_gap):
        self._max_gap = max_gap

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

    def process_other(self, ref, cass, test):
        if self._verbose:
            print("        Finding first cassette end in test sequence")
        (cass_pos_1, cass1_isRC) = self._find_best_cassette_end(cass, test, Ends.CASS_START)
        
        if self._verbose:
            print("        Finding second cassette end in test sequence")
        (cass_pos_2, cass2_isRC) = self._find_best_cassette_end(cass, test, Ends.CASS_END)

        # Both should be RC or normal
        if cass1_isRC != cass2_isRC:
            if self._verbose:
                print("ERROR: Cassette RC mismatch\n")
            return (None, None)

        # Checking results are OK
        if cass_pos_1 is None or cass_pos_2 is None:
            if self._verbose:
                print("ERROR: Cassette ends not found\n")
            return (None, None)

        # Finding cassette-adjacent test sequence in reference
        if self._verbose:
            print("        Finding first cassette-adjacent test sequence in reference sequence")
        (alignment1, isRC1) = self._find_best_target_in_ref(ref, test, cass_pos_1.path, self._num_bases, self._num_bases)
        
        if self._verbose:
            print("        Finding second cassette-adjacent test sequence in reference sequence")
        (alignment2, isRC2) = self._find_best_target_in_ref(ref, test, cass_pos_2.path, self._num_bases, 0)

        if alignment1 is None or alignment2 is None:
            if self._verbose:
                print("ERROR: Test sequence not found in reference\n")
            return (None, None)

        # Both should be RC or normal
        if isRC1 != isRC2:
            if self._verbose:
                print("ERROR: RC mismatch\n")
            return (None, None)

        if isRC1:
            clevage_site_1 = alignment2.path[-1][0]
            clevage_site_2 = alignment1.path[0][0]         
        else:
            clevage_site_1 = alignment1.path[-1][0]
            clevage_site_2 = alignment2.path[0][0]

        # primer_truncation_length(cass, ref, clevage_site_1, clevage_site_2)

        if abs(clevage_site_1 - clevage_site_2) > self._max_gap:
            if self._verbose:
                print("ERROR: Clevage site gap (%i) exceeds maximum permitted\n" % abs(clevage_site_1 - clevage_site_2))
            return (None, None)
            
        return (clevage_site_1, clevage_site_2)
    
    def _find_best_cassette_end(self, cass, test, end):
        cass_pos = self._find_cassette_end(cass, test, end)     
        cass_pos_rc = self._find_cassette_end(cass.reverse_complement(), test, end)

        if cass_pos is None and cass_pos_rc is not None:
            if self._verbose:
                print("                Best score = %f" % cass_pos_rc.score)
            return (cass_pos_rc, True)

        elif cass_pos is not None and cass_pos_rc is None:
            if self._verbose:
                print("                Best score = %f" % cass_pos.score)
            return (cass_pos, False)

        elif cass_pos is None and cass_pos_rc is None:
            return (None, False)

        if cass_pos_rc.score > cass_pos.score:
            if self._verbose:
                print("                Best score = %f" % cass_pos_rc.score)
            return (cass_pos_rc, True)
        else:
            if self._verbose:
                print("                Best score = %f" % cass_pos.score)
            return (cass_pos, False)

        cass1_isRC = cass_pos_rc.score > cass_pos.score
        if cass1_isRC:
            cass_pos = cass_pos_rc
        
    def _find_cassette_end(self, cass, test, end):
        max_alignment = Align.PairwiseAlignment(
            target="", query="", path=((0, 0), (0, 0)), score=0.0)

        if end is Ends.CASS_START:
            # Checking for "full" cassette end
            alignments = self._aligner.align(test, cass[0: self._num_bases])
            max_alignment = self._get_max_alignment(alignments)

        elif end is Ends.CASS_END:
            # Checking for "full" cassette end
            alignments = self._aligner.align(test, cass[-self._num_bases ::])
            max_alignment = self._get_max_alignment(alignments)

        if max_alignment is None:
            return None

        if self._verbose:
            print("            Best score = %f" % max_alignment.score)
                  
        return max_alignment

    def _find_best_target_in_ref(self, ref, test, path, search_length, search_offset):
        max_alignment = Align.PairwiseAlignment(
            target="", query="", path=((0, 0), (0, 0)), score=0.0)

        max_isRC = False
        for en in range(len(path)):
            (alignment, isRC) = self._find_target_in_ref(ref, test, path[en][0]-search_offset, self._num_bases)

            if alignment is None:
                continue

            if alignment.score > max_alignment.score:
                max_alignment = alignment
                max_isRC = isRC

        if max_alignment.score == 0 or get_quality(max_alignment) < self._min_quality:
            return (None, False)

        return (max_alignment, max_isRC)

    def _find_target_in_ref(self, ref, test, pos, search_length):
        test_target = test[pos: pos + search_length]

        # Finding test target in reference sequence
        alignments = self._aligner.align(ref, test_target)
        max_alignment = self._get_max_alignment(alignments)

        # If no matches were found, try the reverse complement of the test target
        alignments_rc = self._aligner.align(ref, test_target.reverse_complement())
        max_alignment_rc = self._get_max_alignment(alignments_rc)

        if max_alignment is None and max_alignment_rc is not None:
            if self._verbose:
                print("                Best score = %f" % max_alignment_rc.score)
            return (max_alignment_rc, True)

        elif max_alignment is not None and max_alignment_rc is None:
            if self._verbose:
                print("                Best score = %f" % max_alignment.score)
            return (max_alignment, False)

        elif max_alignment is None and max_alignment_rc is None:
            return (None, False)

        if max_alignment_rc.score > max_alignment.score:
            if self._verbose:
                print("                Best score = %f" % max_alignment_rc.score)
            return (max_alignment_rc, True)
        else:
            if self._verbose:
                print("                Best score = %f" % max_alignment.score)
            return (max_alignment, False)

    def _get_max_alignment(self, alignments):
        max_alignment = Align.PairwiseAlignment(
            target="", query="", path=((0, 0), (0, 0)), score=0.0)

        for alignment in alignments:
            if alignment.score > max_alignment.score and get_quality(alignment) >= self._min_quality:
                max_alignment = alignment

        if max_alignment.score == 0:
            return None

        return max_alignment

def get_quality(alignment):
    return alignment.score / len(alignment.query)

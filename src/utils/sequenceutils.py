from Bio import Align
from enums.ends import Ends

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

    # Returns:
    #     cleavage_site_b = index of nucleotide immediately 5' of cleavage site on top strand
    #     cleavage_site_t = index of nucleotide immediately 3' of cleavage site on bottom strand
    #     local_seq_1 = local sequence at top strand cleavage site (local_r nucleotides either side of cleavage site)
    #     local_seq_2 = local sequence at bottom strand cleavage site (local_r nucleotides either side of cleavage site)
    def get_cleavage_positions(self, ref, cass, test):
        if self._verbose:
            print("        Finding first cassette end in test sequence:")
        (cass_pos_1, cass1_isRC) = self._find_best_cassette_end(cass, test, Ends.CASS_START)
                
        if self._verbose:
            print("        Finding second cassette end in test sequence:")
        (cass_pos_2, cass2_isRC) = self._find_best_cassette_end(cass, test, Ends.CASS_END)
        
        # Both should be RC or normal.
        if cass1_isRC != cass2_isRC:
            if self._verbose:
                print("ERROR: Cassette RC mismatch\n")
            return (None, None)

        # Checking results are OK
        if cass_pos_1 is None or cass_pos_2 is None:
            if self._verbose:
                print("ERROR: Cassette ends not found\n")
            return (None, None)  

        # Note: cass_pos_1 should always be less than cass_pos_2 due to the way it's searched for - it doesn't matter if it's sense or RC

        # Finding cassette-adjacent test sequence in reference
        if self._verbose:
            print("        Finding first cassette-adjacent test sequence in reference sequence:")
        (alignment1, isRC1, test_cut_1) = self._find_best_target_in_ref(ref, test, cass_pos_1.path, self._num_bases, self._num_bases)
        test_cut_1 = test_cut_1 + self._num_bases
        
        if self._verbose:
            print("        Finding second cassette-adjacent test sequence in reference sequence:")
        (alignment2, isRC2, test_cut_2) = self._find_best_target_in_ref(ref, test, cass_pos_2.path, self._num_bases, 0)
        
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
            cleavage_site_t = alignment1.path[0][0]
            cleavage_site_b = alignment2.path[-1][0]
        else:
            cleavage_site_t = alignment2.path[0][0]
            cleavage_site_b = alignment1.path[-1][0]

        if abs(cleavage_site_t - cleavage_site_b) > self._max_gap:
            if self._verbose:
                print("ERROR: cleavage site gap (%i) exceeds maximum permitted\n" % abs(cleavage_site_b - cleavage_site_t))
            return (None, None)

        # (5' overhangs only) both sequences should match
        if cleavage_site_t < cleavage_site_b:
            diff = cleavage_site_b - cleavage_site_t
            overhang_1 = test[test_cut_1 - diff :test_cut_1]
            overhang_2 = test[test_cut_2: test_cut_2 + diff]
            
            if overhang_1 != overhang_2:
                if self._verbose:
                    print("ERROR: Missmatch in 5' overhang (%s, %s)\n" % (overhang_1, overhang_2))
                return (None, None)
            
        return (cleavage_site_t, cleavage_site_b)
    
    def _find_best_cassette_end(self, cass, test, end):
        # Testing both orientations and taking best result
        cass_pos = self._find_cassette_end(cass, test, end)     
        cass_pos_rc = self._find_cassette_end(cass.reverse_complement(), test, end)

        if cass_pos is None and cass_pos_rc is not None:
            if self._verbose:
                print("            Best score = %.2f (reverse complement)" % cass_pos_rc.score)
            return (cass_pos_rc, True)

        elif cass_pos is not None and cass_pos_rc is None:
            if self._verbose:
                print("            Best score = %.2f (sense)" % cass_pos.score)
            return (cass_pos, False)

        elif cass_pos is None and cass_pos_rc is None:
            return (None, False)

        if cass_pos_rc.score > cass_pos.score:
            if self._verbose:
                print("            Best score = %.2f (reverse complement)" % cass_pos_rc.score)
            return (cass_pos_rc, True)
        else:
            if self._verbose:
                print("            Best score = %.2f (sense)" % cass_pos.score)
            return (cass_pos, False)
        
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
                  
        return max_alignment

    def _find_best_target_in_ref(self, ref, test, path, search_length, search_offset):
        max_alignment = Align.PairwiseAlignment(
            target="", query="", path=((0, 0), (0, 0)), score=0.0)

        max_isRC = False
        max_en = 0
        for en in range(len(path)):
            (alignment, isRC) = self._find_target_in_ref(ref, test, path[en][0]-search_offset, search_length)

            if alignment is None:
                continue

            if alignment.score > max_alignment.score:
                max_alignment = alignment
                max_isRC = isRC
                max_en = en

        if max_alignment.score == 0 or get_quality(max_alignment) < self._min_quality:
            return (None, False, 0)

        if self._verbose:
            if max_isRC:
                print("            Best score = %.2f (reverse complement)" % max_alignment.score)
            else:
                print("            Best score = %.2f (sense)" % max_alignment.score)

        return (max_alignment, max_isRC, path[max_en][0]-search_offset)

    def _find_target_in_ref(self, ref, test, pos, search_length):
        test_target = get_seq(test, pos, pos+search_length)
        if test_target is None:
            if self._verbose:
                print("            No test sequence found")
            return (None, False)

        # Finding test target in reference sequence
        alignments = self._aligner.align(ref, test_target)
        max_alignment = self._get_max_alignment(alignments)

        # If no matches were found, try the reverse complement of the test target
        alignments_rc = self._aligner.align(ref, test_target.reverse_complement())
        max_alignment_rc = self._get_max_alignment(alignments_rc)

        if max_alignment is None and max_alignment_rc is not None:
            return (max_alignment_rc, True)

        elif max_alignment is not None and max_alignment_rc is None:
            return (max_alignment, False)

        elif max_alignment is None and max_alignment_rc is None:
            return (None, False)

        if max_alignment_rc.score > max_alignment.score:
            return (max_alignment_rc, True)
        else:
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

def get_seq(seq, pos1, pos2):
    seq_len = len(seq)

    # Returning None if the output would require cycling over the sequence more than once
    if abs(pos1) > seq_len or abs(pos2) > 2*seq_len:
        return None

    # Getting sequence
    if pos1 >= 0 and pos2 < seq_len:
        # pos1 and pos2 are within the limits of seq
        return seq[pos1:pos2]

    elif pos1 < 0 and pos2 < seq_len:
        # Taking the first part of the sequence from the end
        seq1 = seq[pos1:]
        seq2 = seq[0:pos2]
        return seq1+seq2

    elif pos1 >= 0 and pos2 >= seq_len:
        # Taking the last part of the sequence from the start
        seq1 = seq[pos1::]
        seq2 = seq[0:pos2-seq_len]
        return seq1+seq2

    elif pos1 < 0 and pos2 >= seq_len:
        # Taking the last part of the sequence from the start
        seq1 = seq[pos1:]
        seq2 = seq[0:pos2-seq_len]
        return seq1+seq+seq2

def get_local_sequences(ref, cleavage_site_t, cleavage_site_b, local_r=1):
    if cleavage_site_t is None or cleavage_site_b is None:
        return ("", "")

    local_seq_t = ref[cleavage_site_t - local_r :cleavage_site_t + local_r]
    local_seq_b = ref[cleavage_site_b - local_r : cleavage_site_b + local_r].reverse_complement()
        
    return (local_seq_t, local_seq_b)

def get_quality(alignment):
    return alignment.score / len(alignment.query)

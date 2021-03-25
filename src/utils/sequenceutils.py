from Bio import Align
from matplotlib.pyplot import minorticks_off
from enums.ends import Ends
from enums.orientation import Orientation

import math

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
    def get_cleavage_positions(self, ref, cass, test, error_store=None):
        if self._verbose:
            print("        Finding first cassette end in test sequence:")
        (cass_pos_1, cass1_isRC) = self._find_best_cassette_end(cass, test, Ends.CASS_START, Orientation.BOTH)
                
        if self._verbose:
            print("        Finding second cassette end in test sequence:")
        orientation = Orientation.RC if cass1_isRC else Orientation.SENSE
        (cass_pos_2, cass2_isRC) = self._find_best_cassette_end(cass, test, Ends.CASS_END, orientation)

        # Both should be RC or normal.
        if cass1_isRC != cass2_isRC:
            if self._verbose:
                print("ERROR: Cassette RC mismatch\n")
            if error_store is not None:
                error_store.cassette_ends_mismatch()
            return (None, None, False)

        # Checking results are OK
        if cass_pos_1 is None or cass_pos_2 is None:
            if self._verbose:
                print("ERROR: Cassette ends not found\n")
            if error_store is not None:
                error_store.cassette_not_found_in_test()
            return (None, None, False)  

        # Finding cassette-adjacent test sequence in reference
        if self._verbose:
            print("        Finding first cassette-adjacent test sequence in reference sequence:")
        (alignment1, isRC1, cass_start) = self._find_best_target_in_ref(ref, test, cass_pos_1.path[0:-1], self._num_bases, self._num_bases, self._min_quality)
                
        if self._verbose:
            print("        Finding second cassette-adjacent test sequence in reference sequence:")
        (alignment2, isRC2, cass_end) = self._find_best_target_in_ref(ref, test, cass_pos_2.path[1:], self._num_bases, 0, self._min_quality)
        
        if alignment1 is None or alignment2 is None:
            if self._verbose:
                print("ERROR: Test sequence not found in reference\n")
            if error_store is not None:
                error_store.test_not_found_in_reference()
            return (None, None, False)

        # Both should be RC or normal
        if isRC1 != isRC2:
            if self._verbose:
                print("ERROR: RC mismatch\n")
            if error_store is not None:
                error_store.test_ends_mismatch()
            return (None, None, False)

        if isRC1:            
            cleavage_site_t = alignment1.path[0][0]
            cleavage_site_b = alignment2.path[-1][0]
        else:
            cleavage_site_t = alignment2.path[0][0]
            cleavage_site_b = alignment1.path[-1][0]

        if abs(cleavage_site_t - cleavage_site_b) > self._max_gap:
            if self._verbose:
                print("ERROR: cleavage site gap (%i) exceeds maximum permitted\n" % abs(cleavage_site_b - cleavage_site_t))
            if error_store is not None:
                error_store.max_gap_exceeded()
            return (None, None, False)
            
        if self._verbose:
            print("        Testing for sequence splitting:")
        midpoint_site = self.get_midpoint_position(ref, test, cass_start, cass_end)     
        if midpoint_site is None:
            if self._verbose:
                print("ERROR: Midpoint sequence not found\n")
            if error_store is not None:
                error_store.midpoint_not_found()
            return (None, None, False)

        split = sample_is_split(cleavage_site_t, cleavage_site_b, midpoint_site)

        return (cleavage_site_t, cleavage_site_b, split)

    def get_midpoint_position(self, ref, test, cass_start, cass_end):
        if cass_start < cass_end:
            test_mid_pos = int((len(test)+cass_end+cass_start)/2 % len(test))
        else:
            test_mid_pos = int((cass_start-cass_end)/2 + cass_end)

        # Finding position of midpoint sequence in reference (this is the midpoint position)
        (midpoint, isRC) = self._find_target_in_ref(ref, test, test_mid_pos, self._num_bases, 0.75)
        if midpoint is None:
            return None  

        if self._verbose:
                print("            Best score = %.2f (reverse complement)" % midpoint.score)

        return midpoint.path[0][0]

    def _find_best_cassette_end(self, cass, test, end, orientation=Orientation.BOTH):
        # Testing cassette in sense orientation
        if orientation is Orientation.SENSE or orientation is Orientation.BOTH:
            cass_pos = self._find_cassette_end(cass, test, end)     
        else:
            cass_pos = None
        
        # Testing cassette in reverse complement orientation
        if orientation is Orientation.RC or orientation is Orientation.BOTH:
            cass_pos_rc = self._find_cassette_end(cass.reverse_complement(), test, end)
        else:
            cass_pos_rc = None

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

        # Adding a repetition to the end of the sequence, incase the target sequence spans the ends
        len_test = len(test)
        test = test + get_seq(test, 0, self._num_bases)
        
        if end is Ends.CASS_START:
            # Checking for "full" cassette end
            alignments = self._aligner.align(test, cass[0: self._num_bases])
            max_alignment = self._get_max_alignment(alignments,self._min_quality)

        elif end is Ends.CASS_END:
            # Checking for "full" cassette end
            alignments = self._aligner.align(test, cass[-self._num_bases ::])
            max_alignment = self._get_max_alignment(alignments,self._min_quality)

        if max_alignment is None:
            return None
                  
        max_alignment.path = remove_path_rollover(max_alignment.path,len_test)

        return max_alignment

    def _find_best_target_in_ref(self, ref, test, path, search_length, search_offset, min_quality):
        max_alignment = Align.PairwiseAlignment(
            target="", query="", path=((0, 0), (0, 0)), score=0.0)

        max_isRC = False
        max_en = 0
        for en in range(len(path)): # The last position in the path is definitely the end
            (alignment, isRC) = self._find_target_in_ref(ref, test, path[en][0]-search_offset, search_length, min_quality)

            if alignment is None:
                continue

            if alignment.score > max_alignment.score:
                max_alignment = alignment
                max_isRC = isRC
                max_en = en

        if max_alignment.score == 0 or get_quality(max_alignment) < min_quality:
            return (None, False, 0)

        if self._verbose:
            if max_isRC:
                print("            Best score = %.2f (reverse complement)" % max_alignment.score)
            else:
                print("            Best score = %.2f (sense)" % max_alignment.score)

        return (max_alignment, max_isRC, path[max_en][0])

    def _find_target_in_ref(self, ref, test, pos, search_length, min_quality):
        test_target = get_seq(test, pos, pos+search_length)
        # print(test_target)
        
        # Adding a repetition to the end of the sequence, incase the target sequence spans the ends
        len_ref = len(ref)
        ref = ref + get_seq(ref, 0, search_length)
        
        if test_target is None:
            if self._verbose:
                print("            No test sequence found")
            return (None, False)

        # Finding test target in reference sequence
        alignments = self._aligner.align(ref, test_target)
        max_alignment = self._get_max_alignment(alignments,min_quality)

        # If no matches were found, try the reverse complement of the test target
        alignments_rc = self._aligner.align(ref, test_target.reverse_complement())
        max_alignment_rc = self._get_max_alignment(alignments_rc,min_quality)

        if max_alignment is None and max_alignment_rc is not None:
            max_alignment_rc.path = remove_path_rollover(max_alignment_rc.path,len_ref)
            return (max_alignment_rc, True)

        elif max_alignment is not None and max_alignment_rc is None:
            max_alignment.path = remove_path_rollover(max_alignment.path,len_ref)
            return (max_alignment, False)

        elif max_alignment is None and max_alignment_rc is None:
            return (None, False)

        if max_alignment_rc.score > max_alignment.score:
            max_alignment_rc.path = remove_path_rollover(max_alignment_rc.path,len_ref)
            return (max_alignment_rc, True)
        else:
            max_alignment.path = remove_path_rollover(max_alignment.path,len_ref)
            return (max_alignment, False)

    def _get_max_alignment(self, alignments, min_quality):
        max_alignment = Align.PairwiseAlignment(
            target="", query="", path=((0, 0), (0, 0)), score=0.0)

        for alignment in alignments:
            if alignment.score > max_alignment.score and get_quality(alignment) >= min_quality:
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

    local_seq_t = get_seq(ref,cleavage_site_t - local_r, cleavage_site_t + local_r)
    local_seq_b = get_seq(ref, cleavage_site_b - local_r, cleavage_site_b + local_r).reverse_complement()
        
    return (local_seq_t, local_seq_b)

def get_event_width(cleavage_site_t, cleavage_site_b, split, seq_len):
    if split:
        return seq_len-(max(cleavage_site_b,cleavage_site_t)-min(cleavage_site_b,cleavage_site_t))
    else:
        return max(cleavage_site_b,cleavage_site_t)-min(cleavage_site_b,cleavage_site_t)

def get_quality(alignment):
    return alignment.score / len(alignment.query)

def get_sequence_str(ref, cleavage_site_t, cleavage_site_b, split, extra_nt=0):
    if cleavage_site_b is None or cleavage_site_t is None:
        return (None, None)

    event_width = get_event_width(cleavage_site_t, cleavage_site_b, split, len(ref))

    # Identifying break type
    if (cleavage_site_b < cleavage_site_t and not split) or (cleavage_site_b > cleavage_site_t and split):
        # 3' overhang
        left_seq1 = get_seq(ref,cleavage_site_b - 1 - extra_nt, cleavage_site_b)
        mid_seq1 = get_seq(ref,cleavage_site_b,cleavage_site_b+event_width)
        right_seq1 = get_seq(ref,cleavage_site_t, cleavage_site_t + 1 + extra_nt)

        left_seq2 = left_seq1.complement()
        mid_seq2 = mid_seq1.complement()
        right_seq2 = right_seq1.complement()

        seq1 = str("5'...%s %s|%s...3'" % (left_seq1, mid_seq1, right_seq1))
        seq2 = str("3'...%s|%s %s...5'" % (left_seq2, mid_seq2, right_seq2))

        return (seq1,seq2)
        
    elif cleavage_site_b == cleavage_site_t:
        # Blunt end
        left_seq1 = get_seq(ref, cleavage_site_b - 3 - extra_nt, cleavage_site_b)
        right_seq1 = get_seq(ref, cleavage_site_t, cleavage_site_t + 3 + extra_nt)

        left_seq2 = left_seq1.complement()
        right_seq2 = right_seq1.complement()

        seq1 = str("5'...%s|%s...3'" % (left_seq1, right_seq1))
        seq2 = str("3'...%s|%s...5'" % (left_seq2, right_seq2))
        
        return (seq1, seq2)

    elif (cleavage_site_b > cleavage_site_t and not split) or (cleavage_site_b < cleavage_site_t and split):
        # 5' overhang
        left_seq1 = get_seq(ref, cleavage_site_t - 1 - extra_nt, cleavage_site_t)
        mid_seq1 = get_seq(ref, cleavage_site_t, cleavage_site_t+event_width)
        right_seq1 = get_seq(ref, cleavage_site_b, cleavage_site_b + 1 + extra_nt)

        left_seq2 = left_seq1.complement()
        mid_seq2 = mid_seq1.complement()
        right_seq2 = right_seq1.complement()

        seq1 = str("5'...%s|%s %s...3'" % (left_seq1, mid_seq1, right_seq1))
        seq2 = str("3'...%s %s|%s...5'" % (left_seq2, mid_seq2, right_seq2))

        return (seq1, seq2)

def get_type_str(cleavage_site_t, cleavage_site_b, split):
    if cleavage_site_b is None or cleavage_site_t is None:
        return None

    # Identifying break type
    if (cleavage_site_b < cleavage_site_t and not split) or (cleavage_site_b > cleavage_site_t and split):
        return "3' overhang"
        
    elif cleavage_site_b == cleavage_site_t:
        return "Blunt end"

    elif (cleavage_site_b > cleavage_site_t and not split) or (cleavage_site_b < cleavage_site_t and split):
        return "5' overhang"

def sample_is_split(cleavage_site_t, cleavage_site_b, midpoint_site):
    return (cleavage_site_t < midpoint_site and cleavage_site_b > midpoint_site) or (cleavage_site_b < midpoint_site and cleavage_site_t > midpoint_site)

def remove_path_rollover(path,seq_length):
    path = list(path)
    for i in range(len(path)):
        path[i] = (path[i][0] % seq_length,path[i][1])
    
    return tuple(path)
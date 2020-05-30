from Bio import Align
from ends import Ends


class SequenceSearcher():
    def __init__(self, aligner, verbose=False):
        self._aligner = aligner
        self._verbose = verbose

    def get_aligner(self):
        return self._aligner

    def set_aligner(self, aligner):
        self._aligner = aligner

    def get_verbose(self):
        return self._verbose

    def set_verbose(self, verbose):
        self._verbose = verbose

    def _find_cassette_end(self, cass, test, num_bases=20):
        # Finding cassette ends in test sequence
        if self._verbose:
            print("    Finding cassette end in test sequence")

        max_alignment = Align.PairwiseAlignment(
            target="", query="", path=((0, 0), (0, 0)), score=0.0)
        end = 0

        # Testing against cassette start (RC)
        cass_start_rc = cass[0:num_bases].reverse_complement()
        alignments = self._aligner.align(test, cass_start_rc)
        for alignment in alignments:
            if alignment.score > max_alignment.score:
                max_alignment = alignment
                end = Ends.CASS_START_RC

        # Testing against cassette end
        cass_end = cass[-num_bases::]
        alignments = self._aligner.align(test, cass_end)
        for alignment in alignments:
            if alignment.score > max_alignment.score:
                max_alignment = alignment
                end = Ends.CASS_END

        if self._verbose:
            if end == Ends.CASS_START_RC:
                pos_string = "start RC"
            elif end == Ends.CASS_END:
                pos_string = "end"

            print("        Best match for cassette %s (%s)" %
                (pos_string, max_alignment.query))
            print("        Match score = %0.2f (quality %0.2f)\n" %
                (max_alignment.score, get_quality(max_alignment)))

        return (max_alignment.path[-1][0], end, max_alignment.score)

    def _find_target_in_ref(self, ref, test, end, path, num_bases=20, min_quality=1.0):
        (alignments, isRC) = self._find_all_targets_in_ref(
            ref, test, end, path, num_bases, min_quality)

        # If more than one match was found, increasing the number of bases used in the recognition
        while len(alignments) > 1:
            num_bases = num_bases + 1
            (alignments, isRC) = self._find_target_in_ref(
                ref, test, end, path, num_bases, min_quality)

        if len(alignments) == 0:
            print("\nNo match found\n")
            return (None, False)

        return (alignments[0], isRC)

    def _find_all_targets_in_ref(self, ref, test, end, pos, num_bases, min_quality):
        test_target = test[pos: pos + num_bases]

        # Finding test target in reference sequence
        alignments = self._aligner.align(ref, test_target)
        valid_alignments = self._get_valid_alignments(alignments, min_quality)

        # If no matches were found, try the reverse complement of the test target
        isRC = False
        if len(valid_alignments) == 0:
            if self._verbose:
                print("    Reversing test target sequence")

            isRC = True
            test_target = test_target.reverse_complement()
            alignments = self._aligner.align(ref, test_target)
            valid_alignments = self._get_valid_alignments(
                alignments, min_quality)

        if self._verbose:
            print("    Testing reference with %i bases (%i matches found)" %
                  (num_bases, len(valid_alignments)))

        return (valid_alignments, isRC)

    def _get_valid_alignments(self, alignments, min_quality):
        valid_alignments = []

        for alignment in alignments:
            quality = get_quality(alignment)
            if quality >= min_quality:
                valid_alignments.append(alignment)

        return valid_alignments

    def find_clevage_site(self, ref, cass, test, num_bases=20, min_quality=1.0, break_range=5):
        (pos, end, score) = self._find_cassette_end(cass, test, num_bases=num_bases)
        
        # Getting region of test sequence to match to reference
        if self._verbose:
            print("    Finding cassette-adjacent sequence in reference")
            
        (alignment, isRC) = self._find_target_in_ref(ref, test, end, pos, num_bases=num_bases, min_quality=min_quality)

        if alignment is None:
            return

        if self._verbose:
            print("        Match score = %0.2f (quality %0.2f)" %
                (alignment.score, get_quality(alignment)))
        if isRC:
            clevage_position = alignment.path[-1][0]
        else:
            clevage_position = alignment.path[0][0]

        ref_break_seq_left = ref[clevage_position - break_range:clevage_position]
        ref_break_seq_right = ref[clevage_position: clevage_position + break_range]

        if self._verbose:
            print("        Reference break at position %i" % clevage_position)
        
        if self._verbose:
            print("        Reference break at sequence %s | %s\n" %
                (ref_break_seq_left, ref_break_seq_right))

        return (clevage_position, isRC)
            


def get_quality(alignment):
    return alignment.score / len(alignment.query)

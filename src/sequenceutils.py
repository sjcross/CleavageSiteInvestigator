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

    def find_target_in_ref(self, ref, test, end, path, num_bases=20, min_quality=1.0):
        (alignments, isRC) = self._find_target_in_ref(
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

    def _find_target_in_ref(self, ref, test, end, path, num_bases, min_quality):
        pos = path[-1][0]
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


def get_quality(alignment):
    return alignment.score / len(alignment.query)

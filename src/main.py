# Imports
from Bio import Align
from ends import Ends
from filehandling import FileReader

import sequenceutils as su

# Setting some parameters
ref_seq_path = "C:\\Users\\steph\\Desktop\\Mark Szczelkun\\"
ref_seq_name = "pRMA03+L2L2.dna"
cass_seq_path = "C:\\Users\\steph\\Desktop\\Mark Szczelkun\\"
cass_seq_name = "CAT cassette as amplified by RA101 & 102 from pACYC184.dna"
test_seq_path = "C:\\Users\\steph\\Desktop\\Mark Szczelkun\\Sequencing files\\Run250808-10\\"
test_seq_name = "48_041.ab1"

cass_end_length = 20
break_range = 5
verbose = True

# Creating FileHandler object
filereader = FileReader(verbose=verbose)

# Loading reference, cassette and test sequences
print("Loading sequences from file")
ref = filereader.read_sequence(ref_seq_path, ref_seq_name)
cass = filereader.read_sequence(cass_seq_path, cass_seq_name)
test = filereader.read_sequence(test_seq_path, test_seq_name)
print("\r")

# Creating the PairwiseAligner and SequenceSearcher objects
aligner = Align.PairwiseAligner()
aligner.mode = 'local'
aligner.match_score = 1.0
aligner.mismatch_score = -1.0
aligner.gap_score = -1.0
searcher = su.SequenceSearcher(aligner, verbose=verbose)

# Finding cassette ends in test sequence
print("Finding cassette end in test sequence")
max_alignment = Align.PairwiseAlignment(
    target="", query="", path=((0, 0), (0, 0)), score=0.0)
end = 0

# Testing against cassette start (RC)
cass_start_rc = cass[0:cass_end_length].reverse_complement()
alignments = aligner.align(test, cass_start_rc)
for alignment in alignments:
    if alignment.score > max_alignment.score:
        max_alignment = alignment
        end = Ends.CASS_START_RC

# Testing against cassette end
cass_end = cass[-cass_end_length::]
alignments = aligner.align(test, cass_end)
for alignment in alignments:
    if alignment.score > max_alignment.score:
        max_alignment = alignment
        end = Ends.CASS_END

if end == Ends.CASS_START_RC:
    pos_string = "start RC"
elif end == Ends.CASS_END:
    pos_string = "end"

print("    Best match for cassette %s (%s)" %
      (pos_string, max_alignment.query))
print("    Match score = %0.2f (quality %0.2f)\n" %
      (max_alignment.score, su.get_quality(max_alignment)))

# Getting region of test sequence to match to reference
print("Finding cassette-adjacent sequence in reference")
num_bases = 20
min_quality = 0.75
(alignment, isRC) = searcher.find_target_in_ref(ref, test, end,
                                                max_alignment.path, num_bases=num_bases, min_quality=min_quality)

if alignment is None:
    quit()

print("    Match score = %0.2f (quality %0.2f)" %
      (alignment.score, su.get_quality(alignment)))

if isRC:
    break_position = alignment.path[-1][0]
else:
    break_position = alignment.path[0][0]

ref_break_seq_left = ref[break_position - break_range:break_position]
ref_break_seq_right = ref[break_position: break_position + break_range]

print("    Reference break at position %i" % break_position)
print("    Reference break at sequence %s | %s\n" %
      (ref_break_seq_left, ref_break_seq_right))

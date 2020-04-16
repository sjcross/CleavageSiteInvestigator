# Imports
from Bio import Align
from filereader import FileReader

# Setting some parameters
ref_seq_path = "C:\\Users\\steph\\OneDrive\\Desktop\\Mark Szczelkun\\"
ref_seq_name = "pRMA03+L2L2.dna"
cass_seq_path = "C:\\Users\\steph\\OneDrive\\Desktop\\Mark Szczelkun\\"
cass_seq_name = "CAT cassette as amplified by RA101 & 102 from pACYC184.dna"
test_seq_path = "C:\\Users\\steph\\OneDrive\\Desktop\\Mark Szczelkun\\Sequencing files\\Run250808-04\\"
test_seq_name = "01_008.ab1"

cass_end_length = 20

# Creating FileHandler object
filereader = FileReader()
filereader.set_verbose(True)

# Loading reference, cassette and test sequences
ref = filereader.read_sequence(ref_seq_path, ref_seq_name)
# print("Reference sequence:\n%s" % ref)
print("\n")

cass = filereader.read_sequence(cass_seq_path, cass_seq_name)
# print("Cassette sequence:\n%s" % cass)
print("\n")

test = filereader.read_sequence(test_seq_path, test_seq_name)
# print("Test sequence:\n%s" % test)
print("\n")

# Finding cassette ends in test sequence
# Creating the PairwiseAligner object
aligner = Align.PairwiseAligner()
aligner.mode = 'local'
aligner.match_score = 1.0
aligner.mismatch_score = -1.0
aligner.gap_score = -1.0

max_alignment = Align.PairwiseAlignment(
    target="", query="", path=((0, 0), (0, 0)), score=0.0)
orientation = 0

# Testing against cassette start
cass_start = cass[0:cass_end_length]
alignments = aligner.align(test, cass_start)
for alignment in alignments:
    if alignment.score > max_alignment.score:
        max_alignment = alignment
        orientation = 1

# Testing against cassette start (RC)
cass_start_rc = cass[0:cass_end_length].reverse_complement()
alignments = aligner.align(test, cass_start_rc)
for alignment in alignments:
    if alignment.score > max_alignment.score:
        max_alignment = alignment
        orientation = 2

# Testing against cassette end
cass_end = cass[-cass_end_length::]
alignments = aligner.align(test, cass_end)
for alignment in alignments:
    if alignment.score > max_alignment.score:
        max_alignment = alignment
        orientation = 3

# Testing against cassette end (RC)
cass_end_rc = cass[-cass_end_length::].reverse_complement()
alignments = aligner.align(test, cass_end_rc)
for alignment in alignments:
    if alignment.score > max_alignment.score:
        max_alignment = alignment
        orientation = 4

# SHOULD BE POSSIBLE TO SWAP ORIENTATION TO BE 1 OR -1 TO INDIATE IF WE SHOULD BE
# TAKING THE SEQUENCE AFTER OR BEFORE THE MATCH AS THE END
if orientation == 1:
    print("Best match for cassette start (score = %d)" % max_alignment.score)
elif orientation == 2:
    print("Best match for cassette start RC (score = %d)" % max_alignment.score)
elif orientation == 3:
    print("Best match for cassette end (score = %d)" % max_alignment.score)
elif orientation == 4:
    print("Best match for cassette end RC (score = %d)" % max_alignment.score)

print("Matched sequence: %s" % max_alignment.query)
print("Alignment quality = %f" % (max_alignment.score / cass_end_length))
print("\n")

# Getting region of test sequence to match to reference
if orientation == 3:
    pos = max_alignment.path[1][0]
    test_target = test[pos:pos+20].reverse_complement()

print("Test target: %s" % test_target)

# Finding test target in reference sequence
alignments = aligner.align(ref, test_target)
for alignment in alignments:
    print(alignment.score)

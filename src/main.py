# Imports
from Bio import Align, Seq
from ends import Ends
from filehandling import FileReader

import sequenceutils as su

# Setting some parameters
ref_seq_path = "C:\\Users\\Stephen\\Desktop\\Users\\Mark Szczelkun\\"
ref_seq_name = "pRMA03+L2L2.dna"

cass_seq_path = "C:\\Users\\Stephen\\Desktop\\Users\\Mark Szczelkun\\"
cass_seq_name = "CAT cassette as amplified by RA101 & 102 from pACYC184.dna"

# KpnI (3' overhang) examples
test1_seq_path = "C:\\Users\\Stephen\\Desktop\\Users\\Mark Szczelkun\\Sequencing files\\Run220808-04\\"
test1_seq_name = "45_044.ab1"
test2_seq_path = "C:\\Users\\Stephen\\Desktop\\Users\\Mark Szczelkun\\Sequencing files\\Run250808-04\\"
test2_seq_name = "02_007.ab1"

# # SmaI (blunt ends) examples
# test1_seq_path = "C:\\Users\\Stephen\\Desktop\\Users\\Mark Szczelkun\\Sequencing files\\Run210808-04\\"
# test1_seq_name = "40_033.ab1"
# test2_seq_path = "C:\\Users\\Stephen\\Desktop\\Users\\Mark Szczelkun\\Sequencing files\\Run200808-06\\"
# test2_seq_name = "45_044.ab1"

# XhoI (5' overhang) examples
# test1_seq_path = "C:\\Users\\Stephen\\Desktop\\Users\\Mark Szczelkun\\Sequencing files\\Run250808-10\\"
# test1_seq_name = "44_045.ab1"
# test2_seq_path = "C:\\Users\\Stephen\\Desktop\\Users\\Mark Szczelkun\\Sequencing files\\Run250808-10\\"
# test2_seq_name = "48_041.ab1"

cass_end_length = 20
break_range = 5
num_bases = 20
min_quality = 0.75
verbose = True

# Creating FileHandler object
filereader = FileReader(verbose=verbose)

# Loading reference, cassette and test sequences
print("Loading sequences from file")
ref = filereader.read_sequence(ref_seq_path, ref_seq_name)
cass = filereader.read_sequence(cass_seq_path, cass_seq_name)
test1 = filereader.read_sequence(test1_seq_path, test1_seq_name)
test2 = filereader.read_sequence(test2_seq_path, test2_seq_name)
print("\r")

# Creating the PairwiseAligner and SequenceSearcher objects
aligner = Align.PairwiseAligner()
aligner.mode = 'local'
aligner.match_score = 1.0
aligner.mismatch_score = -1.0
aligner.gap_score = -1.0
searcher = su.SequenceSearcher(aligner, verbose=verbose)

# Finding break in test sequence 1
print("Finding break in test sequence 1")
(clevage_site1, isRC1) = searcher.find_clevage_site(ref, cass, test1, num_bases=num_bases, min_quality=min_quality)

# Finding break in test sequence 2
print("Finding break in test sequence 2")
(clevage_site2, isRC2) = searcher.find_clevage_site(ref, cass, test2, num_bases=num_bases, min_quality=min_quality)

# Only one should be RC
if isRC1 and isRC2:
    print("ERROR: Both identified as reverse complement")
    quit()
elif not isRC1 and not isRC2:
    print("ERROR: Neither identified as reverse complement")
    quit()

# If clevage_site1 is RC, switch clevage sites
if isRC2:
    (clevage_site1, clevage_site2) = (clevage_site2, clevage_site1)
    
# Identifying break type
if clevage_site1 < clevage_site2:
    # 3' overhang
    left_seq1 = ref[clevage_site1 - 1 :clevage_site1]
    mid_seq1 = ref[clevage_site1:clevage_site2]
    right_seq1 = ref[clevage_site2: clevage_site2 + 1]

    left_seq2 = left_seq1.complement()
    mid_seq2 = mid_seq1.complement()
    right_seq2 = right_seq1.complement()

    print("Restriction site (3' overhang):")
    print("    5'...%s %s↓%s...3'\r\n    3'...%s↑%s %s...5'\r\n" % (left_seq1, mid_seq1, right_seq1, left_seq2, mid_seq2, right_seq2))
    
elif (clevage_site1 == clevage_site2):
    # Blunt end
    # 3' overhang
    left_seq1 = ref[clevage_site1 - 3 :clevage_site1]
    right_seq1 = ref[clevage_site2: clevage_site2 + 3]

    left_seq2 = left_seq1.complement()
    right_seq2 = right_seq1.complement()

    print("Restriction site (blunt end):")
    print("    5'...%s↓%s...3'\r\n    3'...%s↑%s...5'\r\n" % (left_seq1, right_seq1, left_seq2, right_seq2))

elif clevage_site1 > clevage_site2:
    # 5' overhang
    left_seq1 = ref[clevage_site2 - 1 :clevage_site2]
    mid_seq1 = ref[clevage_site2:clevage_site1]
    right_seq1 = ref[clevage_site1: clevage_site1 + 1]

    left_seq2 = left_seq1.complement()
    mid_seq2 = mid_seq1.complement()
    right_seq2 = right_seq1.complement()

    print("Restriction site (5' overhang):")
    print("    5'...%s↓%s %s...3'\r\n    3'...%s %s↑%s...5'\r\n" % (left_seq1, mid_seq1, right_seq1, left_seq2, mid_seq2, right_seq2))

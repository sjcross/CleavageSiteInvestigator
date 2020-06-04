### Imports ###
from Bio import Align, Seq
from ends import Ends
from seqtype import Seqtype
from filehandling import FileReader

import os
import sequenceutils as su


### Parameters ###
### ROLLING CIRCLE EXAMPLES ###
# Root folder containing all files
root_folder = "D:\\Stephen\\Users\\Mark Szczelkun\\2020-04-28 New files\\"

# The sequence for the plasmid into which the cassette has been inserted
ref_seq_path = root_folder
ref_seq_name = "pUC19.fa"

# The cassette sequence
cass_seq_path = root_folder
cass_seq_name = "Chloramphenicol Cassette overhang.fa"

# The sequencing result
test1_seq_path = root_folder
test1_seq_name = "pUC_PstI_12.fasta"


### SANGER EXAMPLES
# # Root folder containing all files
# root_folder = "D:\\Stephen\\Users\\Mark Szczelkun\\"

# # The sequence for the plasmid into which the cassette has been inserted
# ref_seq_path = root_folder
# ref_seq_name = "pRMA03+L2L2.dna"

# # The cassette sequence
# cass_seq_path = root_folder
# cass_seq_name = "CAT cassette as amplified by RA101 & 102 from pACYC184.dna"

# # The first sequencing result
# test1_seq_path = root_folder + "Sequencing files\\Run220808-04\\"
# test1_seq_name = "45_044.ab1"

# # (Optional) The second sequencing result (must be opposite primer to first sequencing result)
# test2_seq_path = root_folder + "Sequencing files\\Run250808-04\\"
# test2_seq_name = "02_007.ab1"


seq_type = Seqtype.OTHER # Sequencing type (must be either SANGER or OTHER)
num_bases = 20 # Number of bases to match
min_quality = 0.75 # Minimum match quality ("1" is perfect)
verbose = False # Display messages during execution


### Processing ###
# Creating FileHandler object
filereader = FileReader(verbose=verbose)

# Loading reference, cassette and test sequences
print("INPUT: Loading sequences from file")
ref = filereader.read_sequence(ref_seq_path, ref_seq_name)[0]
cass = filereader.read_sequence(cass_seq_path, cass_seq_name)[0]
tests1 = filereader.read_sequence(test1_seq_path, test1_seq_name)
if seq_type is Seqtype.SANGER:
    tests2 = filereader.read_sequence(test2_seq_path, test2_seq_name)
    tests = zip(tests1, tests2)
else:
    tests = tests1
print("\r")

# Creating the PairwiseAligner and SequenceSearcher objects
aligner = Align.PairwiseAligner()
aligner.mode = 'local'
aligner.match_score = 1.0
aligner.mismatch_score = -1.0
aligner.gap_score = -1.0
searcher = su.SequenceSearcher(aligner, min_quality=min_quality, num_bases=num_bases, verbose=verbose)

if seq_type is Seqtype.SANGER:
    print("PROCESSING: Sanger sequence(s)")
    for count, (test1, test2) in enumerate(tests):
        print("Processing test sequence %i" % (count+1))
        (clevage_site1, clevage_site2) = searcher.process_sanger(ref, cass, test1, test2)
        su.print_clevage_site(ref, clevage_site1, clevage_site2)
    
elif seq_type is Seqtype.OTHER:
    print("PROCESSING: \"Other\" sequence(s)")
    for count, test in enumerate(tests):
        print("Processing test sequence %i" % (count+1))
        (clevage_site1, clevage_site2) = searcher.process_other(ref, cass, test)
        su.print_clevage_site(ref, clevage_site1, clevage_site2)

### Imports ###
import os

from Bio import Align, Seq
from enums.ends import Ends
from enums.seqtype import Seqtype
from utils.fileutils import FileReader

from utils import reportutils as ru
from utils import sequenceutils as su


### Parameters ###
### ROLLING CIRCLE EXAMPLES ###
# Root folder containing all files
# root_folder = "D:\\Stephen\\Users\\Mark Szczelkun\\"
root_folder = "F:\\People\\Mark Szczelkun\\"

# The sequence for the plasmid into which the cassette has been inserted
ref_seq_path = root_folder + "2020-04-28 New files\\"
ref_seq_name = "pUC19.fa"

# The cassette sequence
cass_seq_path = root_folder + "2020-04-28 New files\\"
cass_seq_name = "Chloramphenicol Cassette overhang.fa"

# The sequencing result
test_seq_path = root_folder + "2020-06-04 Mix files\\"
test_seq_name = "XbaI_R2C2_Consensus_fix.fasta"


local_r = 1 # Half width of the local sequences to be extracted at restriction sites
max_gap = 10 # Maximum number of bp between 3' and 5' restriction sites
min_quality = 0.75  # Minimum match quality ("1" is perfect)
num_bases = 20  # Number of bases to match
verbose = False # Display messages during execution

### Processing ###
# Creating FileHandler object
filereader = FileReader()

# Loading reference, cassette and test sequences
print("INPUT: Loading sequences from file")

ref = filereader.read_sequence(ref_seq_path, ref_seq_name)[0]
cass = filereader.read_sequence(cass_seq_path, cass_seq_name)[0]
tests = filereader.read_sequence(test_seq_path, test_seq_name)

print("\r")

# Creating the PairwiseAligner and SequenceSearcher objects
aligner = Align.PairwiseAligner()
aligner.mode = 'local'
aligner.match_score = 1.0
aligner.mismatch_score = -1.0
aligner.gap_score = -1.0
searcher = su.SequenceSearcher(aligner, local_r=local_r,max_gap=max_gap, min_quality=min_quality, num_bases=num_bases, verbose=verbose)

# Dict to store results as dual clevage site tuple
results = {}
error_count = 0

print("PROCESSING: Sequence(s)")

for count, test in enumerate(tests):
    if verbose:
        print("    Processing test sequence %i" % (count + 1))

    (clevage_site_t, clevage_site_b, local_seq_t, local_seq_b) = searcher.process(ref, cass, test)
    
    if clevage_site_t == None:
        error_count = error_count + 1
        continue

    k = (clevage_site_t, clevage_site_b)
    if k not in results:
        results[(clevage_site_t, clevage_site_b)] = 1
    else:
        results[(clevage_site_t, clevage_site_b)] = results[(clevage_site_t, clevage_site_b)] + 1

    if verbose:
        print("        Result:")
        ru.print_position(clevage_site_t, clevage_site_b, offset="        ")
        ru.print_type(clevage_site_t, clevage_site_b, offset="        ")
        ru.print_sequence(ref, clevage_site_t, clevage_site_b, offset="        ")
print("\r")

print("RESULTS:")
# Sorting results by frequency
results = ru.sort_results(results)

# Displaying results
for result in results.keys():
    clevage_site_t = result[0]
    clevage_site_b = result[1]
    count = results.get(result)

    ru.print_position(clevage_site_t, clevage_site_b)
    ru.print_count(count)
    ru.print_type(clevage_site_t, clevage_site_b)
    ru.print_sequence(ref, clevage_site_t, clevage_site_b)

ru.print_error_rate(error_count,len(tests), offset="    ")
print("\r")

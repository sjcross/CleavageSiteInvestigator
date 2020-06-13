### Imports ###
import os

from Bio import Align, Seq
from enums.ends import Ends
from enums.seqtype import Seqtype
from utils.fileutils import FileReader

from utils import plotutils as pu
from utils import reportutils as ru
from utils import sequenceutils as su


### Parameters ###
# Root folder containing all files
root_folder = "F:\\People\\Mark Szczelkun\\"

# The sequence for the plasmid into which the cassette has been inserted
ref_seq_path = root_folder + "2020-04-28 New files\\"
ref_seq_name = "pUC19.fa"

# The cassette sequence
cass_seq_path = root_folder + "2020-04-28 New files\\"
cass_seq_name = "Chloramphenicol Cassette overhang.fa"

# The sequencing result
test_seq_path = root_folder + "2020-06-04 Mix files\\"
test_seq_name = "Merge_test.fasta"


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

# Dict to store results as dual cleavage site tuple
results = {}
error_count = 0

print("PROCESSING: Sequence(s)")

for count, test in enumerate(tests):
    if verbose:
        print("    Processing test sequence %i" % (count + 1))

    (cleavage_site_t, cleavage_site_b, local_seq_t, local_seq_b) = searcher.process(ref, cass, test)
    
    if cleavage_site_t == None:
        error_count = error_count + 1
        continue
    
    results[count] = (cleavage_site_t, cleavage_site_b, local_seq_t, local_seq_b)

    if verbose:
        print("        Result:")
        ru.print_position(cleavage_site_t, cleavage_site_b, offset="        ")
        ru.print_type(cleavage_site_t, cleavage_site_b, offset="        ")
        ru.print_sequence(ref, cleavage_site_t, cleavage_site_b, offset="        ")
print("\r")

print("RESULTS:")
# Reporting full sequence frequency
freq_full = ru.get_full_sequence_frequency(results)
if verbose:
    print("    Full sequence frequencies:\n")
ru.print_full_sequence_frequency(ref, freq_full, offset="")

# Reporting local sequence frequency
print("    Local dinucleotide frequencies:\n")
freq_local = ru.get_local_sequence_frequency(results, ru.StrandMode.BOTH, ru.LocalMode.BOTH)
ru.print_local_sequence_frequency(freq_local, nonzero_only=False, offset="")

print("    Local 5' nucleotide frequencies:\n")
freq_5p = ru.get_local_sequence_frequency(results, ru.StrandMode.BOTH, ru.LocalMode.FIVE_P)
ru.print_local_sequence_frequency(freq_5p, nonzero_only=False, offset="")

print("    Local 3' nucleotide frequencies:\n")
freq_3p = ru.get_local_sequence_frequency(results, ru.StrandMode.BOTH, ru.LocalMode.THREE_P)
ru.print_local_sequence_frequency(freq_3p, nonzero_only=False, offset="")

# Plotting sequence distributions
pu.plotFrequency1D(freq_local, freq_5p, freq_3p)

# Reporting top and bottom sequence co-occurrence
(labels, freq2D) = ru.get_sequence_cooccurrence(results)
pu.plotFrequency2D(labels, freq2D)

# Reporting number of errors
ru.print_error_rate(error_count,len(tests), offset="    ")
print("\r")

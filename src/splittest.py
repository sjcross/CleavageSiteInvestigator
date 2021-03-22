### IMPORTS ###
import os

from Bio import Align
from tqdm import tqdm

from utils import csvutils as cu
from utils import errorstore as es
from utils import fileutils as fu
from utils import reportutils as ru
from utils import sequenceutils as su

### Parameters ###
# Required arguments
cass_path = ".\\resources\\Cassette_100bp.fa"  # The cassette sequence
ref_path = ".\\resources\\Ref_250bp.fa"  # The sequence for the plasmid into which the cassette has been inserted
test_path = ".\\resources\\Tests.fa"  # The sequencing result

# Optional arguments
extra_nt = 3  # Number of additional nucleotides to be displayed either side of the cleavage site
local_r = 1  # Half width of the local sequences to be extracted at restriction sites
max_gap = 10000 # Maximum number of bp between 3' and 5' restriction sites
min_quality = 1  # Minimum match quality ("1" is perfect)
num_bases = 20  # Number of bases to match
verbose = False  # Display messages during execution

### Processing ###
# Loading reference, cassette and test sequences
filereader = fu.FileReader(verbose=verbose)
ref = filereader.read_sequence(ref_path)[0][0][0]
cass = filereader.read_sequence(cass_path)[0][0][0]
tests = filereader.read_sequence(test_path)[0]

# Creating the PairwiseAligner and SequenceSearcher objects
aligner = Align.PairwiseAligner()
aligner.mode = 'local'
aligner.match_score = 1.0
aligner.gap_score = -1.0
searcher = su.SequenceSearcher(aligner, max_gap=max_gap, min_quality=min_quality, num_bases=num_bases, verbose=verbose)

# Creating the results file
root_name = os.path.splitext(test_path)[0]
file = fu.open_file(root_name, '_individual', 'csv', append_dt=False)
csv_writer = cu.CSVWriter(extra_nt=extra_nt,local_r=local_r,append_dt=False,double_line_mode=False)
file.write(csv_writer._get_individual_header_line())
    
for iteration, test in enumerate(tests):
    count = 0
    header = test[1]
    test = test[0]    
    
    for test_offs in tqdm(range(len(test)), smoothing=0.1):
        test = test[1:]+test[0]
        
        for ref_offs in range(len(ref)):
            ref = ref[1:]+ref[0]
            
            (cleavage_site_t,cleavage_site_b,split) = searcher.get_cleavage_positions(ref, cass, test)
            (local_seq_t, local_seq_b) = su.get_local_sequences(ref,cleavage_site_t,cleavage_site_b,local_r=local_r)

            # Adding this result to the Excel file
            result = (cleavage_site_t, cleavage_site_b, split, local_seq_t, local_seq_b, f"{header}_{test_offs}_{ref_offs}")            
            file.write(csv_writer._get_individual_result_line(result, count, ref))
            
            count = count + 1

file.close()
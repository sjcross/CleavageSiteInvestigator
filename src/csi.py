### IMPORTS ###
import argparse
import os
import sys

from argparse import RawTextHelpFormatter
from Bio import Align
from tqdm import tqdm

from utils import csvutils as cu
from utils import errorstore as es
from utils import fileutils as fu
from utils import plotutils as pu
from utils import reportutils as ru
from utils import sequenceutils as su
from utils import eventmapwriter as emw
from utils import heatmapwritercsv as hmwc
from utils import heatmapwritersvg as hmws

### Parameters ###
### DEFAULT PARAMETERS ###
def_repeat_filter = "" # Expression for filtering sequences by number of repeats

def_extra_nt = 0 # Number of additional nucleotides to be displayed either side of the cleavage site

def_local_r = 1 # Half width of the local sequences to be extracted at restriction sites

def_max_gap = 10000 # Maximum number of bp between 3' and 5' restriction sites

def_min_quality = 1.0 # Minimum match quality ("1" is perfect)

def_num_bases = 20 # Number of bases to match


# HARDCODED PARAMETERS
csv_double_line_mode = True # Output CSV files should use double line format


### ARGUMENT PARSING ###
# Creating ArgumentParser
parser = argparse.ArgumentParser(description= "Cleavage Site Identifier (CSI)\nFor detailed information please visit https://github.com/sjcross/CleavageSiteIdentifier\n\n", add_help=True, formatter_class=RawTextHelpFormatter)

# We want required arguments above optional ones in the help documentation, so removing optional argument descriptions for now
optional = parser._action_groups.pop()

# Defining required arguments
required = parser.add_argument_group('required arguments')

required.add_argument("-c", "--cass_path", type=str, required=True, help="Path to cassette sequence file.  This is the sequence of the cassette inserted at the cleavage site.\n\n")

required.add_argument("-r", "--ref_path", type=str, required=True, help="Path to reference sequence file.  This is the sequence which has been digested.\n\n")

required.add_argument("-t", "--test_path", type=str, required=True, help= "Path to test sequence file.  This is the sequence of the cleaved sample with cassette inserted.\n\n")

# Reinserting optional arguments and defining new values
parser._action_groups.append(optional)

optional.add_argument("-rf", "--repeat_filter", type=str, default=def_repeat_filter, help="Expression defining filter for accepted number of repeats.  Uses standard Python math notation, where 'x' is the number of repeats (e.g. 'x>=3' will process all sequences with at least 3 repeats).\n\n")

optional.add_argument("-en", "--extra_nt", type=int, default=def_extra_nt, help="Number of additional nucleotides to be displayed either side of the cleavage site.\n\n")

optional.add_argument("-lr", "--local_r", type=int, default=def_local_r, help="Half width of the local sequences to be extracted at restriction sites.\n\n")

optional.add_argument("-mg", "--max_gap", type=int, default=def_max_gap, help="Maximum number of bp between 3' and 5' restriction sites.\n\n")

optional.add_argument("-mq", "--min_quality", type=float, default=def_min_quality, help="Minimum match quality (\"1\" is perfect).\n\n")

optional.add_argument("-nb", "--num_bases", type=int, default=def_num_bases, help="Number of bases to match.\n\n")

optional.add_argument("-pr", "--print_results", action='store_true',  help="Prints results in terminal as they are generated.\n\n")

optional.add_argument("-sp", "--show_plots", action='store_true', help="Display plots in pyplot windows as they are generated.\n\n")

optional.add_argument("-we", "--write_eventmap", action='store_true', help="Write event map image to SVG file.  Output file will be stored in test file folder with same name as the test file, but with the suffix '_eventmap'.\n\n")

optional.add_argument("-whsa", "--write_heatmap_svg_auto", action='store_true', help="Write heatmap image (only spanning range of identified event positions)  to SVG file.  Output file will be stored in test file folder with same name as the test file, but with the suffix '_heatmap'.\n\n")

optional.add_argument("-whsf", "--write_heatmap_svg_full", action='store_true', help="Write heatmap image (spanning full range of reference sequence)  to SVG file.  Output file will be stored in test file folder with same name as the test file, but with the suffix '_heatmap'.\n\n")

optional.add_argument("-whca", "--write_heatmap_csv_auto", action='store_true', help="Write heatmap image (only spanning range of identified event positions)  to CSV file.  Output file will be stored in test file folder with same name as the test file, but with the suffix '_heatmap'.\n\n")

optional.add_argument("-whcf", "--write_heatmap_csv_full", action='store_true', help="Write heatmap image (spanning full range of reference sequence)  to CSV file.  Output file will be stored in test file folder with same name as the test file, but with the suffix '_heatmap'.\n\n")

optional.add_argument("-wi", "--write_individual", action='store_true', help="Write individual cleavage results to CSV file.  Output file will be stored in test file folder with same name as the test file, but with the suffix '_individual'.\n\n")

optional.add_argument("-ws", "--write_summary", action='store_true', help="Write summary of results to CSV file.  Output file will be stored in test file folder with same name as the test file, but with the suffix '_summary'.\n\n")

optional.add_argument("-wo", "--write_output", action='store_true', help="Write all content displayed in console to a text file.  Output file will be stored in test file folder with same name as the test file, but with the suffix '_output'.\n\n")

optional.add_argument("-ad", "--append_datetime", action='store_true', help="Append time and date to all output filenames (prevents accidental file overwriting).\n\n")

optional.add_argument("-v", "--verbose", action='store_true', help="Display detailed messages during execution.\n\n")

args = parser.parse_args()

# Required arguments
cass_path = args.cass_path  # The cassette sequence
ref_path = args.ref_path  # The sequence for the plasmid into which the cassette has been inserted
test_path = args.test_path  # The sequencing result

# Optional arguments
repeat_filter = args.repeat_filter # Expression for filtering sequences by number of repeats
extra_nt = args.extra_nt  # Number of additional nucleotides to be displayed either side of the cleavage site
local_r = args.local_r  # Half width of the local sequences to be extracted at restriction sites
max_gap = args.max_gap  # Maximum number of bp between 3' and 5' restriction sites
min_quality = args.min_quality  # Minimum match quality ("1" is perfect)
num_bases = args.num_bases  # Number of bases to match
print_results = args.print_results  # Display results in terminal as they are generated
show_plots = args.show_plots  # Display plots in pyplot windows as they are generated
append_dt = args.append_datetime # Append time and date to all output filenames
write_eventmap = args.write_eventmap # Write event map image to SVG file
write_heatmap_svg_auto = args.write_heatmap_svg_auto # Write heatmap to SVG file for identified event range
write_heatmap_svg_full = args.write_heatmap_svg_full # Write heatmap to SVG file for full reference range
write_heatmap_csv_auto = args.write_heatmap_csv_auto # Write heatmap to CSV file for identified event range
write_heatmap_csv_full = args.write_heatmap_csv_full # Write heatmap to CSV file for full reference range
write_individual = args.write_individual # Write cleavage results to CSV file
write_summary = args.write_summary # Write summary of results to CSV file
write_output = args.write_output # Write console output to text file
verbose = args.verbose  # Display messages during execution

# Getting the root filename
root_name = os.path.splitext(test_path)[0]

# If necessary, redirecting the output stream to file
new_out = None
if write_output:
    new_out = ru.StdOut(root_name)
    sys.stdout = new_out

# If not showing full messages, just display the current file name
if not verbose:
    print("Processing: %s" % test_path)

### Processing ###
# Creating FileHandler object
filereader = fu.FileReader(verbose=verbose)

# Loading reference, cassette and test sequences
if verbose:
    print("INPUT: Loading sequences from file")
ref = filereader.read_sequence(ref_path)[0][0]
cass = filereader.read_sequence(cass_path)[0][0]

(tests,(n_acc,n_rej)) = filereader.read_sequence(test_path,repeat_filter=repeat_filter)
if verbose:
    print("        Accepted = %i (%.2f%%), rejected = %i (%.2f%%)" % (n_acc, (100*n_acc/(n_acc+n_rej)), n_rej, (100*n_rej/(n_acc+n_rej))))

if verbose:
    print("\r")

# Creating the PairwiseAligner and SequenceSearcher objects
aligner = Align.PairwiseAligner()
aligner.mode = 'local'
aligner.match_score = 1.0
aligner.mismatch_score = -1.0
aligner.gap_score = -1.0
searcher = su.SequenceSearcher(aligner, max_gap=max_gap, min_quality=min_quality, num_bases=num_bases, verbose=verbose)

# Dict to store results as dual cleavage site tuple
results = {}
error_store = es.ErrorStore()
error_count = 0

if verbose:
    print("PROCESSING: %i sequence(s)" % len(tests))

for iteration, test in enumerate(tqdm(tests, disable=verbose, smoothing=0.1)):
    if verbose:
        print("    Processing test sequence %i" % (iteration + 1))

    (cleavage_site_t,cleavage_site_b) = searcher.get_cleavage_positions(ref, cass, test, error_store=error_store)
    midpoint_site = searcher.get_midpoint_position(ref, cass, test)

    (local_seq_t, local_seq_b) = su.get_local_sequences(ref,cleavage_site_t,cleavage_site_b,local_r=local_r)

    if cleavage_site_t == None:
        error_count = error_count + 1
        continue

    results[iteration] = (cleavage_site_t, cleavage_site_b, local_seq_t,local_seq_b)

    if verbose:
        print("        Result:")
        ru.print_position(cleavage_site_t, cleavage_site_b, offset="        ")
        ru.print_type(cleavage_site_t, cleavage_site_b, midpoint_site, offset="        ")
        ru.print_sequence(ref,cleavage_site_t,cleavage_site_b, midpoint_site,extra_nt=extra_nt,offset="        ")

# Reporting full sequence frequency
freq_full = ru.get_full_sequence_frequency(results)
freq_local = ru.get_local_sequence_frequency(results, ru.StrandMode.BOTH, ru.LocalMode.BOTH, local_r)
freq_5p = ru.get_local_sequence_frequency(results, ru.StrandMode.BOTH, ru.LocalMode.FIVE_P, local_r)
freq_3p = ru.get_local_sequence_frequency(results, ru.StrandMode.BOTH, ru.LocalMode.THREE_P, local_r)

if print_results:
    print("\rRESULTS:")
    print("    Full sequence frequency:\n")
    ru.print_full_sequence_frequency(ref, freq_full, midpoint_site, extra_nt=extra_nt, offset="")

    print("    Local dinucleotide frequencies:\n")
    ru.print_local_sequence_frequency(freq_local, nonzero_only=False, offset="")

    print("    Local 5' nucleotide frequencies:\n")
    ru.print_local_sequence_frequency(freq_5p, nonzero_only=False, offset="")

    print("    Local 3' nucleotide frequencies:\n")
    ru.print_local_sequence_frequency(freq_3p, nonzero_only=False, offset="")

    # Reporting number of errors
    print("    Summary of errors:\n")
    error_store.print_counts(offset="    ")
    ru.print_error_rate(error_count, len(tests), offset="    ")

# Plotting sequence distributions
if show_plots and len(results) > 0:
    pu.plotFrequency1D(freq_local, freq_5p, freq_3p, show_percentages=True)

    # Reporting top and bottom sequence co-occurrence
    (labels, freq2D) = ru.get_sequence_cooccurrence(results, local_r)
    pu.plotFrequency2D(labels, freq2D, show_percentages=True)

if write_eventmap:
    # Showing cleavage event distribution
    eventmap_writer = emw.EventMapWriter()
    eventmap_writer.write_map(root_name+'_eventmap.svg', freq_full, ref=ref, append_dt=append_dt)

if write_heatmap_svg_auto:
    # Showing events as heatmap
    heatmap_writer = hmws.HeatMapWriterSVG(grid_opts=(False,1,"gray",1), grid_label_opts=(True,12,"gray",100,10), event_label_opts=(False,10,"invert",1,True), sum_show=False)
    heatmap_writer.write_map(root_name+'_autoheatmap.svg', freq_full, None, None, append_dt)

if write_heatmap_svg_full:
    # Showing events as heatmap
    heatmap_writer = hmws.HeatMapWriterSVG(grid_opts=(False,1,"gray",1), grid_label_opts=(True,12,"gray",100,10), event_label_opts=(False,10,"invert",1,True), sum_show=False)
    heatmap_writer.write_map(root_name+'_fullheatmap.svg', freq_full, ref, None, append_dt)

if write_heatmap_csv_auto:
    # Showing events as heatmap
    heatmap_writer = hmwc.HeatMapWriterCSV(sum_show=False)
    heatmap_writer.write_map(root_name+'_autoheatmap.csv', freq_full, None, None, append_dt)

if write_heatmap_csv_full:
    # Showing events as heatmap
    heatmap_writer = hmwc.HeatMapWriterCSV(sum_show=False)
    heatmap_writer.write_map(root_name+'_fullheatmap.csv', freq_full, ref, None, append_dt=append_dt)

# Creating the CSVWriter object
csv_writer = cu.CSVWriter(extra_nt=extra_nt,local_r=local_r,append_dt=append_dt,double_line_mode=csv_double_line_mode)
if write_individual:
    csv_writer.write_individual(root_name, results, ref)

if write_summary:
    csv_writer.write_summary(root_name, freq_full, ref, error_count)

print("\r")

# If writing output, shut down file and print redirection
if write_output:
    new_out.shutdown()

import argparse
import sys

from argparse import RawTextHelpFormatter
from enum import Enum
from utils import histogramwritercsv as hw

# Using an enum rather than bool, so default can be set for argparse
class SHOWHIDE(Enum):
    HIDE = 'hide'
    SHOW = 'show'

    def __str__(self):
        return str(self.value)

### DEFAULT PARAMETER VALUES ###   
def_hist_bin_width = 1
def_command_vis = SHOWHIDE.SHOW

### ARGUMENT PARSING ###
# Creating ArgumentParser
parser = argparse.ArgumentParser(description= "Event type histogram CSV exporter for Cleavage Site Identifier (CSI)\nFor detailed information please visit https://github.com/sjcross/CleavageSiteIdentifier\n\n", add_help=True, formatter_class=RawTextHelpFormatter)

# We want required arguments above optional ones in the help documentation, so removing optional argument descriptions for now
optional = parser._action_groups.pop()

# Defining required arguments
required = parser.add_argument_group('required arguments')

required.add_argument("-d", "--data_path", type=str, required=True, help= "Path to .csv results file.  This file contains details of each sequence to be included in the strand linkage plot.\n\n")

required.add_argument("-o", "--out_path", type=str, help="Path to location where .csv file will be written.  This should be a complete path including the filename and extension.\n\n")

# Reinserting optional arguments and defining new values
parser._action_groups.append(optional)

optional.add_argument("-r", "--ref_path", type=str, help="Path to reference sequence file.  This is the sequence which has been digested.\n\n")

optional.add_argument("-ad", "--append_datetime", action='store_true', help="Append time and date to all output filenames (prevents accidental file overwriting)\n\n")

optional.add_argument("-pr", "--pos_range", type=int, default=[0,0], nargs=2, help="Minimum and maximum strand positions within the reference sequence to display.  Specified as two integer numbers in the order minimum maximum (e.g. -pr 100 200).  If unspecified, the full reference range will be used.\n\n")

optional.add_argument("-h_bw", "--hist_bin_width", type=int, default=def_hist_bin_width, help="Number of positions that will be binned into a single bar on the histogram.  Default: \"%.1f\".\n\n" % def_hist_bin_width)

optional.add_argument("-cmv", "--command_vis", type=SHOWHIDE, default=def_command_vis, choices=list(SHOWHIDE), help="Controls whether the command line command is displayed at the top of the map.  Must be either \"show\" or \"hide\" (e.g. -cmv \"show\").  Default: \"%s\".\n\n" % def_command_vis)

args = parser.parse_args()

pos_range = tuple(args.pos_range) if args.pos_range != [0,0] else None
command_show = args.command_vis is SHOWHIDE.SHOW

writer = hw.HistogramWriterCSV(hist_bin_width=args.hist_bin_width)

commandStr = ' '.join(sys.argv[:]) if command_show else ''
writer.write_from_file(args.data_path, args.out_path, ref_path=args.ref_path, pos_range=pos_range, append_dt=args.append_datetime, commandStr=commandStr)

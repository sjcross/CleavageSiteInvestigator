import argparse

from argparse import RawTextHelpFormatter
from enum import Enum
from utils import heatmapwritercsv as hmw


# Using an enum rather than bool, so default can be set for argparse
class SHOWHIDE(Enum):
    HIDE = 'hide'
    SHOW = 'show'

    def __str__(self):
        return str(self.value)


### DEFAULT PARAMETER VALUES ###   
def_event_label_decimal_places = 1
def_sum_vis = SHOWHIDE.SHOW
def_count_vis = SHOWHIDE.SHOW


### ARGUMENT PARSING ###
# Creating ArgumentParser
parser = argparse.ArgumentParser(description= "Event map CSV exporter for Cleavage Site Identifier (CSI)\nFor detailed information please visit https://github.com/sjcross/CleavageSiteIdentifier\n\n", add_help=True, formatter_class=RawTextHelpFormatter)

# We want required arguments above optional ones in the help documentation, so removing optional argument descriptions for now
optional = parser._action_groups.pop()

# Defining required arguments
required = parser.add_argument_group('required arguments')

required.add_argument("-d", "--data_path", type=str, required=True, help= "Path to .csv results file.  This file contains details of each sequence to be included in the event map.\n\n")

required.add_argument("-o", "--out_path", type=str, help="Path to location where .svg file will be written.  This should be a complete path including the filename and extension.\n\n")

# Reinserting optional arguments and defining new values
parser._action_groups.append(optional)

optional.add_argument("-r", "--ref_path", type=str, help="Path to reference sequence file.  This is the sequence which has been digested.\n\n")

optional.add_argument("-ad", "--append_datetime", action='store_true', help="Append time and date to all output filenames (prevents accidental file overwriting)\n\n")

optional.add_argument("-pr", "--pos_range", type=int, default=[0,0,0,0], nargs=4, help="Minimum and maximum top and bottom strand positions within the reference sequence to display.  Specified as a pair of integer numbers in the order minimum_top maximum_top minimum_bottom maximum_bottom (e.g. -pr 100 200 400 500).  If unspecified, the full reference range will be used.\n\n")

optional.add_argument("-eldp", "--event_label_decimal_places", type=int, default=def_event_label_decimal_places, help="Number of decimal places to use when displaying event frequencies.")

optional.add_argument("-sv", "--sum_vis", type=SHOWHIDE, default=def_sum_vis, choices=list(SHOWHIDE), help="Controls whether the sum row and columns are displayed.  These event cells are rendered using the same settings as for standard events (e.g. font size, colour and zeros visibility).  Must be either \"show\" or \"hide\" (e.g. -sv \"show\").  Default: \"%s\".\n\n" % def_sum_vis)

optional.add_argument("-cv", "--count_vis", type=SHOWHIDE, default=def_count_vis, choices=list(SHOWHIDE), help="Controls whether the total number of events is displayed underneath the map.  Must be either \"show\" or \"hide\" (e.g. -cv \"show\").  Default: \"%s\".\n\n" % def_sum_vis)

args = parser.parse_args()

sum_show = args.sum_vis is SHOWHIDE.SHOW
count_show = args.count_vis is SHOWHIDE.SHOW
pos_range = tuple(args.pos_range) if args.pos_range != [0,0,0,0] else None

writer = hmw.HeatMapWriterCSV(event_label_decimal_places = args.event_label_decimal_places, sum_show=sum_show, count_show=count_show)
writer.write_map_from_file(args.data_path, args.out_path, ref_path=args.ref_path, pos_range=pos_range, append_dt=args.append_datetime)

import argparse

from enum import Enum
from utils import svgutils as svu

from argparse import RawTextHelpFormatter


# Using an enum rather than bool, so default can be set for argparse
class SHOWHIDE(Enum):
    HIDE = 'hide'
    SHOW = 'show'

    def __str__(self):
        return str(self.value)


### DEFAULT PARAMETER VALUES ###
def_im_w = 800
def_im_h = 200

def_dna_rel_top = 0.3
def_dna_rel_bottom = 0.9
def_dna_rel_left = 0.05
def_dna_rel_right = 0.95

def_dna_mode = svu.DNA_MODE.LINE
def_dna_size = 2
def_dna_colour = "black"

def_end_label_vis = SHOWHIDE.SHOW
def_end_label_size = 20
def_end_label_colour = "black"
def_end_label_rel_gap = 0.01

def_grid_vis = SHOWHIDE.SHOW
def_grid_size = 1
def_grid_colour = "lightgray"
def_grid_increment = 100

def_grid_label_vis = SHOWHIDE.SHOW
def_grid_label_size = 12
def_grid_label_colour = "lightgray"
def_grid_label_increment = 500

def_event_max_size = 2
def_event_colour_1 = "blue"
def_event_colour_2 = "red"


### ARGUMENT PARSING ###
# Creating ArgumentParser
parser = argparse.ArgumentParser(description= "Event map SVG renderer for Cleavage Site Identifier (CSI)\nFor detailed information please visit https://github.com/sjcross/CleavageSiteIdentifier", add_help=True, formatter_class=RawTextHelpFormatter)

# We want required arguments above optional ones in the help documentation, so removing optional argument descriptions for now
optional = parser._action_groups.pop()

# Defining required arguments
required = parser.add_argument_group('required arguments')

required.add_argument("-d", "--data_path", type=str, required=True, help= "path to .csv results file.  This file contains details of each sequence to be included in the event map.")

# Reinserting optional arguments and defining new values
parser._action_groups.append(optional)

required.add_argument("-o", "--out_path", type=str, help="")

required.add_argument("-r", "--ref_path", type=str, help="path to reference sequence file.  This is the sequence which has been digested.")

optional.add_argument("-ad", "--append_datetime", action='store_true', help="Append time and date to all output filenames (prevents accidental file overwriting)")

optional.add_argument("-pr", "--pos_range", type=int, default=[0,0], nargs=2, help="")

optional.add_argument("-id", "--im_dims", type=int, default=[def_im_w, def_im_h], nargs=2, help="")

optional.add_argument("-drp", "--dna_rel_pos", type=float, default=[def_dna_rel_top, def_dna_rel_bottom, def_dna_rel_left, def_dna_rel_right], nargs=4, help="")

optional.add_argument("-dm", "--dna_mode", type=svu.DNA_MODE, default=def_dna_mode, choices=list(svu.DNA_MODE), help="")

optional.add_argument("-ds", "--dna_size", type=int, default=def_dna_size, help="")

optional.add_argument("-dc", "--dna_colour", type=str, default=def_dna_colour, help="")

optional.add_argument("-elv", "--end_label_vis", type=SHOWHIDE, default=def_end_label_vis, choices=list(SHOWHIDE), help="")

optional.add_argument("-els", "--end_label_size", type=int, default=def_end_label_size, help="")

optional.add_argument("-elc", "--end_label_colour", type=str, default=def_end_label_colour, help="")

optional.add_argument("-elg", "--end_label_gap", type=float, default=def_end_label_rel_gap, help="")

optional.add_argument("-gv", "--grid_vis", type=SHOWHIDE, default=def_grid_vis, choices=list(SHOWHIDE), help="")

optional.add_argument("-gs", "--grid_size", type=int, default=def_grid_size, help="")

optional.add_argument("-gc", "--grid_colour", type=str, default=def_grid_colour, help="")

optional.add_argument("-gi", "--grid_increment", type=int, default=def_grid_increment, help="")

optional.add_argument("-glv", "--grid_label_vis", type=SHOWHIDE, default=def_grid_label_vis, choices=list(SHOWHIDE), help="")

optional.add_argument("-gls", "--grid_label_size", type=int, default=def_grid_label_size, help="")

optional.add_argument("-glc", "--grid_label_colour", type=str, default=def_grid_label_colour, help="")

optional.add_argument("-gli", "--grid_label_increment", type=int, default=def_grid_label_increment, help="")

optional.add_argument("-ems", "--event_max_size", type=int, default=def_event_max_size, help="")

optional.add_argument("-ec1", "--event_colour_1", type=str, default=def_event_colour_1, help="")

optional.add_argument("-ec2", "--event_colour_2", type=str, default=def_event_colour_2, help="")

args = parser.parse_args()

data_path = args.data_path
ref_path = args.ref_path
end_label_show = True if args.end_label_vis is SHOWHIDE.SHOW else False
grid_show = True if args.grid_vis is SHOWHIDE.SHOW else False
grid_label_show = True if args.grid_label_vis is SHOWHIDE.SHOW else False

# Required arguments
pos_range = tuple(args.pos_range) if args.pos_range != [0,0] else None
im_dims = tuple(args.im_dims)
rel_pos = tuple(args.dna_rel_pos)
dna_opts = (args.dna_mode,args.dna_size,args.dna_colour)
end_label_opts = (end_label_show,args.end_label_size,args.end_label_colour,args.end_label_gap)
grid_opts = (grid_show,args.grid_size,args.grid_colour,args.grid_increment)
grid_label_opts = (grid_label_show,args.grid_label_size,args.grid_label_colour,args.grid_label_increment)
event_opts = (args.event_max_size,args.event_colour_1,args.event_colour_2)

writer = svu.EventMapWriter(im_dims=im_dims, rel_pos=rel_pos, dna_opts=dna_opts, end_label_opts=end_label_opts, grid_opts=grid_opts, grid_label_opts=grid_label_opts, event_opts=event_opts)
writer.write_event_map_from_file(data_path, out_path=args.out_path, ref_path=ref_path, pos_range=pos_range, append_dt=args.append_datetime)

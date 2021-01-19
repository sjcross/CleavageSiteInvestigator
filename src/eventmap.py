import argparse

from argparse import RawTextHelpFormatter
from enum import Enum
from utils import svgutils as svu


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
def_dna_size = 2.5
def_dna_colour = "black"

def_end_label_vis = SHOWHIDE.SHOW
def_end_label_size = 20
def_end_label_colour = "black"
def_end_label_rel_gap = 0.01

def_grid_vis = SHOWHIDE.SHOW
def_grid_size = 1
def_grid_colour = "lightgray"
def_grid_interval = 100

def_grid_label_vis = SHOWHIDE.SHOW
def_grid_label_size = 12
def_grid_label_colour = "lightgray"
def_grid_label_interval = 500

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

required.add_argument("-d", "--data_path", type=str, required=True, help= "Path to .csv results file.  This file contains details of each sequence to be included in the event map.")

required.add_argument("-o", "--out_path", type=str, help="Path to location where .svg file will be written.  This should be a complete path including the filename and extension.")

# Reinserting optional arguments and defining new values
parser._action_groups.append(optional)

required.add_argument("-r", "--ref_path", type=str, help="Path to reference sequence file.  This is the sequence which has been digested.")

optional.add_argument("-ad", "--append_datetime", action='store_true', help="Append time and date to all output filenames (prevents accidental file overwriting)")

optional.add_argument("-pr", "--pos_range", type=int, default=[0,0], nargs=2, help="Minimum and maximum positions within the reference sequence to display.  Specified as a pair of integer numbers in the order minimum maximum (e.g. -pr 100 200).  If unspecified, the full reference range will be used.")

optional.add_argument("-id", "--im_dims", type=int, default=[def_im_w, def_im_h], nargs=2, help="Pixel dimensions of the output .svg image.  Specified as a pair of integer numbers in the order width height (e.g. -id 800 200).  Default: \"%i %i\"." % (def_im_w, def_im_h))

optional.add_argument("-drp", "--dna_rel_pos", type=float, default=[def_dna_rel_top, def_dna_rel_bottom, def_dna_rel_left, def_dna_rel_right], nargs=4, help="Relative position of DNA strands in output image, relative to the top-left corner.  Positions correspond directly to strands, so allowances must be made for labels.  Specified as a list of 4 floating-point numbers in the order top bottom left right (e.g. -drp 0.3 0.9 0.05 0.95).  Default: \"%.2f %.2f %.2f %.2f\"." % (def_dna_rel_top, def_dna_rel_bottom, def_dna_rel_left, def_dna_rel_right))

optional.add_argument("-dm", "--dna_mode", type=svu.DNA_MODE, default=def_dna_mode, choices=list(svu.DNA_MODE), help="Rendering mode for DNA strands.  \"line\" displays DNA as a pair of solid lines with width controlled by --dna_size.  \"seq\" displays DNA as a pair of letter sequences with font size controlled by --dna_size (Note: position range must be sufficiently narrow in sequence mode to prevent text overlap).  \"none\" displays nothing.  Default: \"%s\"." % def_dna_mode)

optional.add_argument("-ds", "--dna_size", type=int, default=def_dna_size, help="Width of rendered DNA lines (when in --dna_mode \"line\") or font size of DNA sequence letters (when in --dna_mode \"seq\").  Default: \"%.1f\"." % def_dna_size)

optional.add_argument("-dc", "--dna_colour", type=str, default=def_dna_colour, help="Colour of rendered DNA lines or text sequences.  Can be specified as colour names (e.g. \"black\"), as hex values (e.g. \"#16C3D6\" for a light blue) or as rgb values in the range 0-255 (e.g. \"rgb(128,0,128)\" for purple).  Default: \"%s\"." % def_dna_colour)

optional.add_argument("-elv", "--end_label_vis", type=SHOWHIDE, default=def_end_label_vis, choices=list(SHOWHIDE), help="Controls whether the DNA end labels (\"5'\" and \"3'\") are rendered.  Must be either \"show\" or \"hide\" (e.g. -elv \"show\").  Default: \"%s\"." % def_end_label_vis)

optional.add_argument("-els", "--end_label_size", type=int, default=def_end_label_size, help="Font size of DNA end labels.  Default: \"%.1f\"." % def_end_label_size)

optional.add_argument("-elc", "--end_label_colour", type=str, default=def_end_label_colour, help="Colour of the rendered DNA end labels.  Can be specified as colour names (e.g. \"black\"), as hex values (e.g. \"#16C3D6\" for a light blue) or as rgb values in the range 0-255 (e.g. \"rgb(128,0,128)\" for purple).  Default: \"%s\"." % def_end_label_colour)

optional.add_argument("-elg", "--end_label_gap", type=float, default=def_end_label_rel_gap, help="Relative gap between end labels and the rendered DNA strands.  Specified as a proportion of the total image width.  Default: \"%.2f\"." % def_end_label_rel_gap)

optional.add_argument("-gv", "--grid_vis", type=SHOWHIDE, default=def_grid_vis, choices=list(SHOWHIDE), help="Controls whether grid lines are rendered.  Must be either \"show\" or \"hide\" (e.g. -gv \"show\").  Default: \"%s\"." % def_grid_vis)

optional.add_argument("-gs", "--grid_size", type=int, default=def_grid_size, help="Width of grid lines.  Default: \"%.1f\"." % def_grid_size)

optional.add_argument("-gc", "--grid_colour", type=str, default=def_grid_colour, help="Colour of the rendered grid lines.  Can be specified as colour names (e.g. \"black\"), as hex values (e.g. \"#16C3D6\" for a light blue) or as rgb values in the range 0-255 (e.g. \"rgb(128,0,128)\" for purple).  Default: \"%s\"." % def_grid_colour)

optional.add_argument("-gi", "--grid_interval", type=int, default=def_grid_interval, help="Interval between adjacent grid lines.  Default: \"%i\"." % def_grid_interval)

optional.add_argument("-glv", "--grid_label_vis", type=SHOWHIDE, default=def_grid_label_vis, choices=list(SHOWHIDE), help="Controls whether grid labels are rendered.  Grid labels are always shown above the top strand and are rotated vertically.  Must be either \"show\" or \"hide\" (e.g. -glv \"show\").  Default: \"%s\"." % def_grid_label_vis)

optional.add_argument("-gls", "--grid_label_size", type=int, default=def_grid_label_size, help="Font size of grid labels.  Default: \"%.1f\"." % def_grid_label_size)

optional.add_argument("-glc", "--grid_label_colour", type=str, default=def_grid_label_colour, help="Colour of the rendered grid labels.  Can be specified as colour names (e.g. \"black\"), as hex values (e.g. \"#16C3D6\" for a light blue) or as rgb values in the range 0-255 (e.g. \"rgb(128,0,128)\" for purple).  Default: \"%s\"." % def_grid_label_colour)

optional.add_argument("-gli", "--grid_label_interval", type=int, default=def_grid_label_interval, help="Interval between adjacent grid labels.  Default: \"%i\"." % def_grid_label_interval)


optional.add_argument("-ems", "--event_max_size", type=int, default=def_event_max_size, help="Width of line corresponding to the most frequent cleavage event.  Default: \"%.1f\"." % def_event_max_size)

optional.add_argument("-ec1", "--event_colour_1", type=str, default=def_event_colour_1, help="Colour corresponding to least frequent cleavage event.  Event line colours will be interpolated between this colour and --event_colour_2.  Can be specified as colour names (e.g. \"black\"), as hex values (e.g. \"#16C3D6\" for a light blue) or as rgb values in the range 0-255 (e.g. \"rgb(128,0,128)\" for purple).  Default: \"%s\"." % def_event_colour_1)

optional.add_argument("-ec2", "--event_colour_2", type=str, default=def_event_colour_2, help="Colour corresponding to most frequent cleavage event.  Event line colours will be interpolated between --event_colour_1 and this colour.  Can be specified as colour names (e.g. \"black\"), as hex values (e.g. \"#16C3D6\" for a light blue) or as rgb values in the range 0-255 (e.g. \"rgb(128,0,128)\" for purple).  Default: \"%s\"." % def_event_colour_2)

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
grid_opts = (grid_show,args.grid_size,args.grid_colour,args.grid_interval)
grid_label_opts = (grid_label_show,args.grid_label_size,args.grid_label_colour,args.grid_label_interval)
event_opts = (args.event_max_size,args.event_colour_1,args.event_colour_2)



writer = svu.EventMapWriter(im_dims=im_dims, rel_pos=rel_pos, dna_opts=dna_opts, end_label_opts=end_label_opts, grid_opts=grid_opts, grid_label_opts=grid_label_opts, event_opts=event_opts)
writer.write_event_map_from_file(data_path, args.out_path, ref_path=ref_path, pos_range=pos_range, append_dt=args.append_datetime)

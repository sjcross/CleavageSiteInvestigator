import argparse
from matplotlib import pyplot as plt

from argparse import RawTextHelpFormatter
from enum import Enum
from utils import heatmapwritersvg as hmw


# Using an enum rather than bool, so default can be set for argparse
class SHOWHIDE(Enum):
    HIDE = 'hide'
    SHOW = 'show'

    def __str__(self):
        return str(self.value)


### DEFAULT PARAMETER VALUES ###
def_im_dim = 800

def_map_rel_top = 0.1
def_map_rel_left = 0.1
def_map_rel_size = 0.8

def_border_vis = SHOWHIDE.SHOW
def_border_size = 2.5
def_border_colour = "black"

def_axis_label_vis = SHOWHIDE.SHOW
def_axis_label_size = 16
def_axis_label_colour = "gray"
def_axis_label_gap = 50

def_grid_vis = SHOWHIDE.SHOW
def_grid_size = 1
def_grid_colour = "gray"
def_grid_interval = 1

def_grid_label_vis = SHOWHIDE.SHOW
def_grid_label_size = 12
def_grid_label_colour = "gray"
def_grid_label_interval = 100
def_grid_label_gap = 10

def_event_colourmap = "plasma"

def_event_label_vis = SHOWHIDE.SHOW
def_event_label_size = 12
def_event_label_colour = "invert"
def_event_label_decimal_places = 1
def_event_label_zeros_vis = SHOWHIDE.SHOW
    
def_sum_vis = SHOWHIDE.SHOW


### ARGUMENT PARSING ###
# Creating ArgumentParser
parser = argparse.ArgumentParser(description= "Event map SVG renderer for Cleavage Site Identifier (CSI)\nFor detailed information please visit https://github.com/sjcross/CleavageSiteIdentifier\n\n", add_help=True, formatter_class=RawTextHelpFormatter)

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

optional.add_argument("-id", "--im_dim", type=int, default=def_im_dim, help="Pixel dimensions of the output .svg image.  Event map image is square, so specified as a single integer number.  Default: \"%i\".\n\n" % def_im_dim)

optional.add_argument("-mrp", "--map_rel_pos", type=float, default=[def_map_rel_top, def_map_rel_left, def_map_rel_size], nargs=3, help="Position and size of event map in output image, relative to the top-left corner.  Positions correspond directly to map, so allowances must be made for labels.  Specified as a list of 3 floating-point numbers in the order top height size (e.g. -drp 0.3 0.3 0.6).  Default: \"%.2f %.2f %.2f\".\n\n" % (def_map_rel_top, def_map_rel_left, def_map_rel_size))

optional.add_argument("-bv", "--border_vis", type=SHOWHIDE, default=def_border_vis, choices=list(SHOWHIDE), help="Controls whether the border of the map is rendered as a line.  Must be either \"show\" or \"hide\" (e.g. -bv \"show\").  Default: \"%s\".\n\n" % def_border_vis)

optional.add_argument("-bs", "--border_size", type=int, default=def_border_size, help="Width of rendered map border lines.  Default: \"%.1f\".\n\n" % def_border_size)

optional.add_argument("-bc", "--border_colour", type=str, default=def_border_colour, help="Colour of rendered map border lines.  Can be specified as colour names (e.g. \"black\"), as hex values (e.g. \"#16C3D6\" for a light blue) or as rgb values in the range 0-255 (e.g. \"rgb(128,0,128)\" for purple).  Default: \"%s\".\n\n" % def_border_colour)

optional.add_argument("-alv", "--axis_label_vis", type=SHOWHIDE, default=def_axis_label_vis, choices=list(SHOWHIDE), help="Controls whether the axis labels (\"Top strand\" and \"Bottom strand\") are rendered.  Must be either \"show\" or \"hide\" (e.g. -alv \"show\").  Default: \"%s\".\n\n" % def_axis_label_vis)

optional.add_argument("-als", "--axis_label_size", type=int, default=def_axis_label_size, help="Font size of axis labels.  Default: \"%.1f\".\n\n" % def_axis_label_size)

optional.add_argument("-alc", "--axis_label_colour", type=str, default=def_axis_label_colour, help="Colour of the rendered axis labels.  Can be specified as colour names (e.g. \"black\"), as hex values (e.g. \"#16C3D6\" for a light blue) or as rgb values in the range 0-255 (e.g. \"rgb(128,0,128)\" for purple).  Default: \"%s\".\n\n" % def_axis_label_colour)

optional.add_argument("-alg", "--axis_label_gap", type=float, default=def_axis_label_gap, help="Gap between axis labels and the top/left edges of the map region (allowances must be made for grid labels if enabled).  Specified in integer pixel units.  Default: \"%i\".\n\n" % def_axis_label_gap)

optional.add_argument("-gv", "--grid_vis", type=SHOWHIDE, default=def_grid_vis, choices=list(SHOWHIDE), help="Controls whether grid lines are rendered.  Must be either \"show\" or \"hide\" (e.g. -gv \"show\").  Default: \"%s\".\n\n" % def_grid_vis)

optional.add_argument("-gs", "--grid_size", type=int, default=def_grid_size, help="Width of grid lines.  Default: \"%.1f\".\n\n" % def_grid_size)

optional.add_argument("-gc", "--grid_colour", type=str, default=def_grid_colour, help="Colour of the rendered grid lines.  Can be specified as colour names (e.g. \"black\"), as hex values (e.g. \"#16C3D6\" for a light blue) or as rgb values in the range 0-255 (e.g. \"rgb(128,0,128)\" for purple).  Default: \"%s\".\n\n" % def_grid_colour)

optional.add_argument("-gi", "--grid_interval", type=int, default=def_grid_interval, help="Interval between adjacent grid lines.  Default: \"%i\".\n\n" % def_grid_interval)

optional.add_argument("-glv", "--grid_label_vis", type=SHOWHIDE, default=def_grid_label_vis, choices=list(SHOWHIDE), help="Controls whether grid labels are rendered.  Grid labels are always shown above the top strand and are rotated vertically.  Must be either \"show\" or \"hide\" (e.g. -glv \"show\").  Default: \"%s\".\n\n" % def_grid_label_vis)

optional.add_argument("-gls", "--grid_label_size", type=int, default=def_grid_label_size, help="Font size of grid labels.  Default: \"%.1f\".\n\n" % def_grid_label_size)

optional.add_argument("-glc", "--grid_label_colour", type=str, default=def_grid_label_colour, help="Colour of the rendered grid labels.  Can be specified as colour names (e.g. \"black\"), as hex values (e.g. \"#16C3D6\" for a light blue) or as rgb values in the range 0-255 (e.g. \"rgb(128,0,128)\" for purple).  Default: \"%s\".\n\n" % def_grid_label_colour)

optional.add_argument("-gli", "--grid_label_interval", type=int, default=def_grid_label_interval, help="Interval between adjacent grid labels.  Default: \"%i\".\n\n" % def_grid_label_interval)

optional.add_argument("-glg", "--grid_label_gap", type=float, default=def_grid_label_gap, help="Gap between grid labels and the rendered DNA strands.  Specified in integer pixel units.  Default: \"%i\".\n\n" % def_grid_label_gap)

optional.add_argument("-ec", "--event_colourmap", type=str, default=def_event_colourmap, help="Matplotlib colourmap to use for event plotting.  Must be a Matplotlib colourmap.  Default: \"%s\".\n\n" % def_event_colourmap)

optional.add_argument("-elv", "--event_label_vis", type=SHOWHIDE, default=def_event_label_vis, choices=list(SHOWHIDE), help="Controls whether the event labels (showing percentage of total events) are rendered.  Must be either \"show\" or \"hide\" (e.g. -elv \"show\").  Default: \"%s\".\n\n" % def_event_label_vis)

optional.add_argument("-els", "--event_label_size", type=int, default=def_event_label_size, help="Font size of event labels.  Default: \"%.1f\".\n\n" % def_event_label_size)

optional.add_argument("-elc", "--event_label_colour", type=str, default=def_event_label_colour, help="Colour of the rendered event labels.  Can be specified as colour names (e.g. \"black\"), as hex values (e.g. \"#16C3D6\" for a light blue) or as rgb values in the range 0-255 (e.g. \"rgb(128,0,128)\" for purple).  Event labels have the special colour option \"invert\", which uses the opposite RGB colour to the corresponding event square to improve label contrast.  Default: \"%s\".\n\n" % def_event_label_colour)

optional.add_argument("-eldp", "--event_label_decimal_places", type=int, default=def_event_label_decimal_places, help="Number of decimal places to use when displaying event frequencies.")

optional.add_argument("-elzv", "--event_label_zeros_vis", type=SHOWHIDE, default=def_event_label_zeros_vis, choices=list(SHOWHIDE), help="Controls whether labels for positions with zero events are rendered.  Showing these labels can greatly increase file size - especially when rendering a large position range.  Must be either \"show\" or \"hide\" (e.g. -elzv \"show\").  Default: \"%s\".\n\n" % def_event_label_zeros_vis)

optional.add_argument("-sv", "--sum_vis", type=SHOWHIDE, default=def_sum_vis, choices=list(SHOWHIDE), help="Controls whether the sum row and columns are rendered.  These event cells are rendered using the same settings as for standard events (e.g. font size, colour and zeros visibility).  Must be either \"show\" or \"hide\" (e.g. -elv \"show\").  Default: \"%s\".\n\n" % def_sum_vis)

args = parser.parse_args()

data_path = args.data_path
ref_path = args.ref_path
border_show = True if args.border_vis is SHOWHIDE.SHOW else False
axis_label_show = True if args.axis_label_vis is SHOWHIDE.SHOW else False
grid_show = True if args.grid_vis is SHOWHIDE.SHOW else False
grid_label_show = True if args.grid_label_vis is SHOWHIDE.SHOW else False
event_label_show = True if args.event_label_vis is SHOWHIDE.SHOW else False

# Required arguments
pos_range = tuple(args.pos_range) if args.pos_range != [0,0,0,0] else None
im_dim = args.im_dim
rel_pos = tuple(args.map_rel_pos)
border_opts = (border_show,args.border_size,args.border_colour)
axis_label_opts = (axis_label_show,args.axis_label_size,args.axis_label_colour,args.axis_label_gap)
grid_opts = (grid_show,args.grid_size,args.grid_colour,args.grid_interval)
grid_label_opts = (grid_label_show,args.grid_label_size,args.grid_label_colour,args.grid_label_interval,args.grid_label_gap)
event_colourmap = args.event_colourmap
event_label_opts = (event_label_show,args.event_label_size,args.event_label_colour,args.event_label_decimal_places,args.event_label_zeros_vis)

writer = hmw.HeatMapWriterSVG(im_dim=im_dim, rel_pos=rel_pos, border_opts=border_opts, grid_opts=grid_opts, grid_label_opts=grid_label_opts, event_colourmap=event_colourmap)
writer.write_map_from_file(data_path, args.out_path, ref_path=ref_path, pos_range=pos_range, append_dt=args.append_datetime)

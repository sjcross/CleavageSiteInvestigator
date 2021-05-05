import argparse

from argparse import RawTextHelpFormatter
from enum import Enum
from utils import eventmapwriter as emw


# Using an enum rather than bool, so default can be set for argparse
class SHOWHIDE(Enum):
    HIDE = 'hide'
    SHOW = 'show'

    def __str__(self):
        return str(self.value)


### DEFAULT PARAMETER VALUES ###
def_im_w = 800
def_im_h = 300

def_map_rel_top = 0.4
def_map_rel_height = 0.3
def_map_rel_left = 0.05
def_map_rel_width = 0.8

def_dna_mode = emw.DNA_MODE.LINE
def_dna_size = 2.5
def_dna_colour = "black"

def_end_label_vis = SHOWHIDE.SHOW
def_end_label_size = 20
def_end_label_colour = "black"
def_end_label_rel_gap = 0.01

def_grid_vis = SHOWHIDE.SHOW
def_grid_size = 1
def_grid_colour = "gray"
def_grid_interval = 100

def_grid_label_vis = SHOWHIDE.SHOW
def_grid_label_size = 12
def_grid_label_colour = "gray"
def_grid_label_interval = 500
def_grid_label_rel_gap = 0.02

def_colourbar_vis = SHOWHIDE.SHOW
def_colourbar_rel_left = 0.91
def_colourbar_rel_width = 0.02
def_colourbar_border_size = 1

def_colourbar_label_vis = SHOWHIDE.SHOW
def_colourbar_label_size = 12
def_colourbar_label_colour = "black"
def_colourbar_label_interval = 50

def_event_min_size = 0.5
def_event_max_size = 2
def_event_colourmap = "cool"
def_event_min_range = 0
def_event_max_range = 100
def_event_opacity = 0.2
def_event_stack_order = 1


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

optional.add_argument("-pr", "--pos_range", type=int, default=[0,0], nargs=2, help="Minimum and maximum positions within the reference sequence to display.  Specified as a pair of integer numbers in the order minimum maximum (e.g. -pr 100 200).  If unspecified, the full reference range will be used.\n\n")

optional.add_argument("-id", "--im_dims", type=int, default=[def_im_w, def_im_h], nargs=2, help="Pixel dimensions of the output .svg image.  Specified as a pair of integer numbers in the order width height (e.g. -id 800 200).  Default: \"%i %i\".\n\n" % (def_im_w, def_im_h))

optional.add_argument("-mrp", "--map_rel_pos", type=float, default=[def_map_rel_top, def_map_rel_height, def_map_rel_left, def_map_rel_width], nargs=4, help="Position and size of event map in output image, relative to the top-left corner.  Positions correspond directly to map, so allowances must be made for labels.  Specified as a list of 4 floating-point numbers in the order top height left width (e.g. -mrp 0.3 0.6 0.05 0.9).  Default: \"%.2f %.2f %.2f %.2f\".\n\n" % (def_map_rel_top, def_map_rel_height, def_map_rel_left, def_map_rel_width))

optional.add_argument("-dm", "--dna_mode", type=emw.DNA_MODE, default=def_dna_mode, choices=list(emw.DNA_MODE), help="Rendering mode for DNA strands.  \"line\" displays DNA as a pair of solid lines with width controlled by --dna_size.  \"seq\" displays DNA as a pair of letter sequences with font size controlled by --dna_size (Note: position range must be sufficiently narrow in sequence mode to prevent text overlap).  \"none\" displays nothing.  Default: \"%s\".\n\n" % def_dna_mode)

optional.add_argument("-ds", "--dna_size", type=int, default=def_dna_size, help="Width of rendered DNA lines (when in --dna_mode \"line\") or font size of DNA sequence letters (when in --dna_mode \"seq\").  Default: \"%.1f\".\n\n" % def_dna_size)

optional.add_argument("-dc", "--dna_colour", type=str, default=def_dna_colour, help="Colour of rendered DNA lines or text sequences.  Can be specified as colour names (e.g. \"black\"), as hex values (e.g. \"#16C3D6\" for a light blue) or as rgb values in the range 0-255 (e.g. \"rgb(128,0,128)\" for purple).  Default: \"%s\".\n\n" % def_dna_colour)

optional.add_argument("-elv", "--end_label_vis", type=SHOWHIDE, default=def_end_label_vis, choices=list(SHOWHIDE), help="Controls whether the DNA end labels (\"5'\" and \"3'\") are rendered.  Must be either \"show\" or \"hide\" (e.g. -elv \"show\").  Default: \"%s\".\n\n" % def_end_label_vis)

optional.add_argument("-els", "--end_label_size", type=int, default=def_end_label_size, help="Font size of DNA end labels.  Default: \"%.1f\".\n\n" % def_end_label_size)

optional.add_argument("-elc", "--end_label_colour", type=str, default=def_end_label_colour, help="Colour of the rendered DNA end labels.  Can be specified as colour names (e.g. \"black\"), as hex values (e.g. \"#16C3D6\" for a light blue) or as rgb values in the range 0-255 (e.g. \"rgb(128,0,128)\" for purple).  Default: \"%s\".\n\n" % def_end_label_colour)

optional.add_argument("-elrg", "--end_label_rel_gap", type=float, default=def_end_label_rel_gap, help="Gap between end labels and the rendered DNA strands.  Specified as a fraction of the image width.  Default: \"%i\".\n\n" % def_end_label_rel_gap)

optional.add_argument("-gv", "--grid_vis", type=SHOWHIDE, default=def_grid_vis, choices=list(SHOWHIDE), help="Controls whether grid lines are rendered.  Must be either \"show\" or \"hide\" (e.g. -gv \"show\").  Default: \"%s\".\n\n" % def_grid_vis)

optional.add_argument("-gs", "--grid_size", type=int, default=def_grid_size, help="Width of grid lines.  Default: \"%.1f\".\n\n" % def_grid_size)

optional.add_argument("-gc", "--grid_colour", type=str, default=def_grid_colour, help="Colour of the rendered grid lines.  Can be specified as colour names (e.g. \"black\"), as hex values (e.g. \"#16C3D6\" for a light blue) or as rgb values in the range 0-255 (e.g. \"rgb(128,0,128)\" for purple).  Default: \"%s\".\n\n" % def_grid_colour)

optional.add_argument("-gi", "--grid_interval", type=int, default=def_grid_interval, help="Interval between adjacent grid lines.  Default: \"%i\".\n\n" % def_grid_interval)

optional.add_argument("-glv", "--grid_label_vis", type=SHOWHIDE, default=def_grid_label_vis, choices=list(SHOWHIDE), help="Controls whether grid labels are rendered.  Grid labels are always shown above the top strand and are rotated vertically.  Must be either \"show\" or \"hide\" (e.g. -glv \"show\").  Default: \"%s\".\n\n" % def_grid_label_vis)

optional.add_argument("-gls", "--grid_label_size", type=int, default=def_grid_label_size, help="Font size of grid labels.  Default: \"%.1f\".\n\n" % def_grid_label_size)

optional.add_argument("-glc", "--grid_label_colour", type=str, default=def_grid_label_colour, help="Colour of the rendered grid labels.  Can be specified as colour names (e.g. \"black\"), as hex values (e.g. \"#16C3D6\" for a light blue) or as rgb values in the range 0-255 (e.g. \"rgb(128,0,128)\" for purple).  Default: \"%s\".\n\n" % def_grid_label_colour)

optional.add_argument("-gli", "--grid_label_interval", type=int, default=def_grid_label_interval, help="Interval between adjacent grid labels.  Default: \"%i\".\n\n" % def_grid_label_interval)

optional.add_argument("-glrg", "--grid_label_rel_gap", type=float, default=def_grid_label_rel_gap, help="Gap between grid labels and the rendered DNA strands.  Specified as a fraction of the image height.  Default: \"%i\".\n\n" % def_grid_label_rel_gap)

optional.add_argument("-cv", "--colourbar_vis", type=SHOWHIDE, default=def_colourbar_vis, choices=list(SHOWHIDE), help="Controls whether the colourbar is rendered.  Must be either \"show\" or \"hide\" (e.g. -gv \"show\").  Default: \"%s\".\n\n" % def_colourbar_vis)

optional.add_argument("-crp", "--colourbar_rel_pos", type=float, default=[def_colourbar_rel_left,def_colourbar_rel_width], nargs=2, help="Position and size of colourbar in output image, relative to the top-left corner.  Specified as a list of 2 floating-point numbers in the order left width (e.g. -drp 0.9 0.02).  Default: \"%.2f %.2f\".\n\n" % (def_colourbar_rel_left, def_colourbar_rel_width))

optional.add_argument("-cs", "--colourbar_size", type=str, default=def_colourbar_border_size, help="Width of colourbar border.  Default: \"%s\".\n\n" % def_colourbar_border_size)

optional.add_argument("-clv", "--colourbar_label_vis", type=SHOWHIDE, default=def_colourbar_label_vis, choices=list(SHOWHIDE), help="Controls whether colourbar labels are rendered.  Must be either \"show\" or \"hide\" (e.g. -glv \"show\").  Default: \"%s\".\n\n" % def_colourbar_label_vis)

optional.add_argument("-cls", "--colourbar_label_size", type=int, default=def_colourbar_label_size, help="Font size of colourbar labels.  Default: \"%.1f\".\n\n" % def_colourbar_label_size)

optional.add_argument("-clc", "--colourbar_label_colour", type=str, default=def_colourbar_label_colour, help="Colour of the rendered colourbar labels.  Can be specified as colour names (e.g. \"black\"), as hex values (e.g. \"#16C3D6\" for a light blue) or as rgb values in the range 0-255 (e.g. \"rgb(128,0,128)\" for purple).  Default: \"%s\".\n\n" % def_colourbar_label_colour)

optional.add_argument("-cli", "--colourbar_label_interval", type=int, default=def_colourbar_label_interval, help="Interval between labels shown in colourbar.  Default: \"%.1f\".\n\n" % def_colourbar_label_interval)

optional.add_argument("-emis", "--event_min_size", type=float, default=def_event_min_size, help="Width of line corresponding to the least frequent cleavage event.  Default: \"%.1f\".\n\n" % def_event_min_size)

optional.add_argument("-emas", "--event_max_size", type=float, default=def_event_max_size, help="Width of line corresponding to the most frequent cleavage event.  Default: \"%.1f\".\n\n" % def_event_max_size)

optional.add_argument("-ec", "--event_colourmap", type=str, default=def_event_colourmap, help="Matplotlib colourmap to use for event plotting.  Must be a Matplotlib colourmap.  Default: \"%s\".\n\n" % def_event_colourmap)

optional.add_argument("-eo", "--event_opacity", type=float, default=def_event_opacity, help="Opacity of each event line.  Values in the range 0-1.  Default: \"%.1f\".\n\n" % def_event_opacity)

optional.add_argument("-er", "--event_range", type=int, default=[def_event_min_range,def_event_max_range], nargs=2, help="Range of values colourscale will span (specified as percentage of all events).  For automatic range selection, set both values to -1 (e.g. -er -1 -1).  Default: \"%i %i\".\n\n" % (def_event_min_range,def_event_max_range))

optional.add_argument("-eso", "--event_stack_order", type=int, default=def_event_stack_order, help="Mode for ordering events.  Options: 1 (most frequent at back), 2 (most frequent at front).  Default: \"%.1f\".\n\n" % def_event_stack_order)

args = parser.parse_args()

end_label_show = args.end_label_vis is SHOWHIDE.SHOW
grid_show = args.grid_vis is SHOWHIDE.SHOW
grid_label_show = args.grid_label_vis is SHOWHIDE.SHOW
colourbar_show = args.colourbar_vis is SHOWHIDE.SHOW
colourbar_label_show = args.colourbar_label_vis is SHOWHIDE.SHOW

# Required arguments
pos_range = tuple(args.pos_range) if args.pos_range != [0,0] else None
im_dims = tuple(args.im_dims)
rel_pos = tuple(args.map_rel_pos)
dna_opts = (args.dna_mode,args.dna_size,args.dna_colour)
end_label_opts = (end_label_show,args.end_label_size,args.end_label_colour,args.end_label_rel_gap)
grid_opts = (grid_show,args.grid_size,args.grid_colour,args.grid_interval)
grid_label_opts = (grid_label_show,args.grid_label_size,args.grid_label_colour,args.grid_label_interval,args.grid_label_rel_gap)
colourbar_opts = (colourbar_show,args.colourbar_rel_pos[0],args.colourbar_rel_pos[1],args.colourbar_size)
colourbar_label_opts = (colourbar_label_show,args.colourbar_label_size,args.colourbar_label_colour,args.colourbar_label_interval)
event_opts = (args.event_min_size,args.event_max_size,args.event_colourmap,args.event_range[0],args.event_range[1],args.event_opacity,args.event_stack_order)

writer = emw.EventMapWriter(im_dims=im_dims, rel_pos=rel_pos, dna_opts=dna_opts, end_label_opts=end_label_opts, grid_opts=grid_opts, grid_label_opts=grid_label_opts, colourbar_opts=colourbar_opts, colourbar_label_opts=colourbar_label_opts, event_opts=event_opts)
writer.write_map_from_file(args.data_path, args.out_path, ref_path=args.ref_path, pos_range=pos_range, append_dt=args.append_datetime)

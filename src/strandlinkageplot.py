import argparse

from argparse import RawTextHelpFormatter
from enum import Enum
from utils import strandlinkageplotwriter as slpw


# Using an enum rather than bool, so default can be set for argparse
class SHOWHIDE(Enum):
    HIDE = 'hide'
    SHOW = 'show'

    def __str__(self):
        return str(self.value)

class HISTSPLIT(Enum):
    COMBI = 'combi'
    SPLIT = 'split'

    def __str__(self):
        return str(self.value)

### DEFAULT PARAMETER VALUES ###
def_im_w = 800
def_im_h = 700

def_font = "Arial"

def_map_rel_top = 0.43
def_map_rel_height = 0.14
def_map_rel_left = 0.075
def_map_rel_width = 0.78

def_dna_mode = slpw.DNA_MODE.LINE
def_dna_size = 2.5
def_dna_colour = "black"
def_dna_rel_gap = 0.01

def_end_label_vis = SHOWHIDE.SHOW
def_end_label_size = 20
def_end_label_colour = "black"
def_end_label_rel_gap = 0.01
def_end_label_position = slpw.VPOS.CENTRE

def_grid_vis = SHOWHIDE.SHOW
def_grid_size = 1
def_grid_colour = "lightgray"
def_grid_interval = 100

def_grid_label_vis = SHOWHIDE.SHOW
def_grid_label_size = 12
def_grid_label_colour = "gray"
def_grid_label_interval = 500
def_grid_label_rel_gap = 0.015

def_cbar_vis = SHOWHIDE.SHOW
def_cbar_rel_left = 0.91
def_cbar_rel_width = 0.02
def_cbar_border_size = 1

def_cbar_label_vis = SHOWHIDE.SHOW
def_cbar_label_size = 12
def_cbar_label_colour = "gray"
def_cbar_label_interval = 25
def_cbar_label_rel_gap = 0.02

def_event_min_size = 0.5
def_event_max_size = 2
def_event_colourmap = "cool"
def_event_min_range = 0
def_event_max_range = 100
def_event_outside_range_vis = SHOWHIDE.HIDE
def_event_opacity = 1
def_event_stack_order = 1

def_hist_vis = SHOWHIDE.SHOW
def_hist_min_range = 0
def_hist_max_range = 100
def_hist_bin_width = 1
def_hist_colour = "darkgray"
def_hist_rel_height = 0.06
def_hist_rel_gap = 0.03
def_hist_pc_bar_gap = 0
def_hist_overhang = 0
def_hist_log_vis = SHOWHIDE.HIDE

def_hist_label_vis = SHOWHIDE.SHOW
def_hist_label_size = 12
def_hist_label_colour = "gray"
def_hist_label_interval = 50
def_hist_label_rel_gap = 0.01
def_hist_label_position = slpw.HPOS.LEFT
def_hist_label_zero_vis = SHOWHIDE.SHOW

def_hist_grid_vis = SHOWHIDE.SHOW
def_hist_grid_size = 1
def_hist_grid_colour = "lightgray"
def_hist_grid_interval = 25

def_hist_name_vis = SHOWHIDE.SHOW
def_hist_name_size = 12
def_hist_name_colour = "gray"
def_hist_name_rel_gap = 0.01
def_hist_name_position = slpw.HPOS.RIGHT

def_splithist_vis = SHOWHIDE.SHOW
def_splithist_min_range = 0
def_splithist_max_range = 100
def_splithist_colour = "darkgray"
def_splithist_rel_height = 0.06

def_splithist_label_interval = 50

def_splithist_grid_interval = 25

### ARGUMENT PARSING ###
# Creating ArgumentParser
parser = argparse.ArgumentParser(description= "Strand linkage plot SVG renderer for Cleavage Site Identifier (CSI)\nFor detailed information please visit https://github.com/sjcross/CleavageSiteIdentifier\n\n", add_help=True, formatter_class=RawTextHelpFormatter)

# We want required arguments above optional ones in the help documentation, so removing optional argument descriptions for now
optional = parser._action_groups.pop()

# Defining required arguments
required = parser.add_argument_group('required arguments')

required.add_argument("-d", "--data_path", type=str, required=True, help= "Path to .csv results file.  This file contains details of each sequence to be included in the strand linkage plot.\n\n")

required.add_argument("-o", "--out_path", type=str, help="Path to location where .svg file will be written.  This should be a complete path including the filename and extension.\n\n")

# Reinserting optional arguments and defining new values
parser._action_groups.append(optional)

optional.add_argument("-r", "--ref_path", type=str, help="Path to reference sequence file.  This is the sequence which has been digested.\n\n")

optional.add_argument("-ad", "--append_datetime", action='store_true', help="Append time and date to all output filenames (prevents accidental file overwriting)\n\n")

optional.add_argument("-pr", "--pos_range", type=int, default=[0,0], nargs=2, help="Minimum and maximum positions within the reference sequence to display.  Specified as a pair of integer numbers in the order minimum maximum (e.g. -pr 100 200).  If unspecified, the full reference range will be used.\n\n")

optional.add_argument("-id", "--im_dims", type=int, default=[def_im_w, def_im_h], nargs=2, help="Pixel dimensions of the output .svg image.  Specified as a pair of integer numbers in the order width height (e.g. -id 800 200).  Default: \"%i %i\".\n\n" % (def_im_w, def_im_h))

optional.add_argument("-f", "--font", type=str, default=def_font, help="Font to use for all text.  Can be any font currently installed on this computer.  Default: \"%s\".\n\n" % def_font)

optional.add_argument("-m_rp", "--map_rel_pos", type=float, default=[def_map_rel_top, def_map_rel_height, def_map_rel_left, def_map_rel_width], nargs=4, help="Position and size of strand linkage plot in output image, relative to the top-left corner.  Positions correspond directly to map, so allowances must be made for labels.  Specified as a list of 4 floating-point numbers in the order top height left width (e.g. -mrp 0.3 0.6 0.05 0.9).  Default: \"%.2f %.2f %.2f %.2f\".\n\n" % (def_map_rel_top, def_map_rel_height, def_map_rel_left, def_map_rel_width))

optional.add_argument("-d_m", "--dna_mode", type=slpw.DNA_MODE, default=def_dna_mode, choices=list(slpw.DNA_MODE), help="Rendering mode for DNA strands.  \"line\" displays DNA as a pair of solid lines with width controlled by --dna_size.  \"seq\" displays DNA as a pair of letter sequences with font size controlled by --dna_size (Note: position range must be sufficiently narrow in sequence mode to prevent text overlap).  \"none\" displays nothing.  Default: \"%s\".\n\n" % def_dna_mode)

optional.add_argument("-d_s", "--dna_size", type=int, default=def_dna_size, help="Width of rendered DNA lines (when in --dna_mode \"line\") or font size of DNA sequence letters (when in --dna_mode \"seq\").  Default: \"%.1f\".\n\n" % def_dna_size)

optional.add_argument("-d_c", "--dna_colour", type=str, default=def_dna_colour, help="Colour of rendered DNA lines or text sequences.  Can be specified as colour names (e.g. \"black\"), as hex values (e.g. \"#16C3D6\" for a light blue) or as RGB values in the range 0-255 (e.g. \"rgb(128,0,128)\" for purple).  Default: \"%s\".\n\n" % def_dna_colour)

optional.add_argument("-d_rg", "--dna_rel_gap", type=float, default=def_dna_rel_gap, help="Gap between DNA sequence/line and strandlinkageplot.  Specified as a fraction of the image height.  Default: \"%i\".\n\n" % def_dna_rel_gap)

optional.add_argument("-el_v", "--end_label_vis", type=SHOWHIDE, default=def_end_label_vis, choices=list(SHOWHIDE), help="Controls whether the DNA end labels (\"5′\" and \"3′\") are rendered.  Must be either \"show\" or \"hide\" (e.g. -el_v \"show\").  Default: \"%s\".\n\n" % def_end_label_vis)

optional.add_argument("-el_s", "--end_label_size", type=int, default=def_end_label_size, help="Font size of DNA end labels.  Default: \"%.1f\".\n\n" % def_end_label_size)

optional.add_argument("-el_c", "--end_label_colour", type=str, default=def_end_label_colour, help="Colour of the rendered DNA end labels.  Can be specified as colour names (e.g. \"black\"), as hex values (e.g. \"#16C3D6\" for a light blue) or as RGB values in the range 0-255 (e.g. \"rgb(128,0,128)\" for purple).  Default: \"%s\".\n\n" % def_end_label_colour)

optional.add_argument("-el_rg", "--end_label_rel_gap", type=float, default=def_end_label_rel_gap, help="Gap between end labels and the rendered DNA strands.  Specified as a fraction of the image width.  Default: \"%i\".\n\n" % def_end_label_rel_gap)

optional.add_argument("-el_p", "--end_label_position", type=slpw.VPOS, default=def_end_label_position, choices=list(slpw.VPOS), help="Controls whether end labels are vertically centre-aligned with the DNA or aligned inside (top strand top-aligned with DNA and bottom strand bottom-aligned).  Must be either \"centre\" or \"inside\" (e.g. -elp \"centre\").  Default: \"%s\".\n\n" % def_end_label_position)

optional.add_argument("-g_v", "--grid_vis", type=SHOWHIDE, default=def_grid_vis, choices=list(SHOWHIDE), help="Controls whether grid lines are rendered.  Must be either \"show\" or \"hide\" (e.g. -g_v \"show\").  Default: \"%s\".\n\n" % def_grid_vis)

optional.add_argument("-g_s", "--grid_size", type=int, default=def_grid_size, help="Width of grid lines.  Default: \"%.1f\".\n\n" % def_grid_size)

optional.add_argument("-g_c", "--grid_colour", type=str, default=def_grid_colour, help="Colour of the rendered grid lines.  Can be specified as colour names (e.g. \"black\"), as hex values (e.g. \"#16C3D6\" for a light blue) or as RGB values in the range 0-255 (e.g. \"rgb(128,0,128)\" for purple).  Default: \"%s\".\n\n" % def_grid_colour)

optional.add_argument("-g_i", "--grid_interval", type=int, default=def_grid_interval, help="Interval between adjacent grid lines.  Default: \"%i\".\n\n" % def_grid_interval)

optional.add_argument("-gl_v", "--grid_label_vis", type=SHOWHIDE, default=def_grid_label_vis, choices=list(SHOWHIDE), help="Controls whether grid labels are rendered.  Grid labels are always shown above the top strand and are rotated vertically.  Must be either \"show\" or \"hide\" (e.g. -gl_v \"show\").  Default: \"%s\".\n\n" % def_grid_label_vis)

optional.add_argument("-gl_s", "--grid_label_size", type=int, default=def_grid_label_size, help="Font size of grid labels.  Default: \"%.1f\".\n\n" % def_grid_label_size)

optional.add_argument("-gl_c", "--grid_label_colour", type=str, default=def_grid_label_colour, help="Colour of the rendered grid labels.  Can be specified as colour names (e.g. \"black\"), as hex values (e.g. \"#16C3D6\" for a light blue) or as RGB values in the range 0-255 (e.g. \"rgb(128,0,128)\" for purple).  Default: \"%s\".\n\n" % def_grid_label_colour)

optional.add_argument("-gl_i", "--grid_label_interval", type=int, default=def_grid_label_interval, help="Interval between adjacent grid labels.  Default: \"%i\".\n\n" % def_grid_label_interval)

optional.add_argument("-gl_rg", "--grid_label_rel_gap", type=float, default=def_grid_label_rel_gap, help="Gap between grid labels and the rendered DNA strands.  Specified as a fraction of the image height.  Default: \"%i\".\n\n" % def_grid_label_rel_gap)

optional.add_argument("-c_v", "--cbar_vis", type=SHOWHIDE, default=def_cbar_vis, choices=list(SHOWHIDE), help="Controls whether the colourbar is rendered.  Must be either \"show\" or \"hide\" (e.g. -c_v \"show\").  Default: \"%s\".\n\n" % def_cbar_vis)

optional.add_argument("-c_rp", "--cbar_rel_pos", type=float, default=[def_cbar_rel_left,def_cbar_rel_width], nargs=2, help="Position and size of colourbar in output image, relative to the top-left corner.  Specified as a list of 2 floating-point numbers in the order left width (e.g. -crp 0.9 0.02).  Default: \"%.2f %.2f\".\n\n" % (def_cbar_rel_left, def_cbar_rel_width))

optional.add_argument("-c_s", "--cbar_size", type=str, default=def_cbar_border_size, help="Width of colourbar border.  Default: \"%s\".\n\n" % def_cbar_border_size)

optional.add_argument("-cl_v", "--cbar_label_vis", type=SHOWHIDE, default=def_cbar_label_vis, choices=list(SHOWHIDE), help="Controls whether colourbar labels are rendered.  Must be either \"show\" or \"hide\" (e.g. -cl_v \"show\").  Default: \"%s\".\n\n" % def_cbar_label_vis)

optional.add_argument("-cl_s", "--cbar_label_size", type=int, default=def_cbar_label_size, help="Font size of colourbar labels.  Default: \"%.1f\".\n\n" % def_cbar_label_size)

optional.add_argument("-cl_c", "--cbar_label_colour", type=str, default=def_cbar_label_colour, help="Colour of the rendered colourbar labels.  Can be specified as colour names (e.g. \"black\"), as hex values (e.g. \"#16C3D6\" for a light blue) or as RGB values in the range 0-255 (e.g. \"rgb(128,0,128)\" for purple).  Default: \"%s\".\n\n" % def_cbar_label_colour)

optional.add_argument("-cl_i", "--cbar_label_interval", type=int, default=def_cbar_label_interval, help="Interval between labels shown in colourbar.  Default: \"%.1f\".\n\n" % def_cbar_label_interval)

optional.add_argument("-cl_rg", "--cbar_label_rel_gap", type=float, default=def_cbar_label_rel_gap, help="Gap between colourbar labels and the colourbar itself.  Specified as a fraction of the image height.  Default: \"%i\".\n\n" % def_cbar_label_rel_gap)

optional.add_argument("-e_mis", "--event_min_size", type=float, default=def_event_min_size, help="Width of line corresponding to the least frequent cleavage event.  Default: \"%.1f\".\n\n" % def_event_min_size)

optional.add_argument("-e_mas", "--event_max_size", type=float, default=def_event_max_size, help="Width of line corresponding to the most frequent cleavage event.  Default: \"%.1f\".\n\n" % def_event_max_size)

optional.add_argument("-e_c", "--event_colourmap", type=str, default=def_event_colourmap, help="Matplotlib colourmap to use for event plotting.  Must be a Matplotlib colourmap.  Default: \"%s\".\n\n" % def_event_colourmap)

optional.add_argument("-e_r", "--event_range", type=int, default=[def_event_min_range,def_event_max_range], nargs=2, help="Range of values colourscale will span (specified as percentage of all events).  For automatic range selection, set both values to -1 (e.g. -e_r -1 -1).  Default: \"%i %i\".\n\n" % (def_event_min_range,def_event_max_range))

optional.add_argument("-e_orv", "--event_outside_range_vis", type=SHOWHIDE, default=def_event_outside_range_vis, choices=list(SHOWHIDE), help="Controls whether events not entirely within the displayed range are included in the strandlinkageplot.  Must be either \"show\" or \"hide\" (e.g. -e_orv \"show\").  Default: \"%s\".\n\n" % def_event_outside_range_vis)

optional.add_argument("-e_o", "--event_opacity", type=float, default=def_event_opacity, help="Opacity of each event line.  Values in the range 0-1.  Default: \"%.1f\".\n\n" % def_event_opacity)

optional.add_argument("-e_so", "--event_stack_order", type=int, default=def_event_stack_order, help="Mode for ordering events.  Options: 1 (most frequent at back), 2 (most frequent at front).  Default: \"%.1f\".\n\n" % def_event_stack_order)

optional.add_argument("-h_v", "--hist_vis", type=SHOWHIDE, default=def_hist_vis, choices=list(SHOWHIDE), help="Controls whether histograms are rendered above and below the strandlinkageplot.  Must be either \"show\" or \"hide\" (e.g. -h_v \"show\").  Default: \"%s\".\n\n" % def_hist_vis)

optional.add_argument("-h_r", "--hist_range", type=int, default=[def_hist_min_range,def_hist_max_range], nargs=2, help="Range of values histogram will span (specified as percentage of all events).  For automatic range selection, set both values to -1 (e.g. -er -1 -1).  Default: \"%i %i\".\n\n" % (def_hist_min_range,def_hist_max_range))

optional.add_argument("-h_bw", "--hist_bin_width", type=int, default=def_hist_bin_width, help="Number of positions that will be binned into a single bar on the histogram.  Default: \"%.1f\".\n\n" % def_hist_bin_width)

optional.add_argument("-h_c", "--hist_colour", type=str, default=def_hist_colour, help="Colour of the rendered histogram bars.  Can be specified as colour names (e.g. \"black\"), as hex values (e.g. \"#16C3D6\" for a light blue) or as RGB values in the range 0-255 (e.g. \"rgb(128,0,128)\" for purple).  Default: \"%s\".\n\n" % def_hist_colour)

optional.add_argument("-h_rh", "--hist_rel_height", type=float, default=def_hist_rel_height, help="Height of the histogram plots.  Specified as a fraction of the image height.  Default: \"%i\".\n\n" % def_hist_rel_height)

optional.add_argument("-h_rg", "--hist_rel_gap", type=float, default=def_hist_rel_gap, help="Gap between histogram and the strandlinkageplot.  Specified as a fraction of the image height.  Default: \"%i\".\n\n" % def_hist_rel_gap)

optional.add_argument("-h_pbg", "--hist_pc_bar_gap", type=float, default=def_hist_pc_bar_gap, help="Gap between adjacent histogram bars.  Specified as a percentage of the total bar width (e.g. \"0\" will have no gap and \"50\" will have half-width bars).  Default: \"%i\".\n\n" % def_hist_pc_bar_gap)

optional.add_argument("-h_o", "--hist_overhang", type=float, default=def_hist_overhang, help="Additional width of histogram plot area along x-axis.  This allows the histogram to be aligned with the outer edges of the DNA sequence text.  Specified in pixel units..  Default: \"%i\".\n\n" % def_hist_overhang)

optional.add_argument("-h_lv", "--hist_log_vis", type=SHOWHIDE, default=def_hist_log_vis, choices=list(SHOWHIDE), help="Controls whether histogram y-axis is shown as logarithmic values.  Must be either \"show\" or \"hide\" (e.g. -hl_v \"show\").  Default: \"%s\".\n\n" % def_hist_log_vis)

optional.add_argument("-hl_v", "--hist_label_vis", type=SHOWHIDE, default=def_hist_label_vis, choices=list(SHOWHIDE), help="Controls whether histogram labels are rendered.  Must be either \"show\" or \"hide\" (e.g. -hl_v \"show\").  Default: \"%s\".\n\n" % def_hist_label_vis)

optional.add_argument("-hl_s", "--hist_label_size", type=int, default=def_hist_label_size, help="Font size of histogram labels.  Default: \"%.1f\".\n\n" % def_hist_label_size)

optional.add_argument("-hl_c", "--hist_label_colour", type=str, default=def_hist_label_colour, help="Colour of the rendered histogram labels.  Can be specified as colour names (e.g. \"black\"), as hex values (e.g. \"#16C3D6\" for a light blue) or as RGB values in the range 0-255 (e.g. \"rgb(128,0,128)\" for purple).  Default: \"%s\".\n\n" % def_hist_label_colour)

optional.add_argument("-hl_i", "--hist_label_interval", type=int, default=def_hist_label_interval, help="Interval between labels shown in histogram.  Default: \"%.1f\".\n\n" % def_hist_label_interval)

optional.add_argument("-hl_rg", "--hist_label_rel_gap", type=float, default=def_hist_label_rel_gap, help="Gap between histogram labels and side of the histogram.  Specified as a fraction of the image width.  Default: \"%i\".\n\n" % def_hist_label_rel_gap)

optional.add_argument("-hl_p", "--hist_label_position", type=slpw.HPOS, default=def_hist_label_position, choices=list(slpw.HPOS), help="Controls whether histograms labels are rendered to the left or right of the histogram.  Must be either \"left\" or \"right\" (e.g. -hl_p \"left\").  Default: \"%s\".\n\n" % def_hist_label_position)

optional.add_argument("-hl_zv", "--hist_label_zero_vis", type=SHOWHIDE, default=def_hist_label_zero_vis, choices=list(SHOWHIDE), help="Controls whether the 0 histogram label is rendered (having this hidden allows the histogram to be closer to the DNA sequence).  Must be either \"show\" or \"hide\" (e.g. -hl_zv \"show\").  Default: \"%s\".\n\n" % def_hist_label_zero_vis)

optional.add_argument("-hg_v", "--hist_grid_vis", type=SHOWHIDE, default=def_grid_vis, choices=list(SHOWHIDE), help="Controls whether horizontal histogram grid lines are rendered.  Must be either \"show\" or \"hide\" (e.g. -hg_v \"show\").  Default: \"%s\".\n\n" % def_hist_grid_vis)

optional.add_argument("-hg_s", "--hist_grid_size", type=int, default=def_hist_grid_size, help="Width of histogram grid lines.  Default: \"%.1f\".\n\n" % def_hist_grid_size)

optional.add_argument("-hg_c", "--hist_grid_colour", type=str, default=def_hist_grid_colour, help="Colour of the rendered histogram grid lines.  Can be specified as colour names (e.g. \"black\"), as hex values (e.g. \"#16C3D6\" for a light blue) or as RGB values in the range 0-255 (e.g. \"rgb(128,0,128)\" for purple).  Default: \"%s\".\n\n" % def_hist_grid_colour)

optional.add_argument("-hg_i", "--hist_grid_interval", type=int, default=def_hist_grid_interval, help="Interval between horizontal grid lines shown in histogram.  Default: \"%.1f\".\n\n" % def_hist_grid_interval)

optional.add_argument("-hn_v", "--hist_name_vis", type=SHOWHIDE, default=def_hist_name_vis, choices=list(SHOWHIDE), help="Controls whether histogram names (e.g. \"Sum\", \"5′ OH\", etc.) are rendered to the side of each histogram.  Must be either \"show\" or \"hide\" (e.g. -hn_v \"show\").  Default: \"%s\".\n\n" % def_hist_name_vis)

optional.add_argument("-hn_s", "--hist_name_size", type=int, default=def_hist_name_size, help="Font size of histogram names.  Default: \"%.1f\".\n\n" % def_hist_name_size)

optional.add_argument("-hn_c", "--hist_name_colour", type=str, default=def_hist_name_colour, help="Colour of the rendered histogram names.  Can be specified as colour names (e.g. \"black\"), as hex values (e.g. \"#16C3D6\" for a light blue) or as RGB values in the range 0-255 (e.g. \"rgb(128,0,128)\" for purple).  Default: \"%s\".\n\n" % def_hist_name_colour)

optional.add_argument("-hn_rg", "--hist_name_rel_gap", type=float, default=def_hist_name_rel_gap, help="Gap between histogram names and side of the histogram.  Specified as a fraction of the image width.  Default: \"%i\".\n\n" % def_hist_name_rel_gap)

optional.add_argument("-hn_p", "--hist_name_position", type=slpw.HPOS, default=def_hist_name_position, choices=list(slpw.HPOS), help="Controls whether histograms names are rendered to the left or right of the histogram.  Must be either \"left\" or \"right\" (e.g. -hn_p \"left\").  Default: \"%s\".\n\n" % def_hist_name_position)

optional.add_argument("-sh_v", "--splithist_vis", type=SHOWHIDE, default=def_splithist_vis, choices=list(SHOWHIDE), help="Controls whether split histograms (one each for 3′, 5′ and blunt events) are rendered above and below the strandlinkageplot.  Must be either \"show\" or \"hide\" (e.g. -sh_v \"show\").  Default: \"%s\".\n\n" % def_splithist_vis)

optional.add_argument("-sh_r", "--splithist_range", type=int, default=[def_splithist_min_range,def_splithist_max_range], nargs=2, help="Range of values split histograms will span (specified as percentage of all events).  For automatic range selection, set both values to -1 (e.g. -er -1 -1).  Default: \"%i %i\".\n\n" % (def_splithist_min_range,def_splithist_max_range))

optional.add_argument("-sh_rh", "--splithist_rel_height", type=float, default=def_splithist_rel_height, help="Height of the split histogram plots.  Specified as a fraction of the image height.  Default: \"%i\".\n\n" % def_splithist_rel_height)

optional.add_argument("-shl_i", "--splithist_label_interval", type=int, default=def_splithist_label_interval, help="Interval between labels shown in split histograms.  Default: \"%.1f\".\n\n" % def_splithist_label_interval)

optional.add_argument("-shg_i", "--splithist_grid_interval", type=int, default=def_splithist_grid_interval, help="Interval between horizontal grid lines shown in split histograms.  Default: \"%.1f\".\n\n" % def_splithist_grid_interval)

args = parser.parse_args()

end_label_show = args.end_label_vis is SHOWHIDE.SHOW
grid_show = args.grid_vis is SHOWHIDE.SHOW
grid_label_show = args.grid_label_vis is SHOWHIDE.SHOW
cbar_show = args.cbar_vis is SHOWHIDE.SHOW
cbar_label_show = args.cbar_label_vis is SHOWHIDE.SHOW
event_outside_range_vis = args.event_outside_range_vis is SHOWHIDE.SHOW
hist_show = args.hist_vis is SHOWHIDE.SHOW
hist_log_show = args.hist_log_vis is SHOWHIDE.SHOW
hist_label_show = args.hist_label_vis is SHOWHIDE.SHOW
hist_label_zero_show = args.hist_label_zero_vis is SHOWHIDE.SHOW
hist_name_show = args.hist_name_vis is SHOWHIDE.SHOW
hist_grid_show = args.hist_grid_vis is SHOWHIDE.SHOW
splithist_show = args.splithist_vis is SHOWHIDE.SHOW

# Required arguments
pos_range = tuple(args.pos_range) if args.pos_range != [0,0] else None
im_dims = tuple(args.im_dims)
rel_pos = tuple(args.map_rel_pos)
dna_opts = (args.dna_mode,args.dna_size,args.dna_colour,args.dna_rel_gap)
end_label_opts = (end_label_show,args.end_label_size,args.end_label_colour,args.end_label_rel_gap,args.end_label_position)
grid_opts = (grid_show,args.grid_size,args.grid_colour,args.grid_interval)
grid_label_opts = (grid_label_show,args.grid_label_size,args.grid_label_colour,args.grid_label_interval,args.grid_label_rel_gap)
cbar_opts = (cbar_show,args.cbar_rel_pos[0],args.cbar_rel_pos[1],args.cbar_size)
cbar_label_opts = (cbar_label_show,args.cbar_label_size,args.cbar_label_colour,args.cbar_label_interval,args.cbar_label_rel_gap)
event_opts = (args.event_min_size,args.event_max_size,args.event_colourmap,args.event_range[0],args.event_range[1],event_outside_range_vis,args.event_opacity,args.event_stack_order)
hist_opts = (hist_show,args.hist_range[0],args.hist_range[1],args.hist_bin_width,args.hist_colour,args.hist_rel_height,args.hist_rel_gap,args.hist_pc_bar_gap,args.hist_overhang,hist_log_show)
hist_label_opts = (hist_label_show,args.hist_label_size,args.hist_label_colour,args.hist_label_interval,args.hist_label_rel_gap,args.hist_label_position,hist_label_zero_show)
hist_grid_opts = (hist_grid_show,args.hist_grid_size,args.hist_grid_colour,args.hist_grid_interval)
hist_name_opts = (hist_name_show,args.hist_name_size,args.hist_name_colour,args.hist_name_rel_gap,args.hist_name_position)
splithist_opts = (splithist_show,args.splithist_range[0],args.splithist_range[1],args.splithist_rel_height)
splithist_label_opts = args.splithist_label_interval
splithist_grid_opts = args.splithist_grid_interval


writer = slpw.StrandLinkagePlotWriter(im_dims=im_dims, font=args.font, rel_pos=rel_pos, dna_opts=dna_opts, end_label_opts=end_label_opts, grid_opts=grid_opts, grid_label_opts=grid_label_opts, cbar_opts=cbar_opts, cbar_label_opts=cbar_label_opts, event_opts=event_opts, hist_opts=hist_opts, hist_label_opts=hist_label_opts, hist_grid_opts=hist_grid_opts, hist_name_opts=hist_name_opts, splithist_opts=splithist_opts, splithist_label_opts=splithist_label_opts,splithist_grid_opts=splithist_grid_opts)

writer.write_from_file(args.data_path, args.out_path, ref_path=args.ref_path, pos_range=pos_range, append_dt=args.append_datetime)

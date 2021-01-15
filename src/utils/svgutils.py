import datetime as dt
import math
import matplotlib as mpl
import numpy as np
import os
import svgwrite as svg

from utils import csvutils as cu  
from utils import plotutils as pu
from utils import reportutils as ru

from utils.fileutils import FileReader

class EventMapWriter():
    ## CONSTRUCTOR

    def __init__(self, dims=(800,200), rel_pos=(0.3,0.9,0.05,0.95), dna_opts=(False,2,"black"), end_label_opts=(True,20,"black",0.01), grid_opts=(True,1,"lightgray",100), grid_label_opts=(True,12,"lightgray",500), event_opts=(2,"blue","red")):
        self._doc_w = dims[0]
        self._doc_h = dims[1]

        self._dna_rel_top = rel_pos[0]
        self._dna_rel_bottom = rel_pos[1]
        self._dna_rel_left = rel_pos[2]
        self._dna_rel_right = rel_pos[3]

        self._dna_x1 = self._doc_w*self._dna_rel_left
        self._dna_x2 = self._doc_w*self._dna_rel_right
        self._dna_y1 = self._doc_h*self._dna_rel_top
        self._dna_y2 = self._doc_h*self._dna_rel_bottom
        
        self._dna_seq_show = dna_opts[0]
        self._dna_seq_size = dna_opts[1]
        self._dna_seq_colour = dna_opts[2]

        self._grid_show = grid_opts[0]
        self._grid_size = grid_opts[1]
        self._grid_colour = grid_opts[2]
        self._grid_increment = grid_opts[3]

        self._grid_label_show = grid_label_opts[0]
        self._grid_label_size = grid_label_opts[1]
        self._grid_label_colour = grid_label_opts[2]
        self._grid_label_increment = grid_label_opts[3]

        self._end_label_show = end_label_opts[0]
        self._end_label_size = end_label_opts[1]
        self._end_label_colour = end_label_opts[2]
        self._rel_end_label_gap = end_label_opts[3]
        self._end_label_gap = self._doc_w*self._rel_end_label_gap

        self._event_max_size = event_opts[0]
        self._event_colour_1 = event_opts[1]
        self._event_colour_2 = event_opts[2]

    ## GETTERS AND SETTERS

    def set_dna_seq_show(self, dna_seq_show):
        self._dna_seq_show = dna_seq_show


    ## PUBLIC METHODS

    def write_event_map_from_file(self, csv_path, out_path="", ref_path="", pos_range=None, append_dt=False):
        # Loading results from CSV
        cr = cu.CSVReader()
        results = cr.read_individual(csv_path)
        freq = ru.get_full_sequence_frequency(results)

        # If no output path is specified, use the same as the input csv file
        if out_path == "":
            out_path = csv_path

        # Loading reference sequence if path provided
        if ref_path == "":
            ref = None
        else:
            filereader = FileReader(verbose=False)
            ref = filereader.read_sequence(ref_path)[0][0]
        
        self.write_event_map(out_path, freq, ref=ref, pos_range=pos_range, append_dt=append_dt)

    def write_event_map(self, out_path, freq, ref=None, pos_range=None, append_dt=False):
        if pos_range is None:
            pos_min = 0
            if ref is None:
                # Rounding up to the nearest 10
                pos_max = get_max_event_pos(freq)
                pos_max = math.ceil(pos_max/10)*10
            else:
                pos_max = len(ref)-1
        else:
            pos_min = pos_range[0]
            pos_max = pos_range[1]

        # Checking pos_min and pos_max are different (to prevent divide by zero errors)
        if pos_min == pos_max:
            print("WARNING: Min and max sequence positions must be different")
            return

        # Creating output SVG document
        root_name = os.path.splitext(out_path)[0]
        datetime_str = dt.datetime.now().strftime("_%Y-%m-%d_%H-%M-%S") if append_dt else ""
        outname = root_name+ '_eventmap' + datetime_str + '.' + 'svg'
        dwg = svg.Drawing(outname, size = ("%spx" % self._doc_w, "%spx" % self._doc_h))

        if self._grid_show:
            self._add_grid_lines(dwg, pos_min, pos_max)

        if self._grid_label_show:
            self._add_grid_labels(dwg, pos_min, pos_max)
        
        if self._end_label_show:
            self._add_end_labels(dwg)

        self._add_event_lines(dwg, pos_min, pos_max, freq)
        self._add_dna(dwg, pos_min, pos_max, ref=ref)

        # Writing SVG to file
        dwg.save()

    def _add_grid_lines(self, dwg, pos_min, pos_max):
        # Adding horizontal grid line
        grid_yc = (self._dna_y1+self._dna_y2)/2
        dwg.add(svg.shapes.Line((self._dna_x1, grid_yc), (self._dna_x2, grid_yc), stroke=self._grid_colour, stroke_width=self._grid_size))

        # Getting range of vertical grid lines
        grid_min = self._grid_increment*math.ceil(pos_min/self._grid_increment)
        grid_max = self._grid_increment*math.floor(pos_max/self._grid_increment)
        if grid_max%self._grid_increment == 0:
            grid_max += 1

        # Adding vertical grid lines
        for grid_pos in range(grid_min, grid_max, self._grid_increment):
            grid_x = self._dna_x1 + (self._dna_x2-self._dna_x1)*((grid_pos-pos_min)/(pos_max-pos_min))
            grid_t_y = self._dna_y1+self._dna_seq_size/2
            grid_b_y = self._dna_y2-self._dna_seq_size/2
            dwg.add(svg.shapes.Line((grid_x, grid_t_y), (grid_x, grid_b_y), stroke=self._grid_colour, stroke_width=self._grid_size))

    def _add_grid_labels(self, dwg, pos_min, pos_max):
        # Getting range of grid label positions
        grid_label_min = self._grid_label_increment*math.ceil(pos_min/self._grid_label_increment)
        grid_label_max = self._grid_label_increment*math.floor(pos_max/self._grid_label_increment)
        if grid_label_max%self._grid_label_increment == 0:
            grid_label_max += 1

        # Adding grid labels
        for grid_label_pos in range(grid_label_min, grid_label_max, self._grid_label_increment):
            grid_label_x = self._dna_x1 + (self._dna_x2-self._dna_x1)*((grid_label_pos-pos_min)/(pos_max-pos_min))
            grid_label_y = self._dna_y1-self._dna_seq_size/2-self._end_label_gap
            rot = "rotate(%i,%i,%i)" % (-90,grid_label_x,grid_label_y)
            dwg.add(svg.text.Text(str(grid_label_pos), insert=(grid_label_x,grid_label_y), transform=rot, style="text-anchor:start; baseline-shift:-50%", font_size=self._grid_label_size, fill=self._grid_label_colour))

    def _add_dna(self, dwg, pos_min, pos_max, ref=None):
        if self._dna_seq_show and ref is None:
            print("WARNING: No reference path provided.  DNA rendered as solid lines.")
            dna_seq_show = False

        if self._dna_seq_show:
            # Draw DNA as text sequence
            dna_pos_increment = (self._dna_x2-self._dna_x1)/(pos_max-pos_min)
            ref_rc = ref.reverse_complement()
            for dna_pos in range(pos_min,pos_max+1):
                dna_seq_x = self._dna_x1+dna_pos_increment*(dna_pos-pos_min)
                dna_seq_y1 = self._dna_y1+self._dna_seq_size*0.375
                dna_seq_y2 = self._dna_y2+self._dna_seq_size*0.375
                dwg.add(svg.text.Text(str(ref[dna_pos]), insert=(dna_seq_x, dna_seq_y1), style="text-anchor:middle", font_size=self._dna_seq_size, fill=self._dna_seq_colour))
                dwg.add(svg.text.Text(str(ref_rc[dna_pos]), insert=(dna_seq_x, dna_seq_y2), style="text-anchor:middle", font_size=self._dna_seq_size, fill=self._dna_seq_colour))

        else:
            # Draw DNA as lines
            dwg.add(svg.shapes.Line((self._dna_x1, self._dna_y1), (self._dna_x2, self._dna_y1), stroke=self._dna_seq_colour, stroke_width=self._dna_seq_size, style="stroke-linecap:square"))
            dwg.add(svg.shapes.Line((self._dna_x1, self._dna_y2), (self._dna_x2, self._dna_y2), stroke=self._dna_seq_colour, stroke_width=self._dna_seq_size, style="stroke-linecap:square"))

    def _add_end_labels(self, dwg):
        end_label_x1 = self._dna_x1-self._end_label_gap
        end_label_x2 = self._dna_x2+self._end_label_gap
        end_label_y1 = self._dna_y1+self._end_label_size*0.375
        end_label_y2 = self._dna_y2+self._end_label_size*0.375

        dwg.add(svg.text.Text("5'", insert=(end_label_x1, end_label_y1), style="text-anchor:end", font_size=self._end_label_size, fill=self._end_label_colour))
        dwg.add(svg.text.Text("3'", insert=(end_label_x2, end_label_y1), style="text-anchor:start", font_size=self._end_label_size, fill=self._end_label_colour))
        dwg.add(svg.text.Text("3'", insert=(end_label_x1, end_label_y2), style="text-anchor:end", font_size=self._end_label_size, fill=self._end_label_colour))
        dwg.add(svg.text.Text("5'", insert=(end_label_x2, end_label_y2), style="text-anchor:start", font_size=self._end_label_size, fill=self._end_label_colour))

    def _add_event_lines(self, dwg, pos_min, pos_max, freq):
        event_colour_1 = np.array(mpl.colors.to_rgb(self._event_colour_1))
        event_colour_2 = np.array(mpl.colors.to_rgb(self._event_colour_2))

        # Adding cleavage event lines
        total = sum(freq.values())
        max_events = max(freq.values())
        min_events = min(freq.values())
        diff_events = max_events-min_events

        # Preventing divide by zero errors
        if diff_events == 0:
            diff_events = 1

        for (cleavage_site_t, cleavage_site_b) in freq.keys():        
            norm_count = (freq.get((cleavage_site_t, cleavage_site_b))-min_events)/diff_events
            # norm_count = freq.get((cleavage_site_t, cleavage_site_b)/total
            
            # Adding line (adding width 1 to ensure everything is visible)
            event_width = self._event_max_size*norm_count+1

            event_t_xc = self._dna_x1 + (self._dna_x2-self._dna_x1)*((cleavage_site_t-pos_min)/(pos_max-pos_min))
            event_t_x1 = event_t_xc-event_width
            event_t_x2 = event_t_xc+event_width
            event_t_y = self._dna_y1+self._dna_seq_size/2

            event_b_xc = self._dna_x1 + (self._dna_x2-self._dna_x1)*((cleavage_site_b-pos_min)/(pos_max-pos_min))
            event_b_x1 = event_b_xc-event_width
            event_b_x2 = event_b_xc+event_width
            event_b_y = self._dna_y2-self._dna_seq_size/2

            col = mpl.colors.to_hex((1-norm_count)*event_colour_1 + norm_count*event_colour_2)

            dwg.add(svg.shapes.Polygon(points=[(event_t_x1,event_t_y), (event_t_x2,event_t_y), (event_b_x2,event_b_y), (event_b_x1,event_b_y)], fill=col))

def get_max_event_pos(freq):
    max_pos = 0
    for (cleavage_site_t, cleavage_site_b) in freq.keys():
        max_pos = math.max(max_pos,cleavage_site_t)
        max_pos = math.max(max_pos,cleavage_site_b)
    
    return max_pos
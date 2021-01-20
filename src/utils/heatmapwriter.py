# TODO: Row and column sums
# TODO: Fill empty cells with background colour
# TODO: Have colour normalisation to visible position range only

import datetime as dt
import math
from matplotlib import cm
import os
import svgwrite as svg
import sys

from enum import Enum
from utils import csvutils as cu  
from utils import fileutils as fu
from utils import reportutils as ru


class DNA_MODE(Enum):
    LINE = 'line'
    NONE = 'none'
    SEQUENCE = 'seq'

    def __str__(self):
        return str(self.value)


class HeatMapWriter():
    ## CONSTRUCTOR

    def __init__(self, im_dim=800, rel_pos=(0.1,0.1,0.85), border_opts=(True,1,"black"), grid_opts=(True,1,"lightgray",1), grid_label_opts=(True,12,"lightgray",10,10), event_colourmap="cool", event_label_opts=(False,10,"white")):
        self._im_dim = im_dim
        
        self._map_rel_top = rel_pos[0]
        self._map_rel_left = rel_pos[1]
        self._map_rel_size = rel_pos[2]
        
        self._border_show = border_opts[0]
        self._border_size = border_opts[1]
        self._border_colour = border_opts[2]

        self._grid_show = grid_opts[0]
        self._grid_size = grid_opts[1]
        self._grid_colour = grid_opts[2]
        self._grid_interval = grid_opts[3]

        self._grid_label_show = grid_label_opts[0]
        self._grid_label_size = grid_label_opts[1]
        self._grid_label_colour = grid_label_opts[2]
        self._grid_label_interval = grid_label_opts[3]
        self._grid_label_gap = grid_label_opts[4]

        self._event_colourmap = event_colourmap
        
        self._event_label_show = event_label_opts[0]
        self._event_label_size = event_label_opts[1]
        self._event_label_colour = event_label_opts[2]
               

    ## PUBLIC METHODS

    def write_map_from_file(self, csv_path, out_path, ref_path="", pos_ranges=None, append_dt=False):
        # Loading results from CSV
        cr = cu.CSVReader()
        results = cr.read_individual(csv_path)
        freq = ru.get_full_sequence_frequency(results)

        # Loading reference sequence if path provided
        if ref_path == "":
            ref = None
        else:
            filereader = fu.FileReader(verbose=False)
            ref = filereader.read_sequence(ref_path)[0][0]
            
            self.write_map(out_path, freq, ref=ref, pos_ranges=pos_ranges, append_dt=append_dt)

    def write_map(self, out_path, freq, ref=None, pos_ranges=None, append_dt=False):
        if pos_ranges is None:
            pos_t_min = 0
            pos_b_min = 0
            if ref is None:
                # Rounding up to the nearest 10
                (pos_t_min,pos_t_max,pos_b_min,pos_b_max) = get_event_pos_ranges(freq,round=10)
                
            else:
                pos_t_max = len(ref)-1
                pos_b_max = len(ref)-1
        else:
            pos_t_min = pos_ranges[0]
            pos_t_max = pos_ranges[1]
            pos_b_min = pos_ranges[2]
            pos_b_max = pos_ranges[3]

        # Updating the pos_ranges variable (neater to pass as an argument)
        pos_ranges = (pos_t_min, pos_t_max, pos_b_min, pos_b_max)
         
        # Checking pos_min and pos_max are different (to prevent divide by zero errors)
        if pos_t_min == pos_t_max or pos_b_min == pos_b_max:
            print("WARNING: Min and max sequence positions must be different")
            return

        # Creating output SVG document
        root_name = os.path.splitext(out_path)[0]
        datetime_str = dt.datetime.now().strftime("_%Y-%m-%d_%H-%M-%S") if append_dt else ""
        outname = root_name+ datetime_str + '.' + 'svg'
        dwg = svg.Drawing(outname, size = ("%spx" % self._im_dim, "%spx" % self._im_dim))

        # Defining limits of heat map
        pos_t_range = pos_t_max-pos_t_min
        pos_b_range = pos_b_max-pos_b_min
        max_range = max(pos_t_range, pos_b_range)
        pos_t_rel_size = pos_t_range/max_range
        pos_b_rel_size = pos_b_range/max_range

        map_x1 = self._im_dim*self._map_rel_left
        map_x2 = map_x1 + self._im_dim*self._map_rel_size*pos_t_rel_size
        map_y1 = self._im_dim*self._map_rel_top
        map_y2 = map_y1 + self._im_dim*self._map_rel_size*pos_b_rel_size
        map_xy = (map_x1, map_y1, map_x2, map_y2)

        self._add_events(dwg, pos_ranges, map_xy, freq)

        if self._grid_show:
            self._add_grid_lines(dwg, pos_ranges, map_xy)

        if self._grid_label_show:
            self._add_grid_labels(dwg, pos_ranges, map_xy)
                
        if self._border_show:
            self._add_border(dwg, map_xy)

        # Writing SVG to file
        dwg.save()

    def _add_border(self, dwg, map_xy):
        (map_x1, map_y1, map_x2, map_y2) = map_xy
        map_w = map_x2-map_x1
        map_h = map_y2-map_y1

        # Adding border rectangle
        dwg.add(svg.shapes.Rect(insert=(map_x1,map_y1), size=(map_w,map_h), fill='none', stroke=self._border_colour, stroke_width=self._border_size))
        
    def _add_grid_lines(self, dwg, pos_ranges, map_xy):
        (pos_t_min, pos_t_max, pos_b_min, pos_b_max) = pos_ranges
        (map_x1, map_y1, map_x2, map_y2) = map_xy

        # Getting range of vertical grid lines
        grid_pos_t_min = self._grid_interval*math.ceil(pos_t_min/self._grid_interval)
        grid_pos_t_max = self._grid_interval*math.floor(pos_t_max/self._grid_interval)
        if grid_pos_t_max%self._grid_interval == 0:
            grid_pos_t_max += 1

        # Adding vertical grid lines
        for grid_pos_t in range(grid_pos_t_min, grid_pos_t_max, self._grid_interval):
            grid_x = map_x1 + (map_x2-map_x1)*((grid_pos_t-pos_t_min)/(pos_t_max-pos_t_min+1))
            dwg.add(svg.shapes.Line((grid_x, map_y1), (grid_x, map_y2), stroke=self._grid_colour, stroke_width=self._grid_size))

        # Getting range of horizontal grid lines
        grid_pos_b_min = self._grid_interval*math.ceil(pos_b_min/self._grid_interval)
        grid_pos_b_max = self._grid_interval*math.floor(pos_b_max/self._grid_interval)
        if grid_pos_b_max%self._grid_interval == 0:
            grid_pos_b_max += 1

        # Adding horizontal grid lines
        for grid_pos_b in range(grid_pos_b_min, grid_pos_b_max, self._grid_interval):
            grid_y = map_y1 + (map_y2-map_y1)*((grid_pos_b-pos_b_min)/(pos_b_max-pos_b_min+1))
            dwg.add(svg.shapes.Line((map_x1, grid_y), (map_x2, grid_y), stroke=self._grid_colour, stroke_width=self._grid_size))

    def _add_grid_labels(self, dwg, pos_ranges, map_xy):
        (pos_t_min, pos_t_max, pos_b_min, pos_b_max) = pos_ranges
        (map_x1, map_y1, map_x2, map_y2) = map_xy

        # Getting range of vertical grid label positions
        grid_label_pos_t_min = self._grid_label_interval*math.ceil(pos_t_min/self._grid_label_interval)
        grid_label_pos_t_max = self._grid_label_interval*math.floor(pos_t_max/self._grid_label_interval)
        if grid_label_pos_t_max%self._grid_label_interval == 0:
            grid_label_pos_t_max += 1

        # Adding vertical line grid labels
        for grid_label_pos_t in range(grid_label_pos_t_min, grid_label_pos_t_max, self._grid_label_interval):
            grid_label_x = map_x1 + (map_x2-map_x1)*(((grid_label_pos_t+0.5)-pos_t_min)/(pos_t_max-pos_t_min+1)) + self._grid_label_size*0.375
            grid_label_y = map_y1-self._border_size/2-self._grid_label_gap
            rot = "rotate(%i,%i,%i)" % (-90,grid_label_x,grid_label_y)
            dwg.add(svg.text.Text(str(grid_label_pos_t), insert=(grid_label_x,grid_label_y), transform=rot, style="text-anchor:start", font_size=self._grid_label_size, fill=self._grid_label_colour))

        # Getting range of horizontal grid label positions
        grid_label_pos_b_min = self._grid_label_interval*math.ceil(pos_b_min/self._grid_label_interval)
        grid_label_pos_b_max = self._grid_label_interval*math.floor(pos_b_max/self._grid_label_interval)
        if grid_label_pos_b_max%self._grid_label_interval == 0:
            grid_label_pos_b_max += 1

        # Adding horizontal line grid labels
        for grid_label_pos_t in range(grid_label_pos_b_min, grid_label_pos_b_max, self._grid_label_interval):
            grid_label_x = map_x1-self._border_size/2-self._grid_label_gap
            grid_label_y = map_y1 + (map_y2-map_y1)*(((grid_label_pos_t+0.5)-pos_b_min)/(pos_b_max-pos_b_min+1)) + self._grid_label_size*0.375
            dwg.add(svg.text.Text(str(grid_label_pos_t), insert=(grid_label_x,grid_label_y), style="text-anchor:end", font_size=self._grid_label_size, fill=self._grid_label_colour))

    def _add_events(self, dwg, pos_ranges, map_xy, freq):
        (pos_t_min, pos_t_max, pos_b_min, pos_b_max) = pos_ranges
        (map_x1, map_y1, map_x2, map_y2) = map_xy

        cmap = cm.get_cmap(self._event_colourmap)

        # Adding cleavage event lines
        total = sum(freq.values())
        max_events = max(freq.values())
        min_events = min(freq.values())
        diff_events = max_events-min_events

        event_dim = (map_x2-map_x1)/(pos_t_max-pos_t_min+1)

        # Preventing divide by zero errors
        if diff_events == 0:
            diff_events = 1

        for (cleavage_site_t, cleavage_site_b) in freq.keys():        
            # Checking this event is within the rendered range
            if cleavage_site_t < pos_t_min or cleavage_site_t > pos_t_max or cleavage_site_b < pos_b_min or cleavage_site_b > pos_b_max:
                continue;
                
            norm_count = (freq.get((cleavage_site_t, cleavage_site_b))-min_events)/diff_events
            # norm_count = freq.get((cleavage_site_t, cleavage_site_b)/total
            
            # Adding line (adding width 1 to ensure everything is visible)
            event_x1 = map_x1 + (map_x2-map_x1)*((cleavage_site_t-pos_t_min)/(pos_t_max-pos_t_min+1))
            event_y1 = map_y1 + (map_y2-map_y1)*((cleavage_site_b-pos_b_min)/(pos_b_max-pos_b_min+1))
            rgba = cmap(norm_count)
            col = "rgb(%i,%i,%i)" % (rgba[0]*255,rgba[1]*255,rgba[2]*255)

            dwg.add(svg.shapes.Rect(insert=(event_x1,event_y1), size=(event_dim,event_dim), fill=col, stroke='none'))

            # If selected, adding event percentage labels
            if self._event_label_show:
                # event_label_x = map_x1-self._border_size/2-self._grid_label_gap
                
                event_label_x = event_x1 + event_dim/2
                event_label_y = event_y1 + event_dim/2 + self._event_label_size*0.375 
                event_pc = 100*freq.get((cleavage_site_t, cleavage_site_b))/total
                dwg.add(svg.text.Text("%.1f" % event_pc, insert=(event_label_x,event_label_y), style="text-anchor:middle", font_size=self._event_label_size, fill=self._event_label_colour))


def get_event_pos_ranges(freq, round=1):
    pos_t_min = sys.maxint
    pos_t_max = 0
    pos_b_min = sys.maxint
    pos_b_max = 0

    for (cleavage_site_t, cleavage_site_b) in freq.keys():
        pos_t_min = min(pos_t_min,cleavage_site_t)
        pos_t_max = max(pos_t_max,cleavage_site_t)
        pos_b_min = min(pos_b_min,cleavage_site_b)
        pos_b_max = max(pos_b_max,cleavage_site_b)

    if round != 1:
        pos_t_min = math.floor(pos_t_min/round)*round
        pos_t_max = math.ceil(pos_t_max/round)*round
        pos_b_min = math.floor(pos_b_min/round)*round
        pos_b_max = math.ceil(pos_b_max/round)*round

    print("%i_%i_%i_%i" % (pos_t_min,pos_t_max,pos_b_min,pos_b_max))
    
    return (pos_t_min,pos_t_max,pos_b_min,pos_b_max)

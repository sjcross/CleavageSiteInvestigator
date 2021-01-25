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

    def __init__(self, im_dim=800, rel_pos=(0.1,0.1,0.8), border_opts=(True,1,"black"), axis_label_opts=(True,16,"gray",50), grid_opts=(True,1,"gray",1), grid_label_opts=(True,12,"gray",10,10), event_colourmap="plasma", event_label_opts=(True,10,"invert",True), sum_show=True):
        self._im_dim = im_dim
        
        self._map_rel_top = rel_pos[0]
        self._map_rel_left = rel_pos[1]
        self._map_rel_size = rel_pos[2]
        
        self._border_show = border_opts[0]
        self._border_size = border_opts[1]
        self._border_colour = border_opts[2]

        self._axis_label_show = axis_label_opts[0]
        self._axis_label_size = axis_label_opts[1]
        self._axis_label_colour = axis_label_opts[2]
        self._axis_label_gap = axis_label_opts[3]

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
        self._cmap = cm.get_cmap(event_colourmap)
        
        self._event_label_show = event_label_opts[0]
        self._event_label_size = event_label_opts[1]
        self._event_label_colour = event_label_opts[2]
        self._event_label_zeros_show = event_label_opts[3]

        self._sum_show = sum_show
                       

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
                (pos_t_min,pos_t_max,pos_b_min,pos_b_max) = get_event_pos_ranges(freq)
                
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
        
        dwg = svg.Drawing(outname, size = ("%spx" % self._im_dim, "%spx" % self._im_dim), **{'shape-rendering':'crispEdges'})
        
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

        if self._axis_label_show:
            self._add_axis_labels(dwg, map_xy)

        if self._grid_show:
            self._add_grid_lines(dwg, pos_ranges, map_xy)

        if self._grid_label_show:
            self._add_grid_labels(dwg, pos_ranges, map_xy)
                
        if self._border_show:
            self._add_border(dwg, pos_ranges, map_xy)

        # Writing SVG to file
        dwg.save()

    def _add_axis_labels(self, dwg, map_xy):
        (map_x1, map_y1, map_x2, map_y2) = map_xy

        # Adding top-strand label
        axis_label_x = (map_x2-map_x1)/2 + map_x1
        axis_label_y = map_y1 - self._axis_label_gap
        dwg.add(svg.text.Text("Top strand", insert=(axis_label_x,axis_label_y), style="text-anchor:middle", font_size=self._axis_label_size, fill=self._axis_label_colour))

        # Adding bottom-strand label
        axis_label_x = map_x1 - self._axis_label_gap
        axis_label_y = (map_y2-map_y1)/2 + map_y1
        rot = "rotate(%i,%i,%i)" % (-90,axis_label_x,axis_label_y)
        dwg.add(svg.text.Text("Bottom strand", insert=(axis_label_x,axis_label_y), transform=rot, style="text-anchor:middle", font_size=self._axis_label_size, fill=self._axis_label_colour))

    def _add_border(self, dwg, pos_ranges, map_xy):
        (pos_t_min, pos_t_max, pos_b_min, pos_b_max) = pos_ranges
        (map_x1, map_y1, map_x2, map_y2) = map_xy
        map_w = map_x2-map_x1
        map_h = map_y2-map_y1

        # Adding border rectangle
        dwg.add(svg.shapes.Rect(insert=(map_x1,map_y1), size=(map_w,map_h), fill='none', stroke=self._border_colour, stroke_width=self._border_size))

        # Adding extra border if sum column and row are to be shown
        if self._sum_show:
            # Determine event size
            event_dim = (map_x2-map_x1)/(pos_t_max-pos_t_min+1)

            dwg.add(svg.shapes.Rect(insert=(map_x1+map_w,map_y1), size=(event_dim,map_h), fill='none', stroke=self._border_colour, stroke_width=self._border_size))
            dwg.add(svg.shapes.Rect(insert=(map_x1,map_y1+map_h), size=(map_w,event_dim), fill='none', stroke=self._border_colour, stroke_width=self._border_size))
        
    def _add_grid_lines(self, dwg, pos_ranges, map_xy):
        (pos_t_min, pos_t_max, pos_b_min, pos_b_max) = pos_ranges
        (map_x1, map_y1, map_x2, map_y2) = map_xy

        # Determining extra grid length if sum column and row are to be shown
        event_dim = (map_x2-map_x1)/(pos_t_max-pos_t_min+1) if self._sum_show else 0

        # Getting range of vertical grid lines
        grid_pos_t_min = self._grid_interval*math.ceil(pos_t_min/self._grid_interval)
        grid_pos_t_max = self._grid_interval*math.floor(pos_t_max/self._grid_interval)
        if grid_pos_t_max%self._grid_interval == 0:
            grid_pos_t_max += 1

        # Adding vertical grid lines
        for grid_pos_t in range(grid_pos_t_min, grid_pos_t_max, self._grid_interval):
            grid_x = map_x1 + (map_x2-map_x1)*((grid_pos_t-pos_t_min)/(pos_t_max-pos_t_min+1))
            dwg.add(svg.shapes.Line((grid_x, map_y1), (grid_x, map_y2+event_dim), stroke=self._grid_colour, stroke_width=self._grid_size))

        # Getting range of horizontal grid lines
        grid_pos_b_min = self._grid_interval*math.ceil(pos_b_min/self._grid_interval)
        grid_pos_b_max = self._grid_interval*math.floor(pos_b_max/self._grid_interval)
        if grid_pos_b_max%self._grid_interval == 0:
            grid_pos_b_max += 1

        # Adding horizontal grid lines
        for grid_pos_b in range(grid_pos_b_min, grid_pos_b_max, self._grid_interval):
            grid_y = map_y1 + (map_y2-map_y1)*((grid_pos_b-pos_b_min)/(pos_b_max-pos_b_min+1))
            dwg.add(svg.shapes.Line((map_x1, grid_y), (map_x2+event_dim, grid_y), stroke=self._grid_colour, stroke_width=self._grid_size))

    def _add_grid_labels(self, dwg, pos_ranges, map_xy):
        (pos_t_min, pos_t_max, pos_b_min, pos_b_max) = pos_ranges
        (map_x1, map_y1, map_x2, map_y2) = map_xy

        # Determining extra grid length if sum column and row are to be shown
        event_dim = (map_x2-map_x1)/(pos_t_max-pos_t_min+1) if self._sum_show else 0

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
        for grid_label_pos_b in range(grid_label_pos_b_min, grid_label_pos_b_max, self._grid_label_interval):
            grid_label_x = map_x1-self._border_size/2-self._grid_label_gap
            grid_label_y = map_y1 + (map_y2-map_y1)*(((grid_label_pos_b+0.5)-pos_b_min)/(pos_b_max-pos_b_min+1)) + self._grid_label_size*0.375
            dwg.add(svg.text.Text(str(grid_label_pos_b), insert=(grid_label_x,grid_label_y), style="text-anchor:end", font_size=self._grid_label_size, fill=self._grid_label_colour))

        # Adding sum labels if sum column and row are to be shown
        if self._sum_show:
            # Adding top strand sum label
            grid_label_x = map_x2 + event_dim/2 + self._grid_label_size*0.375
            grid_label_y = map_y1-self._border_size/2-self._grid_label_gap
            rot = "rotate(%i,%i,%i)" % (-90,grid_label_x,grid_label_y)
            dwg.add(svg.text.Text("Sum", insert=(grid_label_x,grid_label_y), transform=rot, style="text-anchor:start", font_size=self._grid_label_size, fill=self._grid_label_colour))

            # Adding bottom strand sum label
            grid_label_x = map_x1-self._border_size/2-self._grid_label_gap
            grid_label_y = map_y2 + event_dim/2 + self._grid_label_size*0.375
            dwg.add(svg.text.Text("Sum", insert=(grid_label_x,grid_label_y), style="text-anchor:end", font_size=self._grid_label_size, fill=self._grid_label_colour))

    def _add_events(self, dwg, pos_ranges, map_xy, freq):
        (pos_t_min, pos_t_max, pos_b_min, pos_b_max) = pos_ranges
        (map_x1, map_y1, map_x2, map_y2) = map_xy
        event_dim = (map_x2-map_x1)/(pos_t_max-pos_t_min+1)

        # Determining total events and maximum number in a cell
        max_events = get_max_events(pos_ranges, freq, self._sum_show)
        sum_events = get_sum_events(pos_ranges, freq)

        # Adding background
        self._add_background(dwg, pos_ranges, map_xy)
        
        # Adding events
        for cleavage_site_t in range(pos_t_min,pos_t_max+1):
            for cleavage_site_b in range(pos_b_min,pos_b_max+1):   
                event_x1 = map_x1 + (map_x2-map_x1)*((cleavage_site_t-pos_t_min)/(pos_t_max-pos_t_min+1))
                event_y1 = map_y1 + (map_y2-map_y1)*((cleavage_site_b-pos_b_min)/(pos_b_max-pos_b_min+1))

                (norm_count, event_pc) = get_event_stats((cleavage_site_t,cleavage_site_b), freq, max_events, sum_events)
                if (cleavage_site_t,cleavage_site_b) in freq.keys() or (self._event_label_show and self._event_label_zeros_show):
                    self._add_event(dwg, event_x1, event_y1, event_dim, norm_count, event_pc)
        
        # If enabled, showing sum row and column
        if self._sum_show:
            (freq_t, freq_b) = get_full_sequence_summed_frequency(freq, pos_ranges)

            for cleavage_site_t in range(pos_t_min,pos_t_max+1):
                event_x1 = map_x1 + (map_x2-map_x1)*((cleavage_site_t-pos_t_min)/(pos_t_max-pos_t_min+1))
                event_y1 = map_y2

                (norm_count, event_pc) = get_event_stats(cleavage_site_t, freq_t, max_events, sum_events)
                if cleavage_site_t in freq_t or (self._event_label_show and self._event_label_zeros_show):
                    self._add_event(dwg, event_x1, event_y1, event_dim, norm_count, event_pc)

            for cleavage_site_b in range(pos_b_min,pos_b_max+1):
                event_x1 = map_x2
                event_y1 = map_y1 + (map_y2-map_y1)*((cleavage_site_b-pos_b_min)/(pos_b_max-pos_b_min+1))

                (norm_count, event_pc) = get_event_stats(cleavage_site_b, freq_b, max_events, sum_events)
                if cleavage_site_b in freq_b or (self._event_label_show and self._event_label_zeros_show):
                    self._add_event(dwg, event_x1, event_y1, event_dim, norm_count, event_pc)
        
    def _add_background(self, dwg, pos_ranges, map_xy):
        (pos_t_min, pos_t_max, pos_b_min, pos_b_max) = pos_ranges
        (map_x1, map_y1, map_x2, map_y2) = map_xy        

        map_w = map_x2-map_x1
        map_h = map_y2-map_y1

        rgba = self._cmap(0)
        col = "rgb(%i,%i,%i)" % (rgba[0]*255,rgba[1]*255,rgba[2]*255)
        dwg.add(svg.shapes.Rect(insert=(map_x1,map_y1), size=(map_w,map_h), fill=col, stroke='none'))

        if self._sum_show:
            event_dim = (map_x2-map_x1)/(pos_t_max-pos_t_min+1)
            dwg.add(svg.shapes.Rect(insert=(map_x1+map_w,map_y1), size=(event_dim,map_h), fill=col, stroke='none'))
            dwg.add(svg.shapes.Rect(insert=(map_x1,map_y1+map_h), size=(map_w,event_dim), fill=col, stroke='none'))            

    def _add_event(self, dwg, event_x1, event_y1, event_dim, norm_count, event_pc):
        rgba = self._cmap(norm_count)
        col = "rgb(%i,%i,%i)" % (rgba[0]*255,rgba[1]*255,rgba[2]*255)
        
        dwg.add(svg.shapes.Rect(insert=(event_x1,event_y1), size=(event_dim,event_dim), fill=col, stroke='none'))

        # If selected, adding event percentage labels
        if self._event_label_show:              
            self._add_event_label(dwg, event_x1, event_y1, event_dim, norm_count, event_pc)

    def _add_event_label(self, dwg, event_x1, event_y1, event_dim, norm_count, event_pc):
        event_label_x = event_x1 + event_dim/2
        event_label_y = event_y1 + event_dim/2 + self._event_label_size*0.375 
        
        if self._event_label_colour == "invert":
            rgba = self._cmap(norm_count)
            col = "rgb(%i,%i,%i)" % (255-(rgba[0]*255),255-(rgba[1]*255),255-(rgba[2]*255))
        else:
            col = self._event_label_colour

        dwg.add(svg.text.Text("%.1f" % event_pc, insert=(event_label_x,event_label_y), style="text-anchor:middle", font_size=self._event_label_size, fill=col))
        
def get_event_pos_ranges(freq, round=1):
    pos_t_min = sys.maxsize
    pos_t_max = 0
    pos_b_min = sys.maxsize
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
    
    return (pos_t_min,pos_t_max,pos_b_min,pos_b_max)

def get_max_events(pos_ranges, freq, sum_show):
    (pos_t_min, pos_t_max, pos_b_min, pos_b_max) = pos_ranges

    max_events = 0
    if sum_show:
        (freq_t, freq_b) = get_full_sequence_summed_frequency(freq, pos_ranges)

        for t in freq_t.keys():
            if t >= pos_t_min and t <= pos_t_max:
                max_events = max(max_events,freq_t.get(t))

        for b in freq_b.keys():
            if b >= pos_b_min and b <= pos_b_max:
                max_events = max(max_events,freq_b.get(b))

    else:
        for (t,b) in freq.keys():
            if t >= pos_t_min and t <= pos_t_max and b >= pos_b_min and b <= pos_b_max:
                max_events = max(max_events,freq.get((t,b)))

    return max_events

def get_sum_events(pos_ranges, freq):
    (pos_t_min, pos_t_max, pos_b_min, pos_b_max) = pos_ranges

    sum_events = 0
    for (t,b) in freq.keys():
            if t >= pos_t_min and t <= pos_t_max and b >= pos_b_min and b <= pos_b_max:
                sum_events = sum_events + freq.get((t,b))

    return sum_events

def get_full_sequence_summed_frequency(freq_full, pos_ranges):
    (pos_t_min, pos_t_max, pos_b_min, pos_b_max) = pos_ranges

    freq_t = {}
    freq_b = {}

    for (t, b) in freq_full.keys():
        freq = freq_full.get((t, b))

        # Adding this frequency to the relevant elements of freq_t and freq_b
        if t >= pos_t_min and t <= pos_t_max and b >= pos_b_min and b <= pos_b_max:
            freq_t[t] = freq_t[t] + freq if t in freq_t else freq
            freq_b[b] = freq_b[b] + freq if b in freq_b else freq

    return (freq_t, freq_b)

def get_event_stats(key, freq, max_events, sum_events):
    if key in freq:
        norm_count = freq.get(key)/max_events
        event_pc = 100*freq.get(key)/sum_events  
    else:
        norm_count = 0
        event_pc = 0

    return (norm_count, event_pc)

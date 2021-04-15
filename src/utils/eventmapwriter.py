# TODO: Have colour normalisation to visible position range only

import datetime as dt
import math
import os
import svgwrite as svg

from enum import Enum
from matplotlib import cm
from utils import abstractmapwriter as amw
from utils import reportutils as ru


class DNA_MODE(Enum):
    LINE = 'line'
    NONE = 'none'
    SEQUENCE = 'seq'

    def __str__(self):
        return str(self.value)


class EventMapWriter(amw.AbstractMapWriter):
    ## CONSTRUCTOR

    def __init__(self, im_dims=(800,200), rel_pos=(0.3,0.6,0.05,0.8), dna_opts=(DNA_MODE.LINE,2,"black"), end_label_opts=(True,20,"black",10), grid_opts=(True,1,"gray",100), 
    grid_label_opts=(True,12,"gray",500,10), colourbar_opts=(True,0.91,0.02,1), colourbar_label_opts=(True,12,"black"), event_opts=(0.2,2,"cool",0.2,True,1)):
        self._im_w = im_dims[0]
        self._im_h = im_dims[1]

        self._map_rel_top = rel_pos[0]
        self._map_rel_height = rel_pos[1]
        self._map_rel_left = rel_pos[2]
        self._map_rel_width = rel_pos[3]
        
        self._dna_mode = dna_opts[0]
        self._dna_size = dna_opts[1]
        self._dna_colour = dna_opts[2]

        self._end_label_show = end_label_opts[0]
        self._end_label_size = end_label_opts[1]
        self._end_label_colour = end_label_opts[2]
        self._end_label_gap = end_label_opts[3]
        
        self._grid_show = grid_opts[0]
        self._grid_size = grid_opts[1]
        self._grid_colour = grid_opts[2]
        self._grid_interval = grid_opts[3]

        self._grid_label_show = grid_label_opts[0]
        self._grid_label_size = grid_label_opts[1]
        self._grid_label_colour = grid_label_opts[2]
        self._grid_label_interval = grid_label_opts[3]
        self._grid_label_gap = grid_label_opts[4]

        self._colourbar_show = colourbar_opts[0]
        self._colourbar_rel_left = colourbar_opts[1]
        self._colourbar_rel_width = colourbar_opts[2]
        self._colourbar_border_size = colourbar_opts[3]

        self._colourbar_label_show = colourbar_label_opts[0]
        self._colourbar_label_size = colourbar_label_opts[1]
        self._colourbar_label_colour = colourbar_label_opts[2]
        
        self._event_min_size = event_opts[0]
        self._event_max_size = event_opts[1]
        self._event_colourmap = event_opts[2]
        self._event_opacity = event_opts[3]
        self._event_fill_range = event_opts[4]
        self._event_stack_order = event_opts[5]


    ## GETTERS AND SETTERS

    def set_dna_mode(self, dna_mode):
        self._dna_mode = dna_mode


    ## PUBLIC METHODS

    def write_map(self, out_path, freq, ref=None, pos_range=None, append_dt=False):
        if self._event_stack_order == 1:
            freq = ru.sort_results(freq,True)
        elif self._event_stack_order == 2:
            freq = ru.sort_results(freq,False)
            
        (pos_min, pos_max) = amw.get_single_pos_range(freq, ref, pos_range)

        # Creating output SVG document
        root_name = os.path.splitext(out_path)[0]
        datetime_str = dt.datetime.now().strftime("_%Y-%m-%d_%H-%M-%S") if append_dt else ""
        outname = root_name+ datetime_str + '.' + 'svg'
        dwg = svg.Drawing(outname, size = ("%spx" % self._im_w, "%spx" % self._im_h))

        # Defining limits of DNA strands
        map_x1 = self._im_w*self._map_rel_left
        map_y1 = self._im_h*self._map_rel_top
        map_x2 = map_x1 + self._im_w*self._map_rel_width
        map_y2 = map_y1 + self._im_h*self._map_rel_height
        map_xy = (map_x1, map_y1, map_x2, map_y2)

        if self._grid_show:
            self._add_grid_lines(dwg, pos_min, pos_max, map_xy)

        if self._grid_label_show:
            self._add_grid_labels(dwg, pos_min, pos_max, map_xy)
        
        if self._end_label_show:
            self._add_end_labels(dwg, map_xy)

        if self._colourbar_show:
            self._add_colourbar(dwg, freq)

        if len(freq) > 0:
            # Getting reference length (or estimating)
            if ref is None:
                pos_ranges = amw.get_event_pos_range(freq, round=1)
                ref_len = max(pos_ranges[1],pos_ranges[3])
            else:
                ref_len = len(ref)

            self._add_event_lines(dwg, pos_min, pos_max, map_xy, freq, ref_len)
            
        self._add_dna(dwg, pos_min, pos_max, map_xy, ref=ref)

        # Writing SVG to file
        dwg.save()

    def _add_grid_lines(self, dwg, pos_min, pos_max, map_xy):
        (map_x1, map_y1, map_x2, map_y2) = map_xy

        # Adding horizontal grid line
        grid_yc = (map_y1+map_y2)/2
        dwg.add(svg.shapes.Line((map_x1, grid_yc), (map_x2, grid_yc), stroke=self._grid_colour, stroke_width=self._grid_size))

        # Getting range of vertical grid lines
        grid_pos_min = self._grid_interval*math.ceil(pos_min/self._grid_interval)
        grid_pos_max = self._grid_interval*math.floor(pos_max/self._grid_interval)
        if grid_pos_max%self._grid_interval == 0:
            grid_pos_max += 1

        # Adding vertical grid lines
        for grid_pos in range(grid_pos_min, grid_pos_max, self._grid_interval):
            grid_x = map_x1 + (map_x2-map_x1)*((grid_pos-pos_min)/(pos_max-pos_min))
            grid_y1 = map_y1+self._dna_size/2
            grid_y2 = map_y2-self._dna_size/2
            dwg.add(svg.shapes.Line((grid_x, grid_y1), (grid_x, grid_y2), stroke=self._grid_colour, stroke_width=self._grid_size))

    def _add_grid_labels(self, dwg, pos_min, pos_max, map_xy):
        (map_x1, map_y1, map_x2, map_y2) = map_xy

        # Getting range of grid label positions
        grid_label_min = self._grid_label_interval*math.ceil(pos_min/self._grid_label_interval)
        grid_label_max = self._grid_label_interval*math.floor(pos_max/self._grid_label_interval)
        if grid_label_max%self._grid_label_interval == 0:
            grid_label_max += 1

        # Adding grid labels
        for grid_label_pos in range(grid_label_min, grid_label_max, self._grid_label_interval):
            grid_label_x = map_x1 + (map_x2-map_x1)*((grid_label_pos-pos_min)/(pos_max-pos_min)) + self._grid_label_size*0.375
            grid_label_y = map_y1-self._dna_size/2-self._grid_label_gap
            rot = "rotate(%i,%i,%i)" % (-90,grid_label_x,grid_label_y)
            dwg.add(svg.text.Text(str(grid_label_pos), insert=(grid_label_x,grid_label_y), transform=rot, style="text-anchor:start", font_size=self._grid_label_size, fill=self._grid_label_colour))

    def _add_dna(self, dwg, pos_min, pos_max, map_xy, ref=None):
        (map_x1, map_y1, map_x2, map_y2) = map_xy
        
        if self._dna_mode is DNA_MODE.SEQUENCE and ref is None:
            print("WARNING: No reference path provided.  DNA rendered as solid lines.")
            self._dna_mode = DNA_MODE.LINE

        if self._dna_mode is DNA_MODE.SEQUENCE:
            # Draw DNA as text sequence
            dna_pos_interval = (map_x2-map_x1)/(pos_max-pos_min)
            ref_rc = ref.reverse_complement()
            for dna_pos in range(pos_min,pos_max+1):
                dna_seq_x = map_x1+dna_pos_interval*(dna_pos-pos_min)
                dna_seq_y1 = map_y1+self._dna_size*0.375
                dna_seq_y2 = map_y2+self._dna_size*0.375
                dwg.add(svg.text.Text(str(ref[dna_pos]), insert=(dna_seq_x, dna_seq_y1), style="text-anchor:middle", font_size=self._dna_size, fill=self._dna_colour))
                dwg.add(svg.text.Text(str(ref_rc[dna_pos]), insert=(dna_seq_x, dna_seq_y2), style="text-anchor:middle", font_size=self._dna_size, fill=self._dna_colour))

        elif self._dna_mode is DNA_MODE.LINE:
            # Draw DNA as lines
            dwg.add(svg.shapes.Line((map_x1, map_y1), (map_x2, map_y1), stroke=self._dna_colour, stroke_width=self._dna_size, style="stroke-linecap:square"))
            dwg.add(svg.shapes.Line((map_x1, map_y2), (map_x2, map_y2), stroke=self._dna_colour, stroke_width=self._dna_size, style="stroke-linecap:square"))

    def _add_end_labels(self, dwg, map_xy):
        (map_x1, map_y1, map_x2, map_y2) = map_xy

        end_label_x1 = map_x1-self._end_label_gap
        end_label_x2 = map_x2+self._end_label_gap
        end_label_y1 = map_y1+self._end_label_size*0.375
        end_label_y2 = map_y2+self._end_label_size*0.375

        dwg.add(svg.text.Text("5'", insert=(end_label_x1, end_label_y1), style="text-anchor:end", font_size=self._end_label_size, fill=self._end_label_colour))
        dwg.add(svg.text.Text("3'", insert=(end_label_x2, end_label_y1), style="text-anchor:start", font_size=self._end_label_size, fill=self._end_label_colour))
        dwg.add(svg.text.Text("3'", insert=(end_label_x1, end_label_y2), style="text-anchor:end", font_size=self._end_label_size, fill=self._end_label_colour))
        dwg.add(svg.text.Text("5'", insert=(end_label_x2, end_label_y2), style="text-anchor:start", font_size=self._end_label_size, fill=self._end_label_colour))

    def _add_colourbar(self, dwg, freq):
        colourbar_x1 = self._im_w*self._colourbar_rel_left
        colourbar_y1 = self._im_h*self._map_rel_top
        colourbar_x2 = colourbar_x1 + self._im_w*self._colourbar_rel_width
        colourbar_y2 = colourbar_y1 + self._im_h*self._map_rel_height
        colourbar_w = self._im_w*self._colourbar_rel_width
        colourbar_h = self._im_h*self._map_rel_height

        # Adding colourbar lines
        cmap = cm.get_cmap(self._event_colourmap)    

        line_w = colourbar_h/256+0.1 # A bit extra width to prevent gaps appearing between lines
        line_offs = colourbar_h/255
        for i in range(0,256):
            y = colourbar_y2-i*line_offs
            rgba = cmap(i/256)
            col = "rgb(%i,%i,%i)" % (rgba[0]*255,rgba[1]*255,rgba[2]*255)
            dwg.add(svg.shapes.Line((colourbar_x1, y), (colourbar_x2, y), stroke=col, stroke_width=line_w, style="stroke-linecap:square"))

        # Adding border around colourbar
        dwg.add(svg.shapes.Rect(insert=(colourbar_x1,colourbar_y1), size=(colourbar_w,colourbar_h), stroke="black", fill="none", stroke_width=self._colourbar_border_size, style="stroke-linecap:square"))

        if self._colourbar_label_show:
            total = sum(freq.values())
            max_events = max(freq.values())
            min_events = min(freq.values())

            if self._event_fill_range:
                min_text = "%.1f%%" % ((min_events/total)*100)
                max_text = "%.1f%%" % ((max_events/total)*100)
            else:
                min_text = "0%"
                max_text = "100%"

            end_label_x = colourbar_x2+self._end_label_gap
            end_label_y1 = colourbar_y1+self._colourbar_label_size*0.375
            end_label_y2 = colourbar_y2+self._colourbar_label_size*0.375

            dwg.add(svg.text.Text(max_text, insert=(end_label_x, end_label_y1), style="text-anchor:start", font_size=self._colourbar_label_size, fill=self._colourbar_label_colour))
            dwg.add(svg.text.Text(min_text, insert=(end_label_x, end_label_y2), style="text-anchor:start", font_size=self._colourbar_label_size, fill=self._colourbar_label_colour))

    def _add_event_lines(self, dwg, pos_min, pos_max, map_xy, freq, ref_len):
        # Adding cleavage event lines
        total = sum(freq.values())
        max_events = max(freq.values())
        min_events = min(freq.values())
        diff_events = max_events-min_events

        # Preventing divide by zero errors
        if diff_events == 0:
            diff_events = 1

        for (cleavage_site_t, cleavage_site_b, split) in freq.keys():   
            # Checking this event is within the rendered range
            if cleavage_site_t < pos_min or cleavage_site_t > pos_max or cleavage_site_b < pos_min or cleavage_site_b > pos_max:
                continue;
                
            if self._event_fill_range:
                norm_count = (freq.get((cleavage_site_t, cleavage_site_b, split))-min_events)/diff_events
            else:
                norm_count = freq.get((cleavage_site_t, cleavage_site_b, split))/total
            
            if split:
                self._add_split_line(dwg, cleavage_site_t, cleavage_site_b, pos_min, pos_max, map_xy, norm_count, ref_len)
            else:
                self._add_continuous_line(dwg, cleavage_site_t, cleavage_site_b, pos_min, pos_max, map_xy, norm_count)
                     

    def _add_continuous_line(self, dwg, cleavage_site_t, cleavage_site_b, pos_min, pos_max, map_xy, norm_count):    
        (map_x1, map_y1, map_x2, map_y2) = map_xy
        cmap = cm.get_cmap(self._event_colourmap)

        # Adding line (adding width 1 to ensure everything is visible)
        event_width = (self._event_max_size-self._event_min_size)*norm_count+self._event_min_size
        
        event_t_x = map_x1 + (map_x2-map_x1)*((cleavage_site_t-pos_min)/(pos_max-pos_min))
        event_t_y = map_y1+self._dna_size/2

        event_b_x = map_x1 + (map_x2-map_x1)*((cleavage_site_b-pos_min)/(pos_max-pos_min))
        event_b_y = map_y2-self._dna_size/2

        rgba = cmap(norm_count)
        col = "rgb(%i,%i,%i)" % (rgba[0]*255,rgba[1]*255,rgba[2]*255)

        dwg.add(svg.shapes.Line((event_t_x, event_t_y), (event_b_x, event_b_y), stroke=col, stroke_width=event_width, style="stroke-linecap:square;stroke-opacity:%f" % self._event_opacity))

    def _add_split_line(self, dwg, cleavage_site_t, cleavage_site_b, pos_min, pos_max, map_xy, norm_count, ref_len):  
        (map_x1, map_y1, map_x2, map_y2) = map_xy
        cmap = cm.get_cmap(self._event_colourmap)

        # Adding line (adding width 1 to ensure everything is visible)
        event_width = (self._event_max_size-self._event_min_size)*norm_count+self._event_min_size
                
        # Processing is dependent on which site is lower
        if cleavage_site_t < cleavage_site_b:
            site_sep = ref_len - (cleavage_site_b-cleavage_site_t)

            event_t_x1 = map_x1
            event_t_y1 = map_y1 + (map_y2-map_y1)*((cleavage_site_t-pos_min)/site_sep)
            event_t_x2 = map_x1 + (map_x2-map_x1)*((cleavage_site_t-pos_min)/(pos_max-pos_min))
            event_t_y2 = map_y1 + self._dna_size/2

            event_b_x1 = map_x1 + (map_x2-map_x1)*((cleavage_site_b-pos_min)/(pos_max-pos_min))
            event_b_y1 = map_y2 - self._dna_size/2
            event_b_x2 = map_x2
            event_b_y2 = map_y2 - (map_y2-map_y1)*((pos_max-cleavage_site_b)/site_sep)

        else:
            site_sep = ref_len - (cleavage_site_t-cleavage_site_b)

            event_t_x1 = map_x1 + (map_x2-map_x1)*((cleavage_site_t-pos_min)/(pos_max-pos_min))
            event_t_y1 = map_y1 + self._dna_size/2
            event_t_x2 = map_x2
            event_t_y2 = map_y1 + (map_y2-map_y1)*((pos_max-cleavage_site_t)/site_sep)

            event_b_x1 = map_x1
            event_b_y1 = map_y2 - (map_y2-map_y1)*((cleavage_site_b-pos_min)/site_sep)
            event_b_x2 = map_x1 + (map_x2-map_x1)*((cleavage_site_b-pos_min)/(pos_max-pos_min))
            event_b_y2 = map_y2 - self._dna_size/2

        rgba = cmap(norm_count)
        col = "rgb(%i,%i,%i)" % (rgba[0]*255,rgba[1]*255,rgba[2]*255)

        dwg.add(svg.shapes.Line((event_t_x1, event_t_y1), (event_t_x2, event_t_y2), stroke=col, stroke_width=event_width, style="stroke-linecap:square;stroke-opacity:%f" % self._event_opacity))
        dwg.add(svg.shapes.Line((event_b_x1, event_b_y1), (event_b_x2, event_b_y2), stroke=col, stroke_width=event_width, style="stroke-linecap:square;stroke-opacity:%f" % self._event_opacity))

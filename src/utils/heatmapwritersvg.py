from matplotlib import cm
from utils import abstractmapwriter as amw

import datetime as dt
import math
import os
import svgwrite as svg


class HeatMapWriterSVG(amw.AbstractMapWriter):
    ## CONSTRUCTOR

    def __init__(self, im_dim=800, font="Arial", rel_pos=(0.1,0.1,0.8), border_opts=(True,1,"black"), axis_label_opts=(True,16,"gray",50), grid_opts=(True,1,"gray",1), grid_label_opts=(True,12,"gray",10,10), event_colourmap="plasma", event_label_opts=(True,10,"invert",1,True), sum_show=True):
        self._im_dim = im_dim

        self._font = font
        
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
        self._event_label_number_format = "%%.%if" % event_label_opts[3]
        self._event_label_zeros_show = event_label_opts[4]

        self._sum_show = sum_show

    def write_map(self, out_path, freq, ref, pos_range, append_dt):
        # Getting pos ranges to plot based on available information
        pos_range = amw.get_double_pos_range(freq, ref, pos_range)
        if pos_range is None:
            return

        # Creating output SVG document
        root_name = os.path.splitext(out_path)[0]
        datetime_str = dt.datetime.now().strftime("_%Y-%m-%d_%H-%M-%S") if append_dt else ""
        outname = root_name+ datetime_str + '.' + 'svg'
        
        dwg = svg.Drawing(outname, size = ("%spx" % self._im_dim, "%spx" % self._im_dim), **{'shape-rendering':'crispEdges'})
        
        # Defining limits of heat map
        pos_t_range = pos_range[1]-pos_range[0]
        pos_b_range = pos_range[3]-pos_range[2]
        max_range = max(pos_t_range, pos_b_range)
        pos_t_rel_size = pos_t_range/max_range
        pos_b_rel_size = pos_b_range/max_range

        map_x1 = self._im_dim*self._map_rel_left
        map_x2 = map_x1 + self._im_dim*self._map_rel_size*pos_t_rel_size
        map_y1 = self._im_dim*self._map_rel_top
        map_y2 = map_y1 + self._im_dim*self._map_rel_size*pos_b_rel_size
        map_xy = (map_x1, map_y1, map_x2, map_y2)

        self._add_events(dwg, pos_range, map_xy, freq)

        if self._axis_label_show:
            self._add_axis_labels(dwg, map_xy)

        if self._grid_show:
            self._add_grid_lines(dwg, pos_range, map_xy)

        if self._grid_label_show:
            self._add_grid_labels(dwg, pos_range, map_xy)
                
        if self._border_show:
            self._add_border(dwg, pos_range, map_xy)

        # Writing SVG to file
        dwg.save()

    def _add_axis_labels(self, dwg, map_xy):
        (map_x1, map_y1, map_x2, map_y2) = map_xy

        # Adding top-strand label
        axis_label_x = (map_x2-map_x1)/2 + map_x1
        axis_label_y = map_y1 - self._axis_label_gap
        dwg.add(svg.text.Text("Top strand", insert=(axis_label_x,axis_label_y), style=f"text-anchor:middle;font-family:{self._font}", font_size=self._axis_label_size, fill=self._axis_label_colour))

        # Adding bottom-strand label
        axis_label_x = map_x1 - self._axis_label_gap
        axis_label_y = (map_y2-map_y1)/2 + map_y1
        rot = "rotate(%i,%i,%i)" % (-90,axis_label_x,axis_label_y)
        dwg.add(svg.text.Text("Bottom strand", insert=(axis_label_x,axis_label_y), transform=rot, style=f"text-anchor:middle;font-family:{self._font}", font_size=self._axis_label_size, fill=self._axis_label_colour))

    def _add_border(self, dwg, pos_range, map_xy):
        (pos_t_min, pos_t_max, pos_b_min, pos_b_max) = pos_range
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
        
    def _add_grid_lines(self, dwg, pos_range, map_xy):
        (pos_t_min, pos_t_max, pos_b_min, pos_b_max) = pos_range
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

    def _add_grid_labels(self, dwg, pos_range, map_xy):
        (pos_t_min, pos_t_max, pos_b_min, pos_b_max) = pos_range
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
            grid_label_x = map_x1 + (map_x2-map_x1)*(((grid_label_pos_t+0.5)-pos_t_min)/(pos_t_max-pos_t_min+1))
            grid_label_y = map_y1-self._border_size/2-self._grid_label_gap
            rot = "rotate(%i,%i,%i)" % (-90,grid_label_x,grid_label_y)
            dwg.add(svg.text.Text(str(grid_label_pos_t), insert=(grid_label_x,grid_label_y), transform=rot, style=f"text-anchor:start;font-family:{self._font};dominant-baseline:central", font_size=self._grid_label_size, fill=self._grid_label_colour))

        # Getting range of horizontal grid label positions
        grid_label_pos_b_min = self._grid_label_interval*math.ceil(pos_b_min/self._grid_label_interval)
        grid_label_pos_b_max = self._grid_label_interval*math.floor(pos_b_max/self._grid_label_interval)
        if grid_label_pos_b_max%self._grid_label_interval == 0:
            grid_label_pos_b_max += 1

        # Adding horizontal line grid labels
        for grid_label_pos_b in range(grid_label_pos_b_min, grid_label_pos_b_max, self._grid_label_interval):
            grid_label_x = map_x1-self._border_size/2-self._grid_label_gap
            grid_label_y = map_y1 + (map_y2-map_y1)*(((grid_label_pos_b+0.5)-pos_b_min)/(pos_b_max-pos_b_min+1))
            dwg.add(svg.text.Text(str(grid_label_pos_b), insert=(grid_label_x,grid_label_y), style=f"text-anchor:end;font-family:{self._font};dominant-baseline:central", font_size=self._grid_label_size, fill=self._grid_label_colour))

        # Adding sum labels if sum column and row are to be shown
        if self._sum_show:
            # Adding top strand sum label
            grid_label_x = map_x2 + event_dim/2
            grid_label_y = map_y1-self._border_size/2-self._grid_label_gap
            rot = "rotate(%i,%i,%i)" % (-90,grid_label_x,grid_label_y)
            dwg.add(svg.text.Text("Sum", insert=(grid_label_x,grid_label_y), transform=rot, style=f"text-anchor:start;font-family:{self._font};dominant-baseline:central", font_size=self._grid_label_size, fill=self._grid_label_colour))

            # Adding bottom strand sum label
            grid_label_x = map_x1-self._border_size/2-self._grid_label_gap
            grid_label_y = map_y2 + event_dim/2
            dwg.add(svg.text.Text("Sum", insert=(grid_label_x,grid_label_y), style=f"text-anchor:end;font-family:{self._font};dominant-baseline:central", font_size=self._grid_label_size, fill=self._grid_label_colour))

    def _add_events(self, dwg, pos_range, map_xy, freq):
        (pos_t_min, pos_t_max, pos_b_min, pos_b_max) = pos_range
        (map_x1, map_y1, map_x2, map_y2) = map_xy
        event_dim = (map_x2-map_x1)/(pos_t_max-pos_t_min+1)

        # Determining total events and maximum number in a cell
        max_events = amw.get_max_events(pos_range, freq, self._sum_show)
        sum_events = amw.get_sum_events(pos_range, freq)

        # Adding background
        self._add_background(dwg, pos_range, map_xy)
        
        # Adding events
        for pos_t in range(pos_t_min,pos_t_max+1):
            for pos_b in range(pos_b_min,pos_b_max+1):  
                event_x1 = map_x1 + (map_x2-map_x1)*((pos_t-pos_t_min)/(pos_t_max-pos_t_min+1))
                event_y1 = map_y1 + (map_y2-map_y1)*((pos_b-pos_b_min)/(pos_b_max-pos_b_min+1))

                norm_count = amw.get_event_norm_count((pos_t,pos_b,True), freq, max_events) + amw.get_event_norm_count((pos_t,pos_b,False), freq, max_events)
                event_pc = amw.get_event_pc((pos_t,pos_b,True), freq, sum_events) + amw.get_event_pc((pos_t,pos_b,False), freq, sum_events)
                if (pos_t,pos_b,True) in freq.keys() or (pos_t,pos_b,False) in freq.keys():                    
                    self._add_event(dwg, event_x1, event_y1, event_dim, norm_count)

                if self._event_label_show and (((pos_t,pos_b,True) in freq.keys() or ((pos_t,pos_b,False) in freq.keys())) or self._event_label_zeros_show):
                    self._add_event_label(dwg, event_x1, event_y1, event_dim, norm_count, event_pc)
                    
        # If enabled, showing sum row and column
        if self._sum_show:
            (freq_t, freq_b) = amw.get_full_sequence_summed_frequency(freq, pos_range)

            for pos_t in range(pos_t_min,pos_t_max+1):
                event_x1 = map_x1 + (map_x2-map_x1)*((pos_t-pos_t_min)/(pos_t_max-pos_t_min+1))
                event_y1 = map_y2

                norm_count = amw.get_event_norm_count(pos_t, freq_t, max_events)
                event_pc = amw.get_event_pc(pos_t, freq_t, sum_events)
                if pos_t in freq_t or self._event_label_zeros_show:                    
                    self._add_event(dwg, event_x1, event_y1, event_dim, norm_count)

                if self._event_label_show and (pos_t in freq_t or self._event_label_zeros_show):
                    self._add_event_label(dwg, event_x1, event_y1, event_dim, norm_count, event_pc)

            for pos_b in range(pos_b_min,pos_b_max+1):
                event_x1 = map_x2
                event_y1 = map_y1 + (map_y2-map_y1)*((pos_b-pos_b_min)/(pos_b_max-pos_b_min+1))

                norm_count = amw.get_event_norm_count(pos_b, freq_b, max_events)
                event_pc = amw.get_event_pc(pos_b, freq_b, sum_events)
                if pos_b in freq_b or self._event_label_zeros_show:
                    self._add_event(dwg, event_x1, event_y1, event_dim, norm_count)

                if self._event_label_show and (pos_b in freq_b or self._event_label_zeros_show):
                    self._add_event_label(dwg, event_x1, event_y1, event_dim, norm_count, event_pc)
        
    def _add_background(self, dwg, pos_range, map_xy):
        (pos_t_min, pos_t_max, pos_b_min, pos_b_max) = pos_range
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

    def _add_event(self, dwg, event_x1, event_y1, event_dim, norm_count):
        rgba = self._cmap(norm_count)
        col = "rgb(%i,%i,%i)" % (rgba[0]*255,rgba[1]*255,rgba[2]*255)
        
        dwg.add(svg.shapes.Rect(insert=(event_x1,event_y1), size=(event_dim,event_dim), fill=col, stroke='none'))

    def _add_event_label(self, dwg, event_x1, event_y1, event_dim, norm_count, event_pc):
        event_label_x = event_x1 + event_dim/2
        event_label_y = event_y1 + event_dim/2
        
        if self._event_label_colour == "invert":
            rgba = self._cmap(norm_count)
            col = "rgb(%i,%i,%i)" % (255-(rgba[0]*255),255-(rgba[1]*255),255-(rgba[2]*255))
        else:
            col = self._event_label_colour

        dwg.add(svg.text.Text(self._event_label_number_format % event_pc, insert=(event_label_x,event_label_y), style=f"text-anchor:middle;font-family:{self._font};dominant-baseline:central", font_size=self._event_label_size, fill=col))
    

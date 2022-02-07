import datetime as dt
import math
import os
import svgwrite as svg

from enum import Enum
from enums.eventtype import Eventtype
from matplotlib import cm
from rich import print
from utils import abstractwriter as aw
from utils import reportutils as ru


class DNA_MODE(Enum):
    LINE = 'line'
    NONE = 'none'
    SEQUENCE = 'seq'

    def __str__(self):
        return str(self.value)


class HPOS(Enum):
    LEFT = 'left'
    RIGHT = 'right'

    def __str__(self):
        return str(self.value)


class VPOS(Enum):
    CENTRE = 'centre'
    INSIDE = 'inside'

    def __str__(self):
        return str(self.value)


class StrandLinkagePlotWriter(aw.AbstractWriter):
    ## CONSTRUCTOR

    def __init__(self,
                 im_dims=(800, 300),
                 font="Arial",
                 rel_pos=(0.42, 0.3, 0.07, 0.78),
                 dna_opts=(DNA_MODE.LINE, 2.5, "black", 0.04),
                 end_label_opts=(True, 20, "black", 0.01, VPOS.CENTRE),
                 grid_opts=(True, 1, "lightgray", 100),
                 grid_label_opts=(True, 12, "gray", 500, 0.04),
                 cbar_opts=(True, 0.91, 0.02, 1),
                 cbar_label_opts=(True, 12, "gray", 25, 0.02),
                 event_opts=(0.5, 2, "cool", 0, 100, True, 0.4, 1),
                 hist_opts=(True, 0, 100, 2, "darkgray", 0.16, 0.07, 20, 0),
                 hist_label_opts=(True, 12, "gray", 25, 0.01, HPOS.LEFT, True),
                 hist_grid_opts=(True, 1, "lightgray", 25),
                 hist_name_opts=(True, 12, "gray", 0.01, HPOS.RIGHT),
                 splithist_opts=(False, 0, 100, "darkgray", 0.16),
                 splithist_label_opts=50,
                 splithist_grid_opts=50):

        self._im_w = im_dims[0]
        self._im_h = im_dims[1]

        self._font = font

        self._map_rel_top = rel_pos[0]
        self._map_rel_height = rel_pos[1]
        self._map_rel_left = rel_pos[2]
        self._map_rel_width = rel_pos[3]

        self._dna_mode = dna_opts[0]
        self._dna_size = dna_opts[1]
        self._dna_colour = dna_opts[2]
        self._dna_rel_gap = dna_opts[3]

        self._end_label_show = end_label_opts[0]
        self._end_label_size = end_label_opts[1]
        self._end_label_colour = end_label_opts[2]
        self._end_label_rel_gap = end_label_opts[3]
        self._end_label_position = end_label_opts[4]

        self._grid_show = grid_opts[0]
        self._grid_size = grid_opts[1]
        self._grid_colour = grid_opts[2]
        self._grid_interval = grid_opts[3]

        self._grid_label_show = grid_label_opts[0]
        self._grid_label_size = grid_label_opts[1]
        self._grid_label_colour = grid_label_opts[2]
        self._grid_label_interval = grid_label_opts[3]
        self._grid_label_rel_gap = grid_label_opts[4]

        self._cbar_show = cbar_opts[0]
        self._cbar_rel_left = cbar_opts[1]
        self._cbar_rel_width = cbar_opts[2]
        self._cbar_border_size = float(cbar_opts[3])

        self._cbar_label_show = cbar_label_opts[0]
        self._cbar_label_size = cbar_label_opts[1]
        self._cbar_label_colour = cbar_label_opts[2]
        self._cbar_label_interval = cbar_label_opts[3]
        self._cbar_label_rel_gap = cbar_label_opts[4]

        self._event_min_size = event_opts[0]
        self._event_max_size = event_opts[1]
        self._event_colourmap = event_opts[2]
        self._event_min_range = event_opts[3]
        self._event_max_range = event_opts[4]
        self._event_outside_range_show = event_opts[5]
        self._event_opacity = event_opts[6]
        self._event_stack_order = event_opts[7]

        self._hist_show = hist_opts[0]
        self._hist_min_range = hist_opts[1]
        self._hist_max_range = hist_opts[2]
        self._hist_bin_width = hist_opts[3]
        self._hist_colour = hist_opts[4]
        self._hist_rel_height = hist_opts[5]
        self._hist_rel_gap = hist_opts[6]
        self._hist_pc_bar_gap = hist_opts[7]
        self._hist_overhang = hist_opts[8]

        self._hist_label_show = hist_label_opts[0]
        self._hist_label_size = hist_label_opts[1]
        self._hist_label_colour = hist_label_opts[2]
        self._hist_label_interval = hist_label_opts[3]
        self._hist_label_rel_gap = hist_label_opts[4]
        self._hist_label_position = hist_label_opts[5]
        self._hist_label_zero_show = hist_label_opts[6]

        self._hist_grid_show = hist_grid_opts[0]
        self._hist_grid_size = hist_grid_opts[1]
        self._hist_grid_colour = hist_grid_opts[2]
        self._hist_grid_interval = hist_grid_opts[3]

        self._hist_name_show = hist_name_opts[0]
        self._hist_name_size = hist_name_opts[1]
        self._hist_name_colour = hist_name_opts[2]
        self._hist_name_rel_gap = hist_name_opts[3]
        self._hist_name_position = hist_name_opts[4]

        self._splithist_show = splithist_opts[0]
        self._splithist_min_range = splithist_opts[1]
        self._splithist_max_range = splithist_opts[2]
        self._splithist_rel_height = splithist_opts[3]

        self._splithist_label_interval = splithist_label_opts

        self._splithist_grid_interval = splithist_grid_opts

    ## GETTERS AND SETTERS

    def set_dna_mode(self, dna_mode):
        self._dna_mode = dna_mode

    ## PUBLIC METHODS

    def write(self,
                  out_path,
                  freq,
                  ref=None,
                  pos_range=None,
                  append_dt=False):
        if self._event_stack_order == 1:
            freq = ru.sort_results(freq, True)
        elif self._event_stack_order == 2:
            freq = ru.sort_results(freq, False)

        (pos_min, pos_max) = aw.get_single_pos_range(freq, ref, pos_range)

        # Creating output SVG document
        root_name = os.path.splitext(out_path)[0]
        datetime_str = dt.datetime.now().strftime(
            "_%Y-%m-%d_%H-%M-%S") if append_dt else ""
        outname = root_name + datetime_str + '.' + 'svg'
        dwg = svg.Drawing(outname,
                          size=("%spx" % self._im_w, "%spx" % self._im_h))

        # Defining limits of DNA strands
        map_x1 = self._im_w * self._map_rel_left
        map_y1 = self._im_h * self._map_rel_top
        map_x2 = map_x1 + self._im_w * self._map_rel_width
        map_y2 = map_y1 + self._im_h * self._map_rel_height
        map_xy = (map_x1, map_y1, map_x2, map_y2)

        if self._grid_show:
            grid_y1 = map_y1 + self._dna_size / 2 + self._dna_rel_gap*self._im_h
            grid_y2 = map_y2 - self._dna_size / 2 - self._dna_rel_gap*self._im_h
            grid_xy = (map_x1,grid_y1,map_x2,grid_y2)
            self._add_grid_lines(dwg, pos_min, pos_max, grid_xy)

        if self._grid_label_show:
            self._add_grid_labels(dwg, pos_min, pos_max, map_xy)

        if self._end_label_show:
            self._add_end_labels(dwg, map_xy)

        if self._cbar_show:
            self._add_cbar(dwg, map_xy, freq)

        n_events = sum(freq.values())

        if self._hist_show:
            (freq_t, freq_b) = ru.get_position_frequency(freq)        
            self._add_sum_hist(dwg, pos_min, pos_max, map_xy, freq_t, n_events, True)
            self._add_sum_hist(dwg, pos_min, pos_max, map_xy, freq_b, n_events, False)

        if self._splithist_show:
            (freq_t_5, freq_b_5) = ru.get_position_frequency(freq, Eventtype.FIVE_P)
            (freq_t_3, freq_b_3) = ru.get_position_frequency(freq, Eventtype.THREE_P)
            (freq_t_b, freq_b_b) = ru.get_position_frequency(freq, Eventtype.BLUNT)        
            self._add_splithists(dwg, pos_min, pos_max, map_xy, freq_t_5, freq_t_3, freq_t_b, n_events, True)
            self._add_splithists(dwg, pos_min, pos_max, map_xy, freq_b_5, freq_b_3, freq_b_b, n_events, False)

        self._add_event_lines(dwg, pos_min, pos_max, map_xy, freq, ref)

        self._add_dna(dwg, pos_min, pos_max, map_xy, ref=ref)

        # Writing SVG to file
        dwg.save()

    def _add_grid_lines(self, dwg, pos_min, pos_max, grid_xy):
        (grid_x1, grid_y1, grid_x2, grid_y2) = grid_xy

        # Getting range of vertical grid lines
        grid_pos_min = self._grid_interval * math.ceil(pos_min / self._grid_interval)
        grid_pos_max = self._grid_interval * math.floor(pos_max / self._grid_interval)
        if grid_pos_max % self._grid_interval == 0:
            grid_pos_max += 1

        # Adding vertical grid lines
        for grid_pos in range(grid_pos_min, grid_pos_max, self._grid_interval):
            grid_x = grid_x1 + (grid_x2 - grid_x1) * ((grid_pos - pos_min) / (pos_max - pos_min))
            dwg.add(svg.shapes.Line((grid_x, grid_y1), (grid_x, grid_y2),
                                stroke=self._grid_colour,
                                stroke_width=self._grid_size))

    def _add_grid_labels(self, dwg, pos_min, pos_max, map_xy):
        (map_x1, map_y1, map_x2, map_y2) = map_xy

        # Getting range of grid label positions
        grid_label_min = self._grid_label_interval * math.ceil(
            pos_min / self._grid_label_interval)
        grid_label_max = self._grid_label_interval * math.floor(
            pos_max / self._grid_label_interval)
        if grid_label_max % self._grid_label_interval == 0:
            grid_label_max += 1

        grid_label_gap = self._grid_label_rel_gap * self._im_h
        if self._hist_show:
            grid_label_gap = grid_label_gap + (self._hist_rel_gap + self._hist_rel_height) * self._im_h
                    
        if self._splithist_show:
            grid_label_gap = grid_label_gap + 3 * (self._hist_rel_gap + self._splithist_rel_height) * self._im_h 
            
        # Adding grid labels
        for grid_label_pos in range(grid_label_min, grid_label_max,
                                    self._grid_label_interval):
            grid_label_x = map_x1 + (map_x2 - map_x1) * (
                (grid_label_pos - pos_min) / (pos_max - pos_min))
            grid_label_y = map_y1 - self._dna_size / 2 - grid_label_gap
            rot = "rotate(%i,%i,%i)" % (-90, grid_label_x, grid_label_y)
            dwg.add(
                svg.text.Text(
                    str(grid_label_pos),
                    insert=(grid_label_x, grid_label_y),
                    transform=rot,
                    style=
                    f"text-anchor:start;font-family:{self._font};dominant-baseline:mathematical",
                    font_size=self._grid_label_size,
                    fill=self._grid_label_colour))

    def _add_dna(self, dwg, pos_min, pos_max, map_xy, ref=None):
        (map_x1, map_y1, map_x2, map_y2) = map_xy

        if self._dna_mode is DNA_MODE.SEQUENCE and ref is None:
            print(
                "WARNING: No reference path provided.  DNA rendered as solid lines."
            )
            self._dna_mode = DNA_MODE.LINE

        if self._dna_mode is DNA_MODE.SEQUENCE:
            # Draw DNA as text sequence
            dna_pos_interval = (map_x2 - map_x1) / (pos_max - pos_min)
            ref_c = ref.complement()
            for dna_pos in range(pos_min, pos_max + 1):
                dna_seq_x = map_x1 + dna_pos_interval * (dna_pos - pos_min)
                dwg.add(
                    svg.text.Text(
                        str(ref[dna_pos - 1]),
                        insert=(dna_seq_x, map_y1),
                        style=
                        f"text-anchor:middle;font-family:{self._font};dominant-baseline:mathematical",
                        font_size=self._dna_size,
                        fill=self._dna_colour))
                dwg.add(
                    svg.text.Text(
                        str(ref_c[dna_pos - 1]),
                        insert=(dna_seq_x, map_y2),
                        style=
                        f"text-anchor:middle;font-family:{self._font};dominant-baseline:mathematical",
                        font_size=self._dna_size,
                        fill=self._dna_colour))

        elif self._dna_mode is DNA_MODE.LINE:
            # Draw DNA as lines
            dwg.add(
                svg.shapes.Line((map_x1, map_y1), (map_x2, map_y1),
                                stroke=self._dna_colour,
                                stroke_width=self._dna_size,
                                style="stroke-linecap:square"))
            dwg.add(
                svg.shapes.Line((map_x1, map_y2), (map_x2, map_y2),
                                stroke=self._dna_colour,
                                stroke_width=self._dna_size,
                                style="stroke-linecap:square"))

    def _add_end_labels(self, dwg, map_xy):
        (map_x1, map_y1, map_x2, map_y2) = map_xy

        end_label_x1 = map_x1 - (self._end_label_rel_gap * self._im_w)
        end_label_x2 = map_x2 + (self._end_label_rel_gap * self._im_w)

        if self._end_label_position == VPOS.INSIDE:
            end_label_y1 = map_y1 - self._dna_size * 0.4  # + self._end_label_size * 0.75 - self._dna_size / 2
            end_label_y2 = map_y2 + self._dna_size * 0.4  # + self._dna_size / 2
            top_baseline = "hanging"
            bottom_baseline = ""
        elif self._end_label_position == VPOS.CENTRE:
            end_label_y1 = map_y1
            end_label_y2 = map_y2
            top_baseline = "mathematical"
            bottom_baseline = "mathematical"

        dwg.add(
            svg.text.Text(
                "5′",
                insert=(end_label_x1, end_label_y1),
                style=
                f"text-anchor:end;font-family:{self._font};dominant-baseline:{top_baseline}",
                font_size=self._end_label_size,
                fill=self._end_label_colour))
        dwg.add(
            svg.text.Text(
                "3′",
                insert=(end_label_x2, end_label_y1),
                style=
                f"text-anchor:start;font-family:{self._font};dominant-baseline:{top_baseline}",
                font_size=self._end_label_size,
                fill=self._end_label_colour))
        dwg.add(
            svg.text.Text(
                "3′",
                insert=(end_label_x1, end_label_y2),
                style=
                f"text-anchor:end;font-family:{self._font};dominant-baseline:{bottom_baseline}",
                font_size=self._end_label_size,
                fill=self._end_label_colour))
        dwg.add(
            svg.text.Text(
                "5′",
                insert=(end_label_x2, end_label_y2),
                style=
                f"text-anchor:start;font-family:{self._font};dominant-baseline:{bottom_baseline}",
                font_size=self._end_label_size,
                fill=self._end_label_colour))

    def _add_cbar(self, dwg, map_xy, freq):
        (map_x1, map_y1, map_x2, map_y2) = map_xy

        cbar_x1 = self._im_w * self._cbar_rel_left + self._cbar_border_size/2
        cbar_x2 = cbar_x1 + self._im_w * self._cbar_rel_width + self._cbar_border_size/2
        cbar_w = self._im_w * self._cbar_rel_width

        cbar_y1 = map_y1 + self._cbar_border_size/2
        cbar_y2 = map_y2 - self._cbar_border_size/2
            
        cbar_h = cbar_y2 - cbar_y1

        # Adding cbar lines
        cmap = cm.get_cmap(self._event_colourmap)

        line_offs = cbar_h / 255
        for i in range(0, 255):
            y1 = cbar_y2 - (i + 1) * line_offs if i == 254 else cbar_y2 - (
                i + 2) * line_offs
            y2 = cbar_y2 - (i) * line_offs
            rgba = cmap(i / 256)
            col = "rgb(%i,%i,%i)" % (rgba[0] * 255, rgba[1] * 255,
                                     rgba[2] * 255)
            dwg.add(
                svg.shapes.Rect(insert=(cbar_x1, y1),
                                size=(cbar_w, y2 - y1),
                                stroke="none",
                                fill=col))

        # Adding border around cbar
        dwg.add(
            svg.shapes.Rect(insert=(cbar_x1, cbar_y1),
                            size=(cbar_w, cbar_h),
                            stroke="black",
                            fill="none",
                            stroke_width=self._cbar_border_size,
                            style="stroke-linecap:square"))

        if self._cbar_label_show:
            if self._event_min_range == -1 and self._event_max_range == -1:
                total = sum(freq.values())
                event_min_range = math.floor(
                    (min(freq.values()) / total) * 100)
                event_max_range = math.ceil((max(freq.values()) / total) * 100)
            else:
                event_min_range = self._event_min_range
                event_max_range = self._event_max_range

            event_range = event_max_range - event_min_range

            # Preventing divide by zero errors
            if event_range == 0:
                event_range = 1

            end_label_x = cbar_x2 + (self._cbar_label_rel_gap * self._im_w)
            for label_value in range(event_min_range, event_max_range,
                                     self._cbar_label_interval):
                end_label_y = cbar_y2 - cbar_h * (
                    label_value - event_min_range) / event_range
                dwg.add(
                    svg.text.Text(
                        "%i%%" % label_value,
                        insert=(end_label_x, end_label_y),
                        style=
                        f"text-anchor:start;font-family:{self._font};dominant-baseline:mathematical",
                        font_size=self._cbar_label_size,
                        fill=self._cbar_label_colour))

            # Adding final label (keeping this separate means we still get it, even if it's not on the interval)
            end_label_y = cbar_y2 - cbar_h * (event_max_range -
                                              event_min_range) / event_range
            dwg.add(
                svg.text.Text(
                    "%i%%" % event_max_range,
                    insert=(end_label_x, end_label_y),
                    style=
                    f"text-anchor:start;font-family:{self._font};dominant-baseline:mathematical",
                    font_size=self._cbar_label_size,
                    fill=self._cbar_label_colour))

    def _add_sum_hist(self, dwg, pos_min, pos_max, map_xy, freq, n_events, is_top):
        (map_x1, map_y1, map_x2, map_y2) = map_xy
        
        hist_x1 = map_x1
        hist_x2 = map_x2
        
        sign = -1 if is_top else 1
        
        hist_y2 = (map_y1 if is_top else map_y2) + sign * (self._hist_rel_gap * self._im_h)
        if (self._splithist_show):
            hist_y2 = hist_y2 + sign * 3 * (self._hist_rel_gap + self._splithist_rel_height) * self._im_h

        hist_y1 = hist_y2 + sign * self._im_h * self._hist_rel_height

        hist_xy = (hist_x1, hist_y1, hist_x2, hist_y2)
        self._add_hist(dwg, pos_min, pos_max, hist_xy, freq, n_events, self._hist_min_range, self._hist_max_range, self._hist_label_interval, self._hist_grid_interval, "Sum", is_top)


    def _add_splithists(self, dwg, pos_min, pos_max, map_xy, freq_5, freq_3, freq_b, n_events, is_top):
        (map_x1, map_y1, map_x2, map_y2) = map_xy
                
        hist_x1 = map_x1
        hist_x2 = map_x2
        
        sign = -1 if is_top else 1

        min_range = self._splithist_min_range
        max_range = self._splithist_max_range
        label_interval = self._splithist_label_interval
        grid_interval = self._splithist_grid_interval
        
        # Split hist 1
        hist_y2 = (map_y1 if is_top else map_y2) + sign * (self._hist_rel_gap * self._im_h)
        hist_y1 = hist_y2 + sign * self._im_h * self._splithist_rel_height
        hist_xy = (hist_x1, hist_y1, hist_x2, hist_y2)
        self._add_hist(dwg, pos_min, pos_max, hist_xy, freq_5, n_events, min_range, max_range, label_interval, grid_interval, "5′ OH", is_top)

        # Split hist 2
        hist_y2 = (map_y1 if is_top else map_y2) + sign * (self._hist_rel_gap * self._im_h)
        hist_y2 = hist_y2 + sign * (self._hist_rel_gap + self._splithist_rel_height) * self._im_h
        hist_y1 = hist_y2 + sign * self._im_h * self._splithist_rel_height
        hist_xy = (hist_x1, hist_y1, hist_x2, hist_y2)
        self._add_hist(dwg, pos_min, pos_max, hist_xy, freq_3, n_events, min_range, max_range, label_interval, grid_interval, "3′ OH", is_top)

        # Split hist 3
        hist_y2 = (map_y1 if is_top else map_y2) + sign * (self._hist_rel_gap * self._im_h)
        hist_y2 = hist_y2 + sign * 2 * (self._hist_rel_gap + self._splithist_rel_height) * self._im_h
        hist_y1 = hist_y2 + sign * self._im_h * self._splithist_rel_height
        hist_xy = (hist_x1, hist_y1, hist_x2, hist_y2)
        self._add_hist(dwg, pos_min, pos_max, hist_xy, freq_b, n_events, min_range, max_range, label_interval, grid_interval, "Blunt", is_top)


    def _add_hist(self, dwg, pos_min, pos_max, hist_xy, freq, n_events, min_range, max_range, label_interval, grid_interval, name, is_top):
        if len(freq) == 0:
            return

        (hist_x1, hist_y1, hist_x2, hist_y2) = hist_xy
        bar_width = (hist_x2 - hist_x1) / (pos_max - pos_min)
        hist_h = abs(hist_y2-hist_y1)

        # Getting the total number of events
        if min_range == -1:
            hist_min_range = math.floor((min(freq.values()) / n_events) * 100)
        else:
            hist_min_range = min_range

        if max_range == -1:
            hist_max_range = math.ceil((max(freq.values()) / n_events) * 100)
        else:
            hist_max_range = max_range
            
        hist_range = hist_max_range - hist_min_range

        # Preventing divide by zero errors
        if hist_range == 0:
            hist_range = 1

        # Used for positioning labels and grid lines
        sign = -1 if is_top else 1
            
        # Adding optional grid lines
        if self._grid_show:
            self._add_grid_lines(dwg, pos_min, pos_max, hist_xy)

        # Adding the x-axis
        grid_x1 = hist_x1 - self._hist_overhang - self._grid_size/2
        grid_x2 = hist_x2 + self._hist_overhang - self._grid_size/2
        grid_y = hist_y2 + self._hist_grid_size / 2 - self._grid_size/2
        dwg.add(
            svg.shapes.Line((grid_x1, grid_y), (grid_x2, grid_y),
                            stroke=self._hist_grid_colour,
                            stroke_width=self._hist_grid_size))

        # # Adding the y-axis
        line_x = ((hist_x1 - self._hist_overhang) if self._hist_label_position == HPOS.LEFT else (hist_x2 + self._hist_overhang)) - self._grid_size/2
                
        dwg.add(
            svg.shapes.Line((line_x, hist_y1), (line_x, hist_y2),
                            stroke=self._hist_grid_colour,
                            stroke_width=self._hist_grid_size))

        # Adding the grid first, so it's at the bottom of the stack
        if self._hist_grid_show:
            for grid_value in range(grid_interval, hist_max_range, grid_interval):

                grid_y = hist_y2 + sign * hist_h * (grid_value - hist_min_range) / hist_range + self._hist_grid_size / 2 - self._grid_size/2
                dwg.add(svg.shapes.Line((grid_x1, grid_y), (grid_x2, grid_y),
                                    stroke=self._hist_grid_colour,
                                    stroke_width=self._hist_grid_size))

            grid_y = hist_y2 + sign * hist_h * (
                hist_max_range -
                hist_min_range) / hist_range + self._hist_grid_size / 2
            dwg.add(
                svg.shapes.Line((grid_x1, grid_y), (grid_x2, grid_y),
                                stroke=self._hist_grid_colour,
                                stroke_width=self._hist_grid_size))

        # Iterating over each position, adding the histogram bar
        for pos in range(pos_min-self._hist_bin_width, pos_max + 1):
            # Only do the plotting once per bin
            if pos % self._hist_bin_width != 0:
                continue
            
            bin_total = 0
            for curr_pos in range(pos, pos + self._hist_bin_width):
                bin_total = bin_total + (0 if freq.get(curr_pos) is None else
                                         freq.get(curr_pos))

            if bin_total == 0:
                continue

            event_pc = (bin_total / n_events) * 100
            norm_count = (event_pc - hist_min_range) / hist_range
            
            bar_offset = bar_width * self._hist_bin_width * self._hist_pc_bar_gap * 0.5 / 100
            bar_x1 = max(hist_x1, hist_x1 + bar_width * (pos - pos_min)) + bar_offset
            bar_x2 = min(hist_x2, hist_x1 + bar_width * (pos - pos_min + self._hist_bin_width)) - bar_offset
            
            if is_top:
                bar_y1 = hist_y2 - norm_count * hist_h
                bar_y2 = hist_y2
            else:
                bar_y1 = hist_y2
                bar_y2 = hist_y2 + norm_count * hist_h
                
            dwg.add(
                svg.shapes.Rect(insert=(bar_x1, bar_y1),
                                size=(bar_x2-bar_x1, bar_y2 - bar_y1),
                                stroke="none",
                                fill=self._hist_colour))
        
        if self._hist_label_show:
            if self._hist_label_position == HPOS.LEFT:
                anchor = "end"
                label_x = hist_x1 - (self._hist_label_rel_gap * self._im_w)
            else:
                anchor = "start"
                label_x = hist_x2 + (self._hist_label_rel_gap * self._im_w)

            hist_range_start = hist_min_range if self._hist_label_zero_show else self.label_interval
            for label_value in range(hist_range_start, hist_max_range, label_interval):
                end_label_y = hist_y2 + sign * hist_h * (label_value - hist_min_range) / hist_range
                dwg.add(svg.text.Text(
                        "%i%%" % label_value,
                        insert=(label_x, end_label_y),
                        style=
                        f"text-anchor:{anchor};font-family:{self._font};dominant-baseline:mathematical",
                        font_size=self._hist_label_size,
                        fill=self._hist_label_colour))

            # Adding final label (keeping this separate means we still get it, even if it's not on the interval)
            end_label_y = hist_y2 + sign * hist_h * (hist_max_range - hist_min_range) / hist_range
            dwg.add(svg.text.Text(
                    "%i%%" % hist_max_range,
                    insert=(label_x, end_label_y),
                    style=
                    f"text-anchor:{anchor};font-family:{self._font};dominant-baseline:mathematical",
                    font_size=self._hist_label_size,
                    fill=self._hist_label_colour))

        if self._hist_name_show:
            if self._hist_name_position == HPOS.LEFT:
                anchor = "end"
                name_x = hist_x1 - (self._hist_name_rel_gap * self._im_w)
            else:
                anchor = "start"
                name_x = hist_x2 + (self._hist_name_rel_gap * self._im_w)

            name_y = hist_y2 + sign * hist_h * 0.5
            dwg.add(svg.text.Text(name, insert=(name_x, name_y), style=f"text-anchor:{anchor};font-family:{self._font};dominant-baseline:mathematical",
                        font_size=self._hist_name_size, fill=self._hist_name_colour))

    def _add_splithist():
        print("TO ADD - SPLITHIST (could just re-use main _add_hist function")

    def _add_event_lines(self, dwg, pos_min, pos_max, map_xy, freq, ref):
        if len(freq) == 0:
            return

        total = sum(freq.values())
        if self._event_min_range == -1 and self._event_max_range == -1:
            event_min_range = math.floor((min(freq.values()) / total) * 100)
            event_max_range = math.ceil((max(freq.values()) / total) * 100)
        else:
            event_min_range = self._event_min_range
            event_max_range = self._event_max_range
        event_range = event_max_range - event_min_range

        # Preventing divide by zero errors
        if event_range == 0:
            event_range = 1

        # Getting reference length (or estimating)
        if ref is None:
            pos_ranges = aw.get_event_pos_range(freq, round=1)
            ref_len = max(pos_ranges[1], pos_ranges[3])
        else:
            ref_len = len(ref)

        for (cleavage_site_t, cleavage_site_b, split) in freq.keys():
            # Checking this event is within the rendered range
            if self._event_outside_range_show:
                if split:
                    # Skipping any split events which don't span the displayed range
                    if (cleavage_site_t < pos_min or cleavage_site_t >= pos_max) and (cleavage_site_b < pos_min or cleavage_site_b >= pos_max):
                        continue
                else:
                    # Skipping any non-split events which don't span the displayed range
                    if (cleavage_site_t < pos_min and cleavage_site_b < pos_min) or (cleavage_site_t >= pos_max and cleavage_site_b >= pos_max):
                        continue
            else:
                # Skipping any events which don't start or end in the displayed range
                if (cleavage_site_t < pos_min or cleavage_site_t >= pos_max) or (cleavage_site_b < pos_min or cleavage_site_b >= pos_max):
                    continue

            event_pc = (freq.get(
                (cleavage_site_t, cleavage_site_b, split)) / total) * 100
            norm_count = (event_pc - event_min_range) / event_range

            if split:
                self._add_split_line(dwg, cleavage_site_t, cleavage_site_b,
                                     pos_min, pos_max, map_xy, norm_count,
                                     ref_len)
            else:
                self._add_continuous_line(dwg, cleavage_site_t,
                                          cleavage_site_b, pos_min, pos_max,
                                          map_xy, norm_count)

    def _add_continuous_line(self, dwg, cleavage_site_t, cleavage_site_b,
                             pos_min, pos_max, map_xy, norm_count):
        (map_x1, map_y1, map_x2, map_y2) = map_xy
        cmap = cm.get_cmap(self._event_colourmap)

        # Adding line (adding width 1 to ensure everything is visible)
        event_width = (self._event_max_size - self._event_min_size
                       ) * norm_count + self._event_min_size

        event_t_x = map_x1 + (map_x2 - map_x1) * (
            (cleavage_site_t + 0.5 - pos_min) / (pos_max - pos_min))
        event_t_y = map_y1 + self._dna_size / 2 + self._dna_rel_gap*self._im_h

        event_b_x = map_x1 + (map_x2 - map_x1) * (
            (cleavage_site_b + 0.5 - pos_min) / (pos_max - pos_min))
        event_b_y = map_y2 - self._dna_size / 2 - self._dna_rel_gap*self._im_h

        (event_t_x, event_t_y, event_b_x, event_b_y) = self._crop_events_to_range(event_t_x, event_t_y, event_b_x, event_b_y, map_xy)

        rgba = cmap(norm_count)
        col = "rgb(%i,%i,%i)" % (rgba[0] * 255, rgba[1] * 255, rgba[2] * 255)

        dwg.add(
            svg.shapes.Line((event_t_x, event_t_y), (event_b_x, event_b_y),
                            stroke=col,
                            stroke_width=event_width,
                            style="stroke-linecap:square;stroke-opacity:%f" %
                            self._event_opacity))

    def _add_split_line(self, dwg, cleavage_site_t, cleavage_site_b, pos_min,
                        pos_max, map_xy, norm_count, ref_len):
        (map_x1, map_y1, map_x2, map_y2) = map_xy
        cmap = cm.get_cmap(self._event_colourmap)

        # Adding line (adding width 1 to ensure everything is visible)
        event_width = (self._event_max_size - self._event_min_size
                       ) * norm_count + self._event_min_size

        # Processing is dependent on which site is lower
        if cleavage_site_t < cleavage_site_b:
            site_sep = ref_len - (cleavage_site_b - cleavage_site_t)

            event_t_x1 = map_x1
            event_t_y1 = map_y1 + (map_y2 - map_y1 - self._dna_rel_gap*self._im_h*2) * (
                (cleavage_site_t + 0.5 - pos_min) / site_sep)
            event_t_x2 = map_x1 + (map_x2 - map_x1) * (
                (cleavage_site_t + 0.5 - pos_min) / (pos_max - pos_min))
            event_t_y2 = map_y1 + self._dna_size / 2 + self._dna_rel_gap*self._im_h

            event_b_x1 = map_x1 + (map_x2 - map_x1) * (
                (cleavage_site_b + 0.5 - pos_min) / (pos_max - pos_min))
            event_b_y1 = map_y2 - self._dna_size / 2 - self._dna_rel_gap*self._im_h
            event_b_x2 = map_x2
            event_b_y2 = map_y2 - (map_y2 - map_y1 - self._dna_rel_gap*self._im_h*2) * (
                (pos_max - (cleavage_site_b + 0.5)) / site_sep)

        else:
            site_sep = ref_len - (cleavage_site_t - cleavage_site_b)

            event_t_x1 = map_x1 + (map_x2 - map_x1) * (
                (cleavage_site_t + 0.5 - pos_min) / (pos_max - pos_min))
            event_t_y1 = map_y1 + self._dna_size / 2 + self._dna_rel_gap*self._im_h
            event_t_x2 = map_x2
            event_t_y2 = map_y1 + (map_y2 - map_y1 - self._dna_rel_gap*self._im_h*2) * (
                (pos_max - (cleavage_site_t + 0.5)) / site_sep)

            event_b_x1 = map_x1
            event_b_y1 = map_y2 - (map_y2 - map_y1 - self._dna_rel_gap*self._im_h*2) * (
                (cleavage_site_b + 0.5 - pos_min) / site_sep)
            event_b_x2 = map_x1 + (map_x2 - map_x1) * (
                (cleavage_site_b + 0.5 - pos_min) / (pos_max - pos_min))
            event_b_y2 = map_y2 - self._dna_size / 2 - self._dna_rel_gap*self._im_h

        rgba = cmap(norm_count)
        col = "rgb(%i,%i,%i)" % (rgba[0] * 255, rgba[1] * 255, rgba[2] * 255)

        dwg.add(
            svg.shapes.Line((event_t_x1, event_t_y1), (event_t_x2, event_t_y2),
                            stroke=col,
                            stroke_width=event_width,
                            style="stroke-linecap:square;stroke-opacity:%f" %
                            self._event_opacity))
        dwg.add(
            svg.shapes.Line((event_b_x1, event_b_y1), (event_b_x2, event_b_y2),
                            stroke=col,
                            stroke_width=event_width,
                            style="stroke-linecap:square;stroke-opacity:%f" %
                            self._event_opacity))

    def _crop_events_to_range(self, t_x_in, t_y_in, b_x_in, b_y_in, map_xy):
        (map_x1, map_y1, map_x2, map_y2) = map_xy

        t_x_out = t_x_in
        t_y_out = t_y_in
        b_x_out = b_x_in
        b_y_out = b_y_in

        if t_x_in < map_x1:
            t_x_out = map_x1
            t_y_out = ((map_x1-t_x_in)/abs(b_x_in-t_x_in))*(b_y_in-t_y_in) + t_y_in
            
        if t_x_in > map_x2:
            t_x_out = map_x2
            t_y_out = ((t_x_in-map_x2)/abs(b_x_in-t_x_in))*(b_y_in-t_y_in) + t_y_in

        if b_x_in < map_x1:
            b_x_out = map_x1
            b_y_out = ((t_x_in-map_x1)/abs(b_x_in-t_x_in))*(b_y_in-t_y_in) + t_y_in
            
        if b_x_in > map_x2:
            b_x_out = map_x2
            b_y_out = ((map_x2-t_x_in)/abs(b_x_in-t_x_in))*(b_y_in-t_y_in) + t_y_in

        return (t_x_out, t_y_out, b_x_out, b_y_out)
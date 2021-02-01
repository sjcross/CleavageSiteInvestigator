from matplotlib import cm
from utils import abstractmapwriter as amw

import datetime as dt
import math
import os


class HeatMapWriterSVG(amw.AbstractMapWriter):
    def __init__(self, im_dim=800, sum_show=True):
        self._im_dim = im_dim
        
        self._sum_show = sum_show

    def write_map(self, out_path, freq, ref, pos_range, append_dt):
        # Getting pos ranges to plot based on available information
        pos_range = amw.get_pos_range(freq, ref, pos_range)
        if pos_range is None:
            return

        # Creating output CSV document
        root_name = os.path.splitext(out_path)[0]
        datetime_str = dt.datetime.now().strftime("_%Y-%m-%d_%H-%M-%S") if append_dt else ""
        outname = root_name+ datetime_str + '.' + 'csv'
                
        # Defining limits of heat map
        pos_t_range = pos_range[1]-pos_range[0]
        pos_b_range = pos_range[3]-pos_range[2]
               

    
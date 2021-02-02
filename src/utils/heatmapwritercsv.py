from utils import abstractmapwriter as amw
from utils import fileutils as fu

import os


class HeatMapWriterCSV(amw.AbstractMapWriter):
    ## CONSTRUCTOR
    
    def __init__(self, im_dim=800, event_label_decimal_places = 4, sum_show=True, count_show=True):
        self._im_dim = im_dim
        
        self._number_format = "%%.%if" % event_label_decimal_places
        self._sum_show = sum_show
        self._count_show = count_show

    def write_map(self, out_path, freq, ref, pos_range, append_dt):
        # Getting pos ranges to plot based on available information
        pos_range = amw.get_pos_range(freq, ref, pos_range)
        if pos_range is None:
            return

        # Creating output CSV document
        root_name = os.path.splitext(out_path)[0]
        file = fu.open_file(root_name, '', 'csv', append_dt)
                
        # Determining total events in the given position range
        sum_events = amw.get_sum_events(pos_range, freq)
        
        # If showing sum, calculating relevant statistics
        (freq_t, freq_b) = amw.get_full_sequence_summed_frequency(freq, pos_range)

        # Initialising string
        str = self._get_header_row(pos_range[0],pos_range[1])
        file.write(str)

        # Iterating over each bottom-strand position and adding as a new row
        for pos_b in range(pos_range[2], pos_range[3]+1):
            str = self._get_row(freq, pos_b, pos_range[0], pos_range[1], sum_events, freq_b)
            file.write(str)

        # If necessary, showing the sum row
        if self._sum_show:
            str = self._get_sum_row(pos_range[0], pos_range[1], sum_events, freq_t)
            file.write(str)

        # If necessary, showing the totala event count
        if self._count_show:
            # Adding a blank line
            file.write("\n")
            file.write(f"N = {sum_events}")

        file.close()
               
    def _get_header_row(self, pos_t_min, pos_t_max):
        # The first element shows the axes
        row = "B V : T >"

        # Iterating over each top strand position
        for pos_t in range(pos_t_min, pos_t_max+1):
            row = row + "," + str(pos_t)

        # If showing sum row and column, add sum heading
        if self._sum_show:
            row = row + ",Sum"

        # Returning with an end line character
        return row + "\n"
    
    def _get_row(self, freq, pos_b, pos_t_min, pos_t_max, sum_events, freq_b=None):
        # The first element of the row is the bottom strand index
        row = str(pos_b)

        # Iterating over each top-strand position, adding as a new column
        for pos_t in range(pos_t_min, pos_t_max+1):
            event_pc = amw.get_event_pc((pos_t,pos_b), freq, sum_events)
            row = row + "," + str(self._number_format % event_pc)

        # Adding sum if necessary
        if self._sum_show:
            event_pc = amw.get_event_pc(pos_b, freq_b, sum_events)
            row = row + "," + str(self._number_format % event_pc)

        return row + "\n"

    def _get_sum_row(self, pos_t_min, pos_t_max, sum_events, freq_t):
        # The first element of the row is the word "Sum"
        row = "Sum"

        # Adding the sum for each columns
        for pos_t in range(pos_t_min, pos_t_max+1):
            event_pc = amw.get_event_pc(pos_t, freq_t, sum_events)
            row = row + "," + str(self._number_format % event_pc)

        return row + "\n"
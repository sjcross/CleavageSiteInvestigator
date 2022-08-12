from matplotlib.pyplot import hist
from enums.eventtype import Eventtype

from utils import abstractwriter as aw
from utils import fileutils as fu
from utils import reportutils as ru

import os


class HistogramWriterCSV(aw.AbstractWriter):
    ## CONSTRUCTOR

    def __init__(self, hist_bin_width=1):
        self._hist_bin_width = hist_bin_width
    
    def write(self, out_path, freq, ref, pos_range, append_dt, commandStr):
        # Getting pos range to plot based on available information
        (pos_min, pos_max) = aw.get_single_pos_range(freq, ref, pos_range)
        
        # Creating output CSV document
        root_name = os.path.splitext(out_path)[0]
        file = fu.open_file(root_name, '', 'csv', append_dt)

        # Getting the frequencies of each event type
        (freq_t, freq_b) = ru.get_position_frequency(freq)  
        (freq_t_5, freq_b_5) = ru.get_position_frequency(freq, Eventtype.FIVE_P)
        (freq_t_3, freq_b_3) = ru.get_position_frequency(freq, Eventtype.THREE_P)
        (freq_t_b, freq_b_b) = ru.get_position_frequency(freq, Eventtype.BLUNT)  
                
        # Initialising string
        str = "Position, Freq. all (top), Freq. 5' (top), Freq. Blunt (top), Freq. 3' (top), Freq. all (bottom), Freq. 5' (bottom), Freq. Blunt (bottom), Freq. 3' (bottom)\n"
        file.write(str)

        # Iterating over each position, adding the histogram bar
        for pos in range(pos_min-self._hist_bin_width+1, pos_max + 1):
            # Only do the plotting once per bin
            if pos % self._hist_bin_width != 0:
                continue
            
            str = "%i," % pos
            str = str + "%i," % self._get_bin(freq_t, pos)
            str = str + "%i," % self._get_bin(freq_t_5, pos)
            str = str + "%i," % self._get_bin(freq_t_b, pos)
            str = str + "%i," % self._get_bin(freq_t_3, pos)
            str = str + "%i," % self._get_bin(freq_b, pos)
            str = str + "%i," % self._get_bin(freq_b_5, pos)
            str = str + "%i," % self._get_bin(freq_b_b, pos)
            str = str + "%i," % self._get_bin(freq_b_3, pos)
            str = str + "\n"

            file.write(str)

        # Adding argument line
        if (commandStr != ''):
            file.write('\n')
            file.write('Command: '+commandStr+'\n')

        file.close()
        
    def _get_bin(self, freq, pos):
        bin_total = 0
        for curr_pos in range(pos, pos + self._hist_bin_width):
            bin_total = bin_total + (0 if freq.get(curr_pos) is None else freq.get(curr_pos))

        return bin_total

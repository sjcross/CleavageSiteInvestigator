import math
import sys

from abc import abstractmethod
from utils import csvutils as cu  
from utils import fileutils as fu
from utils import reportutils as ru


class AbstractMapWriter():
    ## PUBLIC METHODS

    def write_map_from_file(self, csv_path, out_path, ref_path="", pos_range=None, append_dt=False):
        # Loading results from CSV
        cr = cu.CSVReader()

        file_type = cu.get_file_type(csv_path)
        if file_type == cu.FileTypes.INDIVIDUAL:
            results = cr.read_individual(csv_path)
            freq = ru.get_full_sequence_frequency(results)
        elif file_type == cu.FileTypes.SUMMARY:
            freq = cr.read_summary(csv_path)

        # Loading reference sequence if path provided
        if ref_path == "":
            ref = None
        else:
            filereader = fu.FileReader(verbose=False)
            ref = filereader.read_sequence(ref_path)[0][0][0]
            
            self.write_map(out_path, freq, ref=ref, pos_range=pos_range, append_dt=append_dt)

    @abstractmethod
    def write_map(self, out_path, freq, ref=None, pos_range=None, append_dt=False):
        pass

def get_single_pos_range(freq, ref, pos_range):
    if pos_range is None:
        if ref is None:
            # Rounding up to the nearest 10
            pos_ranges = get_event_pos_range(freq, round=10)
            pos_min = min(pos_ranges[0],pos_ranges[2])
            pos_max = max(pos_ranges[1],pos_ranges[3])
        else:
            pos_min = 0
            pos_max = len(ref)-1
    else:
        pos_min = pos_range[0]
        pos_max = pos_range[1]

    # Checking pos_min and pos_max are different (to prevent divide by zero errors)
    if pos_min == pos_max:
        print("WARNING: Min and max sequence positions must be different")
        return

    return (pos_min, pos_max)

def get_double_pos_range(freq, ref, pos_range):
    if pos_range is None:        
        if ref is None:
            # Rounding up to the nearest 10
            (pos_t_min,pos_t_max,pos_b_min,pos_b_max) = get_event_pos_range(freq, round=10)
        else:
            pos_t_min = 0
            pos_b_min = 0
            pos_t_max = len(ref)-1
            pos_b_max = len(ref)-1
    else:
        pos_t_min = pos_range[0]
        pos_t_max = pos_range[1]
        pos_b_min = pos_range[2]
        pos_b_max = pos_range[3]
        
    # Checking pos_min and pos_max are different (to prevent divide by zero errors)
    if pos_t_min == pos_t_max or pos_b_min == pos_b_max:
        print("WARNING: Min and max sequence positions must be different")
        return

    return (pos_t_min, pos_t_max, pos_b_min, pos_b_max)

def get_event_pos_range(freq, round=1):
    pos_t_min = sys.maxsize
    pos_t_max = 0
    pos_b_min = sys.maxsize
    pos_b_max = 0

    for (cleavage_site_t, cleavage_site_b,split) in freq.keys():
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

def get_max_events(pos_range, freq, sum_show):
    (pos_t_min, pos_t_max, pos_b_min, pos_b_max) = pos_range

    max_events = 0
    if sum_show:
        (freq_t, freq_b) = get_full_sequence_summed_frequency(freq, pos_range)

        for t in freq_t.keys():
            if t >= pos_t_min and t <= pos_t_max:
                max_events = max(max_events,freq_t.get(t))

        for b in freq_b.keys():
            if b >= pos_b_min and b <= pos_b_max:
                max_events = max(max_events,freq_b.get(b))

    else:
        for (t,b,split) in freq.keys():
            if t >= pos_t_min and t <= pos_t_max and b >= pos_b_min and b <= pos_b_max:
                max_events = max(max_events,freq.get((t,b,split)))

    return max_events

def get_sum_events(pos_range, freq):
    (pos_t_min, pos_t_max, pos_b_min, pos_b_max) = pos_range

    sum_events = 0
    for (t,b,split) in freq.keys():
            if t >= pos_t_min and t <= pos_t_max and b >= pos_b_min and b <= pos_b_max:
                sum_events = sum_events + freq.get((t,b,split))

    return sum_events

def get_full_sequence_summed_frequency(freq_full, pos_range):
    (pos_t_min, pos_t_max, pos_b_min, pos_b_max) = pos_range

    freq_t = {}
    freq_b = {}

    for (t, b, split) in freq_full.keys():        
        freq = freq_full.get((t, b, split))
        
        # Adding this frequency to the relevant elements of freq_t and freq_b
        if t >= pos_t_min and t <= pos_t_max and b >= pos_b_min and b <= pos_b_max:
            freq_t[t] = freq_t[t] + freq if t in freq_t else freq
            freq_b[b] = freq_b[b] + freq if b in freq_b else freq

    return (freq_t, freq_b)

def get_event_norm_count(key, freq, max_events):
    if key in freq:
        norm_count = freq.get(key)/max_events
    else:
        norm_count = 0

    return norm_count

def get_event_pc(key, freq, sum_events):
    if key in freq:
        event_pc = 100*freq.get(key)/sum_events  
    else:
        event_pc = 0

    return event_pc

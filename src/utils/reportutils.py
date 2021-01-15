import datetime as dt
import io
import math
import numpy as np
import sys

from enum import Enum
from utils import fileutils as fu
from utils import sequenceutils as su

class StrandMode(Enum):
    TOP = 1
    BOTTOM = 2
    BOTH = 3

class LocalMode(Enum):
    BOTH = 1
    FIVE_P = 2
    THREE_P = 3

class StdOut(object):
    def __init__(self,root_name,append_dt=False):
        self._out_file = fu.open_file(root_name, '_output', 'txt', append_dt=append_dt)
        self._sys_out = sys.stdout

    def shutdown(self):
       self._out_file.close()
       sys.stdout = self._sys_out

    def flush(self):
        self._sys_out.flush()

    def write(self,string):
        self._out_file.write(string)
        self._sys_out.write(string)

def get_full_sequence_frequency(results):
    freq = {}

    for (cleavage_site_t, cleavage_site_b, local_site_t, local_site_b) in results.values():
        k = (cleavage_site_t, cleavage_site_b)
        if k not in freq:
            freq[(cleavage_site_t, cleavage_site_b)] = 1
        else:
            freq[(cleavage_site_t, cleavage_site_b)] = freq[(cleavage_site_t, cleavage_site_b)] + 1

    # Sorting results by frequency
    freq = sort_results(freq)

    return freq

def get_local_sequence_frequency(results, strand_mode, local_mode, local_r):
    n_nt = 2*local_r

    # Initialise frequency store
    if local_mode is not LocalMode.BOTH:
        n_nt = math.floor(n_nt / 2)
        
    freq = _init_frequency_dict(n_nt)
    
    for (cleavage_site_t, cleavage_site_b, local_site_t, local_site_b) in results.values():
        if local_mode is LocalMode.FIVE_P:
            local_site_t = local_site_t[0:n_nt]
            local_site_b = local_site_b[0:n_nt]
        elif local_mode is LocalMode.THREE_P:
            local_site_t = local_site_t[-n_nt:]
            local_site_b = local_site_b[-n_nt:]

        if strand_mode is StrandMode.TOP:
            freq[local_site_t] = freq[local_site_t] + 1
        elif strand_mode is StrandMode.BOTTOM:
            freq[local_site_b] = freq[local_site_b] + 1
        elif strand_mode is StrandMode.BOTH:
            freq[local_site_t] = freq[local_site_t] + 1
            freq[local_site_b] = freq[local_site_b] + 1

    return freq

def get_sequence_cooccurrence(results, local_r):
    n_nt = 2*local_r

    # Create list of sequences (this is in the order we will use)
    labels = list(_init_frequency_dict(n_nt).keys())
        
    # Creating a 2D Numpy array
    freq = np.zeros((pow(4, n_nt), pow(4, n_nt)))

    # Iterating over each result, adding one to the matrix
    for (cleavage_site_t, cleavage_site_b, local_site_t, local_site_b) in results.values():
        idx_top = labels.index(local_site_t)
        idx_bottom = labels.index(local_site_b)

        freq[idx_top, idx_bottom] = freq[idx_top, idx_bottom] + 1
        
    return (labels,freq)    

def print_full_sequence_frequency(ref, freq, extra_nt = 0, offset=""):
    total = sum(list(freq.values()))
    
    # Displaying sequence frequency results
    for cleavage_sites in freq.keys():
        cleavage_site_t = cleavage_sites[0]
        cleavage_site_b = cleavage_sites[1]
        count = freq.get(cleavage_sites)

        print_position(cleavage_site_t, cleavage_site_b, offset=offset)
        print_count(count, total, offset=offset)
        print_type(cleavage_site_t, cleavage_site_b, offset=offset)
        print_sequence(ref, cleavage_site_t, cleavage_site_b, extra_nt=extra_nt, offset=offset)
        
    print("\n")

def print_local_sequence_frequency(freq, nonzero_only=False, offset=""):   
    for local_seq in freq.keys():
        if not nonzero_only or freq.get(local_seq) > 0:
            print("    %s: %i" % (local_seq, freq.get(local_seq)))

    print("\n")

def print_position(cleavage_site_t, cleavage_site_b, offset=""):
    if cleavage_site_b is None or cleavage_site_t is None:
        return

    print("%s    Position:    %i, %i" % (offset, cleavage_site_t, cleavage_site_b))

def print_count(count, total, offset=""):
    print("%s    Count:       %i/%i (%.1f%% of events)" % (offset, count, total, 100*(count/total)))

def print_type(cleavage_site_t, cleavage_site_b, offset=""):
    type_str = su.get_type_str(cleavage_site_t, cleavage_site_b)
    if type_str is None:
        return

    print("%s    Type:        %s" % (offset, type_str))

def print_sequence(ref, cleavage_site_t, cleavage_site_b, extra_nt=0, offset=""):
    if cleavage_site_b is None or cleavage_site_t is None:
        return

    (seq1, seq2) = su.get_sequence_str(ref, cleavage_site_t, cleavage_site_b, extra_nt=extra_nt)

    print("%s    Sequence:    %s\r\n                 %s%s\r\n" % (offset, seq1, offset, seq2))
        
def print_error_rate(error_count, sample_count, offset=""):
    error_rate = 100*error_count/sample_count
    print("%sCompleted with %i errors (%f%%)" % (offset, error_count, error_rate))
    
def sort_results(results):
    sorted_results = sorted(results.items(), key=lambda x: x[1], reverse=True)
    
    return dict(sorted_results)

def _init_frequency_dict(n_nt):
    freq = {}

    for i in range(pow(4,n_nt)):
        st = ""
        for j in range(n_nt-1,-1,-1):
            val = math.floor((i/pow(4,j))%4)
            if val == 0:
                st = st + "A"
            elif val == 1:
                st = st + "T"
            elif val == 2:
                st = st + "G"
            elif val == 3:
                st = st + "C"
        freq[st] = 0
    
    return freq

import math
import numpy as np

from enum import Enum

class StrandMode(Enum):
    TOP = 1
    BOTTOM = 2
    BOTH = 3

class LocalMode(Enum):
    BOTH = 1
    FIVE_P = 2
    THREE_P = 3

def get_full_sequence_frequency(results):
    freq = {}

    for (clevage_site_t, clevage_site_b, local_site_t, local_site_b) in results.values():
        k = (clevage_site_t, clevage_site_b)
        if k not in freq:
            freq[(clevage_site_t, clevage_site_b)] = 1
        else:
            freq[(clevage_site_t, clevage_site_b)] = freq[(clevage_site_t, clevage_site_b)] + 1

    # Sorting results by frequency
    freq = sort_results(freq)

    return freq

def get_local_sequence_frequency(results, strand_mode, local_mode):
    # Initialise frequency store based on number of nucleotides in first sequence
    if len(results) == 0:
        # If the results table doesn't include any data, assume a 2nt local sequence
        n_nt = 2
    else:
        n_nt = len(results[0][2])

    if local_mode is not LocalMode.BOTH:
        n_nt = math.floor(n_nt / 2)
        
    freq = _init_frequency_dict(n_nt)

    for (clevage_site_t, clevage_site_b, local_site_t, local_site_b) in results.values():
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

def get_sequence_cooccurrence(results):
    # Initialise frequency store based on number of nucleotides in first sequence
    if len(results) == 0:
        # If the results table doesn't include any data, assume a 2nt local sequence
        n_nt = 2
    else:
        n_nt = len(results[0][2])

    # Create list of sequences (this is in the order we will use)
    labels = list(_init_frequency_dict(n_nt).keys())
        
    # Creating a 2D Numpy array
    freq = np.zeros((pow(4, n_nt), pow(4, n_nt)))

    # Iterating over each result, adding one to the matrix
    for (clevage_site_t, clevage_site_b, local_site_t, local_site_b) in results.values():
        idx_top = labels.index(local_site_t)
        idx_bottom = labels.index(local_site_b)

        freq[idx_top, idx_bottom] = freq[idx_top, idx_bottom] + 1
        
    return (labels,freq)    

def print_full_sequence_frequency(ref, freq, offset=""):
    # Displaying sequence frequency results
    for clevage_sites in freq.keys():
        clevage_site_t = clevage_sites[0]
        clevage_site_b = clevage_sites[1]
        count = freq.get(clevage_sites)

        print_position(clevage_site_t, clevage_site_b,offset=offset)
        print_count(count,offset=offset)
        print_type(clevage_site_t, clevage_site_b,offset=offset)
        print_sequence(ref, clevage_site_t, clevage_site_b, offset=offset)
        
    print("\n")

def print_local_sequence_frequency(freq, nonzero_only=False, offset=""):   
    for local_seq in freq.keys():
        if not nonzero_only or freq.get(local_seq) > 0:
            print("    %s: %i" % (local_seq, freq.get(local_seq)))

    print("\n")

def print_position(clevage_site_t, clevage_site_b, offset=""):
    if clevage_site_b is None or clevage_site_t is None:
        return

    print("%s    Position:    %i, %i" % (offset, clevage_site_b, clevage_site_t))

def print_count(count, offset=""):
    print("%s    Count:       %i" % (offset, count))

def print_type(clevage_site_t, clevage_site_b, offset=""):
    if clevage_site_b is None or clevage_site_t is None:
        return

    # Identifying break type
    if clevage_site_b < clevage_site_t:
        print("%s    Type:        3' overhang" % offset)
        
    elif (clevage_site_b == clevage_site_t):
        print("%s    Type:        Blunt end" % offset)

    elif clevage_site_b > clevage_site_t:
        print("%s    Type:        5' overhang" % offset)

def print_sequence(ref, clevage_site_t, clevage_site_b, offset=""):
    if clevage_site_b is None or clevage_site_t is None:
        return

    # Identifying break type
    if clevage_site_b < clevage_site_t:
        # 3' overhang
        left_seq1 = ref[clevage_site_b - 1 :clevage_site_b]
        mid_seq1 = ref[clevage_site_b:clevage_site_t]
        right_seq1 = ref[clevage_site_t: clevage_site_t + 1]

        left_seq2 = left_seq1.complement()
        mid_seq2 = mid_seq1.complement()
        right_seq2 = right_seq1.complement()

        print("%s    Sequence:    5'...%s %s↓%s...3'\r\n                 %s3'...%s↑%s %s...5'\r\n" % (offset, left_seq1, mid_seq1, right_seq1, offset, left_seq2, mid_seq2, right_seq2))
        
    elif clevage_site_b == clevage_site_t:
        # Blunt end
        left_seq1 = ref[clevage_site_b - 3 :clevage_site_b]
        right_seq1 = ref[clevage_site_t: clevage_site_t + 3]

        left_seq2 = left_seq1.complement()
        right_seq2 = right_seq1.complement()

        print("%s    Sequence:    5'...%s↓%s...3'\r\n                 %s3'...%s↑%s...5'\r\n" % (offset, left_seq1, right_seq1, offset, left_seq2, right_seq2))

    elif clevage_site_b > clevage_site_t:
        # 5' overhang
        left_seq1 = ref[clevage_site_t - 1 :clevage_site_t]
        mid_seq1 = ref[clevage_site_t:clevage_site_b]
        right_seq1 = ref[clevage_site_b: clevage_site_b + 1]

        left_seq2 = left_seq1.complement()
        mid_seq2 = mid_seq1.complement()
        right_seq2 = right_seq1.complement()

        print("%s    Sequence:    5'...%s↓%s %s...3'\r\n                 %s3'...%s %s↑%s...5'\r\n" % (offset, left_seq1, mid_seq1, right_seq1, offset, left_seq2, mid_seq2, right_seq2))

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

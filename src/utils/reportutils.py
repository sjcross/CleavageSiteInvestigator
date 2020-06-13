def sort_results(results):
    sorted_results = sorted(results.items(), key=lambda x: x[1], reverse=True)
    
    return dict(sorted_results)

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
    
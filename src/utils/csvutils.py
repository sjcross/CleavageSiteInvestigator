from utils import sequenceutils as su

def write_csv(filename, results, ref, extra_nt):
    file = open(filename+".csv", 'w')

    # Initialising string
    str = create_header_line()
    file.write(str)

    # Iterating over each result, adding it as a new line
    for index, result in enumerate(results.values()):
        file.write(create_result_line(result, index, ref, extra_nt))

    # Writing csv to file
    file.close()

def create_header_line():
    return "INDEX, TYPE, TOP_POS, BOTTOM_POS, TOP_LOCAL_SEQ, BOTTOM_LOCAL_SEQ, TOP_SEQ, BOTTOM_SEQ"

def create_result_line(result, index, ref, extra_nt):
    # Initialising string for this line with the result index
    new_str = str(index) + ','

    if result is None:
        return new_str + '\n'

    print(result)
    (cleavage_site_t, cleavage_site_b, local_site_t, local_site_b) = result

    # Adding cleavage type
    new_str = new_str + su.get_type_str(cleavage_site_t, cleavage_site_b) + ','

    # Adding top and bottom cleavage positions
    new_str = new_str + str(cleavage_site_t) + ','
    new_str = new_str + str(cleavage_site_b) + ','

    # Adding top and bottom cleavage local sequences
    new_str = new_str + local_site_t + ','
    new_str = new_str + local_site_b + ','

    # Adding sequences
    (seq1, seq2) = su.get_sequence_str(ref, cleavage_site_t, cleavage_site_b, extra_nt=extra_nt)
    new_str = new_str + seq1 + ','
    new_str = new_str + seq2

    return new_str + '\n'
    
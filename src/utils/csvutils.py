import datetime as dt
import io

from utils import sequenceutils as su

def write_individual(root_name, results, ref, extra_nt, double_line_mode=False):
    # Checking if the file is available for writing
    outname = root_name+"_individual.csv"
    
    datetime_str = dt.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    outname = root_name+"_individual_" + datetime_str + ".csv"
    file = io.open(outname, "w", encoding="utf-8")

    # Initialising string
    str = get_header_line(double_line_mode)
    file.write(str)

    # Iterating over each result, adding it as a new line
    for index, result in enumerate(results.values()):
        str = get_result_line(result, index, ref, extra_nt, double_line_mode)
        file.write(str)

    file.close()

def get_header_line(double_line_mode):
    str = "INDEX, TYPE, TOP_POS, BOTTOM_POS, TOP_LOCAL_SEQ, BOTTOM_LOCAL_SEQ, "

    if double_line_mode:
        return str + " SEQUENCE\n"
    else:
        return str + "TOP_SEQUENCE, BOTTOM_SEQUENCE\n"

def get_result_line(result, index, ref, extra_nt, double_line_mode):
    # Initialising string for this line with the result index
    new_str = str(index+1) + ','

    if result is None:
        return new_str + '\n'

    (cleavage_site_t, cleavage_site_b, local_site_t, local_site_b) = result

    # Adding cleavage type
    new_str = new_str + su.get_type_str(cleavage_site_t, cleavage_site_b) + ','

    # Adding top and bottom cleavage positions
    new_str = new_str + str(cleavage_site_t) + ','
    new_str = new_str + str(cleavage_site_b) + ','

    # Adding top and bottom cleavage local sequences
    new_str = new_str + str(local_site_t) + ','
    new_str = new_str + str(local_site_b) + ','

    # Adding sequences
    (seq1, seq2) = su.get_sequence_str(ref, cleavage_site_t, cleavage_site_b, extra_nt=extra_nt)
    if double_line_mode:
        new_str = new_str + seq1 + '\n,,,,,,' + seq2
    else:
        new_str = new_str + seq1 + ','
        new_str = new_str + seq2

    return new_str + '\n'
    
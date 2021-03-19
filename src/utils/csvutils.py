from enum import Enum
from utils import fileutils as fu
from utils import sequenceutils as su

class FileTypes(Enum):
    INDIVIDUAL = 0
    SUMMARY = 1

class CSVReader():
    def read_individual(self, filename):
        results = {}
        double_line_mode = _get_double_line_mode(filename)

        with open(filename, newline='\n') as file:
            # Dumping the headings row
            next(file)

            # Iterating over each row, reading the content
            for iteration, row in enumerate(file):
                results[iteration] = self._read_individual_result_line(row)

                # If in double_line_mode, skipping the next row
                if double_line_mode:
                    next(file)

        return results

    def read_summary(self, filename):
        freq = {}
        double_line_mode = _get_double_line_mode(filename)

        with open(filename, newline='\n') as file:
            # Dumping the headings row
            next(file)

            # Iterating over each row, reading the content
            for row in file:
                (key, value) = self._read_summary_result_line(row)

                # If the key and value are None, we've reached the end of the file
                if key is None and value is None:
                    return freq
                    
                freq[key] = value

                # If in double_line_mode, skipping the next row
                if double_line_mode:
                    next(file)

        return freq

    def _read_individual_result_line(self, row):
        contents = row.split(',')
        return (int(contents[2]), int(contents[3]), contents[4], contents[5])

    def _read_summary_result_line(self, row):
        contents = row.split(',')

        # The final line of the summary file is a record of the number of failed sequences
        if contents[0] == "Error":
            return (None, None)

        key = (int(contents[3]), int(contents[4]))
        value = int(contents[1])

        return (key, value)

class CSVWriter():
    def __init__(self, extra_nt, local_r=1, append_dt=False, double_line_mode=False):
        self._extra_nt = extra_nt
        self._local_r = local_r
        self._append_dt = append_dt
        self._double_line_mode = double_line_mode

    def write_individual(self, root_name, results, ref):
        file = fu.open_file(root_name, '_individual', 'csv', append_dt=self._append_dt)

        # Initialising string
        row = self._get_individual_header_line()
        file.write(row)

        # Iterating over each result, adding it as a new line
        for index, result in enumerate(results.values()):
            row = self._get_individual_result_line(result, index, ref)
            file.write(row)

        file.close()

    def _get_individual_header_line(self):
        row = "INDEX,HEADER,TYPE,TOP_POS,BOTTOM_POS,SPLIT_SEQ,TOP_LOCAL_SEQ,BOTTOM_LOCAL_SEQ,"

        if self._double_line_mode:
            return row + "SEQUENCE\n"
        else:
            return row + "TOP_SEQUENCE,BOTTOM_SEQUENCE\n"

    def _get_individual_result_line(self, result, index, ref):
        # Initialising string for this line with the result index
        new_row = str(index+1) + ','

        if result is None:
            return new_row + '\n'

        (cleavage_site_t, cleavage_site_b, split, local_seq_t, local_seq_b, header) = result

        # Adding header text
        new_row = new_row + header + ","

        # Adding cleavage type
        new_row = new_row + su.get_type_str(cleavage_site_t, cleavage_site_b, split) + ','

        # Adding top and bottom cleavage positions
        new_row = new_row + str(cleavage_site_t) + ','
        new_row = new_row + str(cleavage_site_b) + ','

        # Adding split sequence boolean
        new_row = new_row + str(split) + ','

        # Adding top and bottom cleavage local sequences
        new_row = new_row + str(local_seq_t) + ','
        new_row = new_row + str(local_seq_b) + ','

        # Adding sequences
        (seq1, seq2) = su.get_sequence_str(ref, cleavage_site_t, cleavage_site_b, split, extra_nt=self._extra_nt)
        if self._double_line_mode:
            new_row = new_row + seq1 + '\n,,,,,,,,' + seq2
        else:
            new_row = new_row + seq1 + ','
            new_row = new_row + seq2

        return new_row + '\n'
        
    def write_summary(self, root_name, freq, ref, error_count):
        file = fu.open_file(root_name, '_summary', 'csv', append_dt=self._append_dt)

        if freq is None:
            file.close()
            return

        total = sum(list(freq.values())) 

        # Initialising string
        row = self._get_summary_header_line()
        file.write(row)

        # Iterating over each result, adding it as a new line
        for cleavage_sites in freq.keys():        
            cleavage_site_t = cleavage_sites[0]
            cleavage_site_b = cleavage_sites[1]
            split = cleavage_sites[2]
            count = freq.get(cleavage_sites)
            event_pc = 100*count/total

            row = self._get_summary_result_line(cleavage_site_t, cleavage_site_b, split, count, event_pc, ref)
            file.write(row)

        # Adding a line for the number of errors
        row = self._get_summary_error_line(error_count)
        file.write(row)

        file.close()

    def _get_summary_header_line(self):
        row = "TYPE,COUNT,EVENT_%,TOP_POS,BOTTOM_POS,SPLIT_SEQ,TOP_LOCAL_SEQ,BOTTOM_LOCAL_SEQ,"

        if self._double_line_mode:
            return row + "SEQUENCE\n"
        else:
            return row + "TOP_SEQUENCE,BOTTOM_SEQUENCE\n"

    def _get_summary_result_line(self, cleavage_site_t, cleavage_site_b, split, count, event_pc, ref):
        if cleavage_site_b is None or cleavage_site_t is None:
            return '\n'

        # Initialising string for this line with the cleavage type, count and event %
        type = su.get_type_str(cleavage_site_t,cleavage_site_b,split)
        new_row = type + ',' + str(count) + ',' + str(event_pc) + ','

        # Adding top and bottom cleavage positions
        new_row = new_row + str(cleavage_site_t) + ','
        new_row = new_row + str(cleavage_site_b) + ','

        # Adding split sequence boolean
        new_row = new_row + str(split) + ','

        # Adding top and bottom cleavage local sequences
        (local_seq_t, local_seq_b) = su.get_local_sequences(ref,cleavage_site_t,cleavage_site_b,local_r=self._local_r)
        new_row = new_row + str(local_seq_t) + ','
        new_row = new_row + str(local_seq_b) + ','

        # Adding sequences
        (seq1, seq2) = su.get_sequence_str(ref, cleavage_site_t, cleavage_site_b, split, extra_nt=self._extra_nt)
        if self._double_line_mode:
            new_row = new_row + seq1 + '\n,,,,,,,,' + seq2
        else:
            new_row = new_row + seq1 + ','
            new_row = new_row + seq2

        return new_row + '\n'

    def _get_summary_error_line(self, count):
        return 'Error,' + str(count) + '\n'

def get_file_type(filename):
    with open(filename, newline='\n') as file:
        # Getting the headings row
        headings = next(file).split(',')

    if headings[0] == "TYPE":
        return FileTypes.SUMMARY
    elif headings[0] == "INDEX":
        return FileTypes.INDIVIDUAL

def _get_double_line_mode(filename):
        with open(filename, newline='\n') as csvfile:
            headings = next(csvfile).split(',')
            return "TOP" not in headings[6]
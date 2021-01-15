import csv

from utils import fileutils as fu
from utils import sequenceutils as su

class CSVReader():
    def read_individual(self, filename):
        results = {}
        double_line_mode = self._get_double_line_mode(filename)

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
                

    def _get_double_line_mode(self, filename):
        with open(filename, newline='\n') as csvfile:
            headings = next(csvfile).split(',')
            return "TOP" not in headings[6]

    def _read_individual_result_line(self, row):
        contents = row.split(',')
        return (int(contents[2]), int(contents[3]), contents[4], contents[5])


class CSVWriter():
    def __init__(self, extra_nt, local_r=1, append_dt=False, double_line_mode=False):
        self._extra_nt = extra_nt
        self._local_r = local_r
        self._append_dt = append_dt
        self._double_line_mode = double_line_mode

    def write_individual(self, root_name, results, ref):
        file = fu.open_file(root_name, '_individual', 'csv', append_dt=self._append_dt)

        # Initialising string
        str = self._get_individual_header_line()
        file.write(str)

        # Iterating over each result, adding it as a new line
        for index, result in enumerate(results.values()):
            str = self._get_individual_result_line(result, index, ref)
            file.write(str)

        file.close()

    def _get_individual_header_line(self):
        str = "INDEX,TYPE,TOP_POS,BOTTOM_POS,TOP_LOCAL_SEQ,BOTTOM_LOCAL_SEQ,"

        if self._double_line_mode:
            return str + "SEQUENCE\n"
        else:
            return str + "TOP_SEQUENCE,BOTTOM_SEQUENCE\n"

    def _get_individual_result_line(self, result, index, ref):
        # Initialising string for this line with the result index
        new_str = str(index+1) + ','

        if result is None:
            return new_str + '\n'

        (cleavage_site_t, cleavage_site_b, local_seq_t, local_seq_b) = result

        # Adding cleavage type
        new_str = new_str + su.get_type_str(cleavage_site_t, cleavage_site_b) + ','

        # Adding top and bottom cleavage positions
        new_str = new_str + str(cleavage_site_t) + ','
        new_str = new_str + str(cleavage_site_b) + ','

        # Adding top and bottom cleavage local sequences
        new_str = new_str + str(local_seq_t) + ','
        new_str = new_str + str(local_seq_b) + ','

        # Adding sequences
        (seq1, seq2) = su.get_sequence_str(ref, cleavage_site_t, cleavage_site_b, extra_nt=self._extra_nt)
        if self._double_line_mode:
            new_str = new_str + seq1 + '\n,,,,,,' + seq2
        else:
            new_str = new_str + seq1 + ','
            new_str = new_str + seq2

        return new_str + '\n'
        
    def write_summary(self, root_name, freq, ref, error_count):
        file = fu.open_file(root_name, '_summary', 'csv', append_dt=self._append_dt)

        if freq is None:
            file.close()
            return

        total = sum(list(freq.values())) + error_count

        # Initialising string
        str = self._get_summary_header_line()
        file.write(str)

        # Iterating over each result, adding it as a new line
        for cleavage_sites in freq.keys():        
            cleavage_site_t = cleavage_sites[0]
            cleavage_site_b = cleavage_sites[1]
            count = freq.get(cleavage_sites)
            event_pc = 100*count/total

            str = self._get_summary_result_line(cleavage_site_t, cleavage_site_b, count, event_pc, ref)
            file.write(str)

        # Adding a line for the number of errors
        error_pc = 100*error_count/total
        str = self._get_summary_error_line(error_count, error_pc)
        file.write(str)

        file.close()

    def _get_summary_header_line(self):
        str = "TYPE,COUNT,EVENT_%,TOP_POS,BOTTOM_POS,TOP_LOCAL_SEQ,BOTTOM_LOCAL_SEQ,"

        if self._double_line_mode:
            return str + "SEQUENCE\n"
        else:
            return str + "TOP_SEQUENCE,BOTTOM_SEQUENCE\n"

    def _get_summary_result_line(self, cleavage_site_t, cleavage_site_b, count, event_pc, ref):
        if cleavage_site_b is None or cleavage_site_t is None:
            return '\n'

        # Initialising string for this line with the cleavage type, count and event %
        type = su.get_type_str(cleavage_site_t,cleavage_site_b)
        new_str = type + ',' + str(count) + ',' + str(event_pc) + ','

        # Adding top and bottom cleavage positions
        new_str = new_str + str(cleavage_site_t) + ','
        new_str = new_str + str(cleavage_site_b) + ','

        # Adding top and bottom cleavage local sequences
        (local_seq_t, local_seq_b) = su.get_local_sequences(ref,cleavage_site_t,cleavage_site_b,local_r=self._local_r)
        new_str = new_str + str(local_seq_t) + ','
        new_str = new_str + str(local_seq_b) + ','

        # Adding sequences
        (seq1, seq2) = su.get_sequence_str(ref, cleavage_site_t, cleavage_site_b, extra_nt=self._extra_nt)
        if self._double_line_mode:
            new_str = new_str + seq1 + '\n,,,,,,,' + seq2
        else:
            new_str = new_str + seq1 + ','
            new_str = new_str + seq2

        return new_str + '\n'

    def _get_summary_error_line(self, count, event_pc):
        return 'Error,' + str(count) + ',' + str(event_pc) + '\n'

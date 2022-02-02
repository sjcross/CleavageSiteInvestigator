import datetime as dt
import os
import re

from Bio.Seq import Seq
from rich import print


class FileReader():
    def __init__(self, verbose=True):
        self._verbose = verbose

    def read_sequence(self, path, repeat_filter=""):
        name = os.path.basename(path)
        
        if self._verbose:
            print("    Loading file \"%s\"" % name)

        # Get extension and run appropriate reader
        rootname, ext = os.path.splitext(path)

        if ext == ".ab1":
            return self._read_ab1(path)
        elif ext == ".dna":
            return self._read_dna(path)
        elif ext == ".fa" or ext == ".fasta":
            return self._read_fasta(path, repeat_filter=repeat_filter)
        elif ext == ".txt":
            return self._read_txt(path)
        elif ext == ".seq":
            return self._read_seq(path)

    def get_verbose(self):
        return self._verbose

    def set_verbose(self, verbose):
        self._verbose = verbose

    def _read_ab1(self, path, repeat_filter=""):
        if repeat_filter != "":
            print("            Repeat filtering not possible for .ab1 files")

        if self._verbose:
            print("        Reading as \".ab1\" format")

        file = open(path, "rb")
        full_text = file.read()

        pattern = re.compile("term{1}([ACGT]+)".encode())
        sequence_string = pattern.findall(full_text)[0].decode()

        sequence_string = self._remove_repeats(sequence_string)
        sequence_string = convert_to_upper_case(sequence_string)

        return ([(Seq(sequence_string),"")], 1, 1)

    def _read_dna(self, path, repeat_filter=""):
        if repeat_filter != "":
            print("            Repeat filtering not possible for .dna files")

        if self._verbose:
            print("        Reading as \".dna\" format")

        file = open(path, "rb")
        full_text = file.read()

        pattern = re.compile("([ACGTacgt]+)".encode())
        instances = pattern.findall(full_text)

        sequence_string = self._get_longest_sequence(instances).decode()
        sequence_string = convert_to_upper_case(sequence_string)

        return ([(Seq(sequence_string),"")], 1, 1)

    def _read_fasta(self, path, repeat_filter=""):
        if self._verbose:
            print("        Reading as \".fasta\" format")

        # Initialising counters for accepted and rejected sequences based on the number of repeats
        n_acc = 0
        n_rej = 0

        # Initialising the sequence list
        seqs = []

        file = open(path, "r")
        full_text = file.read()

        pattern = re.compile("(>.+\n)([ACGTacgt\n]+)")
        instances = pattern.findall(full_text)

        # Creating header pattern
        header_pattern = re.compile(">.+_[\d.]+_\d+_(\d+)_\d+\n")

        for instance in instances:
            header = instance[0]
            seq = instance[1]

            # Checking number of repeats
            if repeat_filter != "":
                # Matching C3P0a format in header
                header_instances = header_pattern.findall(header)

                if len(header_instances) == 1:
                    repeat_filter_curr = repeat_filter.replace("x", header_instances[0])

                    # If repeat filter fails for current number of repeats, skip this sequence and increase n_rej by 1
                    if not eval(repeat_filter_curr):                        
                        n_rej = n_rej + 1
                        continue

            # Removing linebreaks in header and sequence
            header = header.replace("\n", "")
            header = header.replace(">", "")
            header = header.replace(",", ";")
            seq = seq.replace("\n", "")
            seq = convert_to_upper_case(seq)
            seqs.append((Seq(seq),header))

            n_acc = n_acc + 1

        if self._verbose:
            print("        Loaded %i sequence(s)" % len(seqs))

        return (seqs, n_acc, n_rej)

    def _read_seq(self, path, repeat_filter=""):
        if repeat_filter != "":
            print("            Repeat filtering not possible for .seq files")

        if self._verbose:
            print("        Reading as \".seq\" format")

        file = open(path, "r")
        sequence_string = file.read()

        # Removing linebreaks
        sequence_string = sequence_string.replace("\n", "")
        sequence_string = convert_to_upper_case(sequence_string)

        return ([(Seq(sequence_string),"")], 1, 1)

    def _read_txt(self, path, repeat_filter=""):
        if repeat_filter != "":
            print("            Repeat filtering not possible for .txt files")

        if self._verbose:
            print("        Reading as \".txt\" format")

        file = open(path, "r")
        return ([(Seq(convert_to_upper_case(file.read())),"")], 1, 1)

    def _get_longest_sequence(self, instances):
        max_len = 0
        max_pos = 0
        for i, instance in enumerate(instances):
            if len(instance) > max_len:
                max_len = len(instance)
                max_pos = i

        return instances[max_pos]

    def _remove_repeats(self, seq_string, n=10):
        # Getting intro sequence
        intro = seq_string[0:n]

        # Find all instances of this string
        pattern = re.compile(intro)
        instances = pattern.split(seq_string)

        if self._verbose:
            print(
                "        Found %i instances of repeated sequence (loading first)"
                % (len(instances) - 1))

        return instances[1]


def convert_to_upper_case(seq_string):
    seq_string = seq_string.replace("g", "G")
    seq_string = seq_string.replace("c", "C")
    seq_string = seq_string.replace("a", "A")
    seq_string = seq_string.replace("t", "T")

    return seq_string

def open_file(root_name, suffix, extension, append_dt):
    datetime_str = dt.datetime.now().strftime("_%Y-%m-%d_%H-%M-%S") if append_dt else ""
    outname = root_name+ suffix + datetime_str + '.' + extension

    try:
        return open(outname, "w", encoding="utf-8")
    except:
        outname_orig = outname
        datetime_str = dt.datetime.now().strftime("_%Y-%m-%d_%H-%M-%S")
        outname = root_name+ suffix + datetime_str + '.' + extension

        print('WARNING: File "%s" unavailable for writing, storing to "%s" instead' % (outname_orig,outname))

        return open(outname, "w", encoding="utf-8")

        
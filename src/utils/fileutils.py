import os
import re

from Bio.Seq import Seq


class FileReader():
    def __init__(self, verbose=True):
        self._verbose = verbose

    def read_sequence(self, path):
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
            return self._read_fasta(path)
        elif ext == ".txt":
            return self._read_txt(path)
        elif ext == ".seq":
            return self._read_seq(path)

    def get_verbose(self):
        return self._verbose

    def set_verbose(self, verbose):
        self._verbose = verbose

    def _read_ab1(self, path):
        if self._verbose:
            print("        Reading as \".ab1\" format")

        file = open(path, "rb")
        full_text = file.read()

        pattern = re.compile("term{1}([ACGT]+)".encode())
        sequence_string = pattern.findall(full_text)[0].decode()

        sequence_string = self._remove_repeats(sequence_string)

        return [Seq(sequence_string)]

    def _read_dna(self, path):
        if self._verbose:
            print("        Reading as \".dna\" format")

        file = open(path, "rb")
        full_text = file.read()

        pattern = re.compile("([ACGTacgt]+)".encode())
        instances = pattern.findall(full_text)

        sequence_string = self._get_longest_sequence(instances).decode()

        return [Seq(sequence_string)]

    def _read_fasta(self, path):
        if self._verbose:
            print("        Reading as \".fasta\" format")

        file = open(path, "r")
        full_text = file.read()
        
        pattern = re.compile(">.+\n([ACGTacgt\n]+)")
        instances = pattern.findall(full_text)

        for i, instance in enumerate(instances):
            # Removing linebreaks
            instance = instance.replace("\n", "")
            instances[i] = Seq(instance)

        if self._verbose:
            print("        Loaded %i sequence(s)" % len(instances))

        return instances

    def _read_seq(self, path):
        if self._verbose:
            print("        Reading as \".seq\" format")

        file = open(path, "r")
        sequence_string = file.read()

        # Removing linebreaks
        sequence_string = sequence_string.replace("\n", "")

        return [Seq(sequence_string)]

    def _read_txt(self, path):
        if self._verbose:
            print("        Reading as \".txt\" format")

        file = open(path, "r")
        return [Seq(file.read())]

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
            print("        Found %i instances of repeated sequence (loading first)" % (
                len(instances)-1))

        return instances[1]

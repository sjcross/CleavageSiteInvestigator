import os

from utils import fileutils as fu
from utils import plotutils as pu
from utils import reportutils as ru

verbose = False
ref_path = "C:\\Users\\steph\\Desktop\\Oscar\\pUC19.fa"
test_path = "C:\\Users\\steph\\Desktop\\Oscar\\test.fasta"

root_name = os.path.splitext(test_path)[0]

filereader = fu.FileReader(verbose=verbose)
ref = filereader.read_sequence(ref_path)[0][0]

freq_full = {}
freq_full[(450,480)] = 21
freq_full[(450,460)] = 68
freq_full[(450,453)] = 9
freq_full[(452,456)] = 43
freq_full[(453,462)] = 29
freq_full[(455,468)] = 31
freq_full[(458,472)] = 25
freq_full[(448,464)] = 72

freq_full[(421,410)] = 59
freq_full[(421,412)] = 10
freq_full[(419,407)] = 19
freq_full[(419,405)] = 43
freq_full[(425,405)] = 43

freq_full = ru.sort_results(freq_full)

pos_min_zb = 0
pos_max_zb = len(ref)
pos_min_zb = 400
pos_max_zb = 500
pu.plotEventDistribution(root_name, ref, freq_full, pos_min_zb, pos_max_zb)
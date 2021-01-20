import os

from utils import fileutils as fu
from utils import heatmapwriter as hmw
from utils import reportutils as ru

verbose = False
ref_path = "C:\\Users\\steph\\OneDrive\\Desktop\\Oscar\\pUC19.fa"
out_path = "C:\\Users\\steph\\OneDrive\\Desktop\\Oscar\\test.svg"

filereader = fu.FileReader(verbose=verbose)
ref = filereader.read_sequence(ref_path)[0][0]

freq_full = {}
freq_full[(450,440)] = 21
freq_full[(450,459)] = 68
freq_full[(450,453)] = 9
freq_full[(452,456)] = 43
freq_full[(453,442)] = 29
freq_full[(455,448)] = 31
freq_full[(458,452)] = 25
freq_full[(448,444)] = 72

freq_full[(461,450)] = 59
freq_full[(461,452)] = 10
freq_full[(459,447)] = 19
freq_full[(459,445)] = 43
freq_full[(445,445)] = 43

freq = ru.sort_results(freq_full)

writer = hmw.HeatMapWriter(grid_label_opts=(True,12,"gray",1,10), event_colourmap="plasma")
writer.write_map(out_path, freq, ref=ref, pos_ranges=(440,460,440,460), append_dt=True)
from utils import fileutils as fu
from utils import eventmapwriter as emw
from utils import heatmapwriter as hmw
from utils import reportutils as ru

verbose = False
ref_path = "C:\\Users\\steph\\OneDrive\\Desktop\\Oscar\\pUC19.fa"
out_path = "C:\\Users\\steph\\OneDrive\\Desktop\\Oscar\\"

filereader = fu.FileReader(verbose=verbose)
ref = filereader.read_sequence(ref_path)[0][0]

freq_full[(460,450)] = 59
freq_full[(461,452)] = 10
freq_full[(459,447)] = 8
freq_full[(459,445)] = 23
freq_full[(465,445)] = 43

freq = ru.sort_results(freq_full)

writer = emw.EventMapWriter(dna_opts=(emw.DNA_MODE.SEQUENCE,8,"black"), grid_opts=(True,1,"lightgray",5), grid_label_opts=(True,12,"gray",20,10), event_opts=(2, "plasma"))
writer.write_map(out_path+"Example eventmap.svg", freq_full, ref=ref, pos_range=(400,500))

writer = hmw.HeatMapWriter(grid_opts=(True,1,"lightgray",1), grid_label_opts=(True,12,"gray",1,10), event_colourmap="plasma")
writer.write_map(out_path+"Example heatmap.svg", freq, ref=ref, pos_ranges=(440,460,440,460))

from utils import fileutils as fu
from utils import eventmapwriter as emw
from utils import heatmapwritercsv as hmwc
from utils import heatmapwritersvg as hmws
from utils import reportutils as ru

verbose = False
ref_path = "C:\\Users\\steph\\OneDrive\\Desktop\\Oscar\\pUC19.fa"
out_path = "C:\\Users\\steph\\OneDrive\\Desktop\\Oscar\\"

filereader = fu.FileReader(verbose=verbose)
ref = filereader.read_sequence(ref_path)[0][0]

freq_full = {}
freq_full[(450,440)] = 21
freq_full[(450,459)] = 68
freq_full[(450,453)] = 9
freq_full[(452,456)] = 43
freq_full[(453,442)] = 29
freq_full[(455,448)] = 31
freq_full[(458,452)] = 4
freq_full[(448,444)] = 72

freq_full[(460,450)] = 59
freq_full[(461,452)] = 10
freq_full[(459,447)] = 8
freq_full[(459,445)] = 23
freq_full[(465,445)] = 43

freq = ru.sort_results(freq_full)

# writer = emw.EventMapWriter(dna_opts=(emw.DNA_MODE.SEQUENCE,8,"black"), grid_opts=(True,1,"lightgray",5), grid_label_opts=(True,12,"gray",20,10), event_opts=(2, "plasma"))
# writer.write_map(out_path+"Example eventmap.svg", freq_full, ref=ref, pos_range=(400,500))

writer = hmws.HeatMapWriterSVG(grid_opts=(True,1,"lightgray",1), grid_label_opts=(True,12,"gray",1,10), event_colourmap="plasma")
writer.write_map(out_path+"Example heatmapb.svg", freq, ref, (440,460,440,460), False)

writer = hmwc.HeatMapWriterCSV()
writer.write_map(out_path+"Example heatmapb.csv", freq, ref, (440,460,440,460), False)
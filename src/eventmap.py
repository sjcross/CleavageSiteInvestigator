from utils import svgutils as svu

filename = "C:\\Users\\steph\\Desktop\\Oscar\\mix_individual.csv"
ref_path = "C:\\Users\\steph\\Desktop\\Oscar\\pUC19.fa" # Can be "" (default), but won't allow sequence to be rendered

writer = svu.EventMapWriter(grid_opts=(True,1,"lightgray",5), grid_label_opts=(True,12,"lightgray",10), dna_opts=(True,8,"black"))
writer.write_event_map_from_file(filename, ref_path=ref_path, pos_range=(400,500), append_dt=True)

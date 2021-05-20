# $data_path = "C:\Users\steph\Desktop\Oscar\mix_individual.csv"
# $ref_path = "C:\Users\steph\Desktop\Oscar\pUC19.fa" # Can be "" (default), but won't allow sequence to be rendered

# # Creating a full-length sequence, showing the DNA as lines
# $out_path = "C:\Users\steph\Desktop\Oscar\eventmap_full.svg"
# python .\src\eventmap.py -d $data_path -o $out_path -r $ref_path -dm 'line'

# Creating a zoomed-in sequence from 400-500, showing DNA as text sequence
# $out_path = "C:\Users\steph\Desktop\Oscar\eventmap_subset_400-500.svg"
# python .\src\eventmap.py -d $data_path -o $out_path -r $ref_path -pr 400 500 -dm 'seq' -ds 8 -gi 5 -gli 20 -eo 1 -er 10 60 -cli 10 -ad

# $data_path = "D:\Stephen\People\Oscar\2020-02-20 Sequencing analysis\2021-02-23 Circularity test data\X data sets\A\A_summary.csv"
# $ref_path = "D:\Stephen\People\Oscar\2020-02-20 Sequencing analysis\2021-02-23 Circularity test data\X data sets\A\puc19_2.fa" # Can be "" (default), but won't allow sequence to be rendered

# Creating a full-length sequence, showing the DNA as lines
# $out_path = "D:\Stephen\People\Oscar\2020-02-20 Sequencing analysis\2021-02-23 Circularity test data\X data sets\A\A_split.svg"
# python .\src\eventmap.py -d $data_path -o $out_path -r $ref_path -dm 'line' -emis 1 -eo 0.5 -eso 2

# Creating a zoomed-in sequence from 400-500, showing DNA as text sequence
# $out_path = "C:\Users\steph\Desktop\Oscar\eventmap_subset_400-500.svg"
# python .\src\eventmap.py -d $data_path -o $out_path -r $ref_path -pr 400 500 -dm 'seq' -ds 8 -gi 5 -gli 20 -dc "#16C3D6"

$data_path = "C:\Users\steph\Desktop\Oscar\McrBC63_5min_summary.csv"
$ref_path = "C:\Users\steph\Desktop\Oscar\pUC19.fa" # Can be "" (default), but won't allow sequence to be rendered
$out_path = "C:\Users\steph\Desktop\Oscar\McrBC63_5min_summary_eventmap.svg"

python .\src\eventmap.py -d $data_path -o $out_path -r $ref_path -pr 2000 2050 -gi 10 -dm 'seq' -ds 12 -hlzv 'hide' -hpbg 20 -eorv 'show' -ho 4
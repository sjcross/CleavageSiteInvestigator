$data_path = "D:\Stephen\People\Oscar\2020-02-20 Sequencing analysis\2020-04-28 New files\mix_individual.csv"
$ref_path = "D:\Stephen\People\Oscar\2020-02-20 Sequencing analysis\2020-04-28 New files\pUC19.fa" # Can be "" (default), but won't allow sequence to be rendered

# Creating a full-length sequence, showing the DNA as lines
$out_path = "D:\Stephen\People\Oscar\2020-02-20 Sequencing analysis\2020-04-28 New files\eventmap_full.svg"
python .\src\eventmap.py -d $data_path -o $out_path -r $ref_path -dm 'line' -ad

# Creating a zoomed-in sequence from 400-500, showing DNA as text sequence
$out_path = "D:\Stephen\People\Oscar\2020-02-20 Sequencing analysis\2020-04-28 New files\eventmap_subset_400-500.svg"
python .\src\eventmap.py -d $data_path -o $out_path -r $ref_path -pr 400 500 -dm 'seq' -ds 8 -gi 5 -gli 20 -ad -dc "#16C3D6"
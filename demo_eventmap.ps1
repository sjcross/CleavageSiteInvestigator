$data_path = "C:\\Users\\steph\\Desktop\\Oscar\\mix_individual.csv"
$ref_path = "C:\\Users\\steph\\Desktop\\Oscar\\pUC19.fa" # Can be "" (default), but won't allow sequence to be rendered

# Creating a full-length sequence, showing the DNA as lines
$out_path = "C:\\Users\\steph\\Desktop\\Oscar\\eventmap_full.svg"
python .\src\eventmap.py -d $data_path -o $out_path -r $ref_path -dm 'line' -ad

# Creating a zoomed-in sequence from 400-500, showing DNA as text sequence
$out_path = "C:\\Users\\steph\\Desktop\\Oscar\\eventmap_subset_400-500.svg"
python .\src\eventmap.py -d $data_path -o $out_path -r $ref_path -pr 400 500 -dm 'seq' -ds 8 -gi 5 -gli 20 -ad -dc "#16C3D6"
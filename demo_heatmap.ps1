$data_path = "C:\\Users\\steph\\Desktop\\Oscar\\mix_individual.csv"
$ref_path = "C:\\Users\\steph\\Desktop\\Oscar\\pUC19.fa" # Can be "" (default), but won't allow sequence to be rendered

# Creating a full-length sequence, showing the DNA as lines
# $out_path = "C:\\Users\\steph\\Desktop\\Oscar\\heatmap_full.svg"
# python .\src\heatmap.py -d $data_path -o $out_path -r $ref_path -ad

# Creating a zoomed-in sequence from 400-500, showing DNA as text sequence
# $out_path = "C:\\Users\\steph\\Desktop\\Oscar\\heatmap_subset.svg"
# python .\src\heatmapsvg.py -d $data_path -o $out_path -r $ref_path -pr 420 440 420 440 -bs 2 -gi 1 -gc "gray" -gli 1 -glc "gray" -ad -ec "plasma" 

$out_path = "C:\\Users\\steph\\Desktop\\Oscar\\heatmap_subset2.csv"
python .\src\heatmapcsv.py -d $data_path -o $out_path -r $ref_path -pr 420 440 420 440

$out_path = "C:\\Users\\steph\\Desktop\\Oscar\\heatmap_subset2.svg"
python .\src\heatmapsvg.py -d $data_path -o $out_path -r $ref_path -pr 420 440 420 440
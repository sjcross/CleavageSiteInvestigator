# For neatness setting these values as variables
$cass_path = "\\LUCA\People\Oscar\2020-02-20 Sequencing analysis\2020-04-28 New files\Chloramphenicol Cassette overhang.fa"
$ref_path = "\\LUCA\People\Oscar\2020-02-20 Sequencing analysis\2020-04-28 New files\pUC19.fa"

# Running a single file
python .\src\csi.py -c $cass_path -r $ref_path -t "\\LUCA\People\Oscar\2020-02-20 Sequencing analysis\2020-04-28 New files\mix.fasta" -en 3 -pr -we -wo -wi -ws -whsf

# # Running a folder full of files (will only process .fa or .fasta files)
# $test_path = "D:\People\CSI\2020-06-04 Mix files"
# Get-ChildItem $test_path -include ('*.fa', '*.fasta') -recurse | 
# Foreach-Object {
#     python .\src\main.py -c $cass_path -r $ref_path -t $_.FullName -sr -en 3
# }

# $cass_path = "C:\Users\steph\Desktop\Oscar\Chloramphenicol Cassette overhang.fa"
# $ref_path = "C:\Users\steph\Desktop\Oscar\pUC19.fa"

# # Running a single file
# python .\src\csi.py -c $cass_path -r $ref_path -t "C:\Users\steph\Desktop\Oscar\mix.fasta" -en 3 -rf "x>=3" -wi -we

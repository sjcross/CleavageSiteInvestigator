# For neatness setting these values as variables
$cass_path = "C:\Users\steph\Desktop\Oscar\Chloramphenicol Cassette overhang.fa"
$ref_path = "C:\Users\steph\Desktop\Oscar\pUC19.fa"

# Running a single file
python .\src\main.py -c $cass_path -r $ref_path -t "C:\Users\steph\Desktop\Oscar\mix.fasta" -sr -sc -en 3 -v -rf "x>=3"

# # Running a folder full of files (will only process .fa or .fasta files)
# $test_path = "D:\People\CSI\2020-06-04 Mix files"
# Get-ChildItem $test_path -include ('*.fa', '*.fasta') -recurse | 
# Foreach-Object {
#     python .\src\main.py -c $cass_path -r $ref_path -t $_.FullName -sr -en 3
# }

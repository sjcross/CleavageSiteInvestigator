# For neatness setting these values as variables
$cass_path = "Z:\Oscar\2020-02-20 Sequencing analysis\2020-11-25 Problem file\Splint1TA.fa"
$ref_path = "Z:\Oscar\2020-02-20 Sequencing analysis\2020-11-25 Problem file\pUC19.fa"

# Running a single file
python .\src\main.py -c $cass_path -r $ref_path -t "Z:\Oscar\2020-02-20 Sequencing analysis\2020-11-25 Problem file\SspITA.fasta" -sr -en 3 -v -rf "x>=3"

# # Running a folder full of files (will only process .fa or .fasta files)
# $test_path = "D:\People\CSI\2020-06-04 Mix files"
# Get-ChildItem $test_path -include ('*.fa', '*.fasta') -recurse | 
# Foreach-Object {
#     python .\src\main.py -c $cass_path -r $ref_path -t $_.FullName -sr -en 3
# }

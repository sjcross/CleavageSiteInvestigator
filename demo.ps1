# For neatness setting these values as variables
$cass_path = "D:\People\CSI\2020-04-28 New files\Chloramphenicol Cassette overhang.fa"
$ref_path = "D:\People\CSI\2020-04-28 New files\pUC19.fa"

# Running a single file
python .\src\main.py -c $cass_path -r $ref_path -t "D:\People\CSI\2020-06-04 Mix files\PstI_R2C2_Consensus.fasta" -sr -en 3

# Running a folder full of files (will only process .fa or .fasta files)
$test_path = "D:\People\CSI\2020-06-04 Mix files"
Get-ChildItem $test_path -include ('*.fa', '*.fasta') -recurse | 
Foreach-Object {
    python .\src\main.py -c $cass_path -r $ref_path -t $_.FullName -sr -en 3
}
$cass_path = "C:\\Users\\steph\\Desktop\\Oscar\\testcass.fa"
$ref_path = "C:\\Users\\steph\\Desktop\\Oscar\\testref.fa" # Can be "" (default), but won't allow sequence to be rendered
$test_path = "C:\\Users\\steph\\Desktop\\Oscar\\testtest.fa"

python .\src\csi.py -c $cass_path -r $ref_path -t $test_path -nb 4 -v
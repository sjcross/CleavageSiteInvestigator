# Imports
from filereader import FileReader

# Setting some parameters
ref_seq_path = "C:\\Users\\steph\\OneDrive\\Desktop\\Mark Szczelkun\\"
ref_seq_name = "pRMA03+L2L2.dna"
cass_seq_path = "C:\\Users\\steph\\OneDrive\\Desktop\\Mark Szczelkun\\"
cass_seq_name = "CAT cassette as amplified by RA101 & 102 from pACYC184.dna"
test_seq_path = "C:\\Users\\steph\\OneDrive\\Desktop\\Mark Szczelkun\\Sequencing files\\Run250808-04\\"
test_seq_name = "02_007.ab1"

cass_end_length = 20
cass_end_offset = 30

# Creating FileHandler object
filereader = FileReader()
filereader.set_verbose(True)

# Loading reference, cassette and test sequences
ref_seq = filereader.read_sequence(ref_seq_path, ref_seq_name)
# print("Reference sequence:\n%s" % ref_seq.get_sequence_string())
print("\n")

cass_seq = filereader.read_sequence(cass_seq_path, cass_seq_name)
# print("Cassette sequence:\n%s" % cass_seq.get_sequence_string())
print("\n")

test_seq = filereader.read_sequence(test_seq_path, test_seq_name)
# print("Test sequence:\n%s" % test_seq.get_sequence_string())
print("\n")

## Finding cassette ends in test sequence
# Finding start sequence
cass_start_seq = cass_seq.get_start_sequence(cass_end_length,offset=cass_end_offset)
# cass_start_r_seq = cass_start_seq.get_reverse_sequence()
# cass_start_c_seq = cass_start_seq.get_complement_sequence()
cass_start_rc_seq = cass_start_seq.get_reverse_complement_sequence()

print("Start sequence: %s (%i matches)" % (cass_start_seq.get_sequence_string(), test_seq.find_sequence_match(cass_start_seq)))
# print("Start sequence (R): %s (%i matches)" % (cass_start_r_seq.get_sequence_string(), test_seq.find_sequence_match(cass_start_r_seq)))
# print("Start sequence (C): %s (%i matches)" % (cass_start_c_seq.get_sequence_string(), test_seq.find_sequence_match(cass_start_c_seq)))
print("Start sequence (RC): %s (%i matches)" % (cass_start_rc_seq.get_sequence_string(), test_seq.find_sequence_match(cass_start_rc_seq)))

cass_end_seq = cass_seq.get_end_sequence(cass_end_length,offset=cass_end_offset)
# cass_end_r_seq = cass_end_seq.get_reverse_sequence()
# cass_end_c_seq = cass_end_seq.get_complement_sequence()
cass_end_rc_seq = cass_end_seq.get_reverse_complement_sequence()

print("End sequence: %s (%i matches)" % (cass_end_seq.get_sequence_string(), test_seq.find_sequence_match(cass_end_seq)))
# print("End sequence (R): %s (%i matches)" % (cass_end_r_seq.get_sequence_string(), test_seq.find_sequence_match(cass_end_r_seq)))
# print("End sequence (C): %s (%i matches)" % (cass_end_c_seq.get_sequence_string(), test_seq.find_sequence_match(cass_end_c_seq)))
print("End sequence (RC): %s (%i matches)" % (cass_end_rc_seq.get_sequence_string(), test_seq.find_sequence_match(cass_end_rc_seq)))

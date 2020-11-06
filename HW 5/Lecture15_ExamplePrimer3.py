import argparse, subprocess

parser = argparse.ArgumentParser(description = 'This program picks PCR primers for any sequence')
parser.add_argument("sequence", type=str, help="Sequence to pick primers")
args = parser.parse_args()

with open('primer3_command_file_python.txt', 'w') as f:
	print('SEQUENCE_ID=example', file = f)
	print('SEQUENCE_TEMPLATE=' + args.sequence, file = f)
	print('PRIMER_THERMODYNAMIC_PARAMETERS_PATH=/Users/vicenteblat/anaconda3/share/primer3/primer3_config/', file = f)
	print('=', file = f)

#subprocess.call(['primer3_core', 'primer3_command_file_python.txt'], stdout = open('primer3_output.txt', 'w'))
output = subprocess.run(['primer3_core', 'primer3_command_file_python.txt'], capture_output = True)

print(output)

#with open('primer3_output.txt') as f:
#	for line in f:
#		blahblahbl

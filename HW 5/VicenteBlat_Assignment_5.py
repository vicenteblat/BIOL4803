import argparse, subprocess, pyfaidx

parser = argparse.ArgumentParser(description = 'This program runs primer3 to identify the best two primers to amplify a given sequence')
parser.add_argument('fasta_file', type=str, help='Enter a genome file containing DNA sequence')
parser.add_argument('chromosome', type=str, help='Enter the chromosome you want amplified')
parser.add_argument('position', type=str, help='Enter the position you want amplified')
args = parser.parse_args()

fasta = pyfaidx.Fasta(args.fasta_file)
sequence = str(fasta[args.chromosome][int(args.position) - 500:int(args.position) + 500])

with open('primer3_commands.txt', 'w') as f:
    print('SEQUENCE_ID=HW5', file = f)
    print('SEQUENCE_TEMPLATE=' + sequence, file = f)
    print('PRIMER_THERMODYNAMIC_PARAMETERS_PATH=/Users/vicenteblat/anaconda3/share/primer3/primer3_config/', file = f)
    print('PRIMER_PRODUCT_SIZE_RANGE=600-800', file = f)
    print('SEQUENCE_EXCLUDED_REGION=400,200', file = f)
    print('=', file = f)

o = open('output.txt', 'w')
output = subprocess.run(['primer3_core', 'primer3_commands.txt'], stdout=o)

file = open('output.txt', 'r')
for line in file:
    if line.find('PRIMER_LEFT_0_SEQUENCE') != -1:
        primer = line.split('=')[1]
        print('Left primer = ' + primer)
    elif line.find('PRIMER_RIGHT_0_SEQUENCE') != -1:
        primer = line.split('=')[1]
        print('Right primer = ' + primer)

file.close()

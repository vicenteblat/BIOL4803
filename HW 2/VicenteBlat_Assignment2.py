import sys
import textwrap

def valid_DNA_sequence(DNA):
    #This function takes in a string containing possible DNA sequence and returns True if all nucleotides are valid (A,a,C,c,G,g,T,t) and False if any nucleotide is invalid
    while True:
        for nucleotide in DNA:
            if nucleotide not in 'AaGgCcTt':
                b = False
                break
            else:
                b = True
        break
    return b

def print_DNA_sequence(DNA, mode):
    # This function takes in a string containing valid DNA sequence and a mode that specifies the output format. Prints info to a screen, maximum characters per line
    for i in range(0,3):
        #This loop is used to create the three possible frames on the forward strand
        #Create DNA sequence for current frame. Ensure it is divisible by 3
        if len(DNA) % 3 == 0:
            end = len(DNA) - 1
        elif len(DNA) % 3 == 1:
            end = len(DNA) - 2
        else:
            end = len(DNA) - 3
        DNA_sequence = DNA[i:end+1]
        # if len(DNA_sequence) < 3:
        #     break
        translated_sequence = translate(DNA_sequence, mode)
        # Print out direction (5' to 3') and frame to screen
        print("5' to 3' Frame: " + str(i))
        # Print out translated_sequence. If nucleotide sequence mode selected, print nucleotide sequence and amino acid sequence, 60 nucleotides per line until entire sequence is printed out
        if translated_sequence and mode == 'DNA':
            if i > 0:
                DNA_sequence = DNA_sequence[:i-3]
            while DNA_sequence:
                print(DNA_sequence[:60])
                print(translated_sequence[:60])
                DNA_sequence = DNA_sequence[60:]
                translated_sequence = translated_sequence[60:]
        else:
            print(translated_sequence)


    rev_DNA_sequence = reverse_complement(DNA)
    for j in range(0,3):
        #This loop is used to create the three possible frames on the reverse strand
        #Create DNA sequence for current frame. Ensure it is divisible by 3
        if len(rev_DNA_sequence) % 3 == 0:
            end = len(rev_DNA_sequence) - 1
        elif len(rev_DNA_sequence) % 3 == 1:
            end = len(rev_DNA_sequence) - 2
        else:
            end = len(rev_DNA_sequence) - 3
        DNA_sequence = rev_DNA_sequence[j:end+1]
        # if len(DNA_sequence) < 3:
        #     break
        translated_sequence = translate(DNA_sequence, mode)
        # Print out direction (3' to 5') and frame to screen
        print("3' to 5' Frame: " + str(j))
        # Print out translated_sequence. If nucleotide sequence mode selected, print nucleotide sequence and amino acid sequence, 60 nucleotides per line until entire sequence is printed out
        if translated_sequence and mode == 'DNA':
            if j > 0:
                DNA_sequence = DNA_sequence[:j-3]
            while DNA_sequence:
                print(DNA_sequence[:60])
                print(translated_sequence[:60])
                DNA_sequence = DNA_sequence[60:]
                translated_sequence = translated_sequence[60:]
        else:
            print(translated_sequence)

def translate(DNA_sequence, mode):
    # Create dictionaries that translate codons into amino acids of appropriate format
    codontable = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'-', 'TAG':'-',
    'TGC':'C', 'TGT':'C', 'TGA':'-', 'TGG':'W',
    }
    # If VERBOSE, modify the start and stop codons and modify all codons to add space after

    # If DNA, modify all codons to add space before and after

    # Loop through DNA sequence codons
    if mode == 'DNA':
        out_seq = ''
    else:
        out_seq = ''
    for i in range(0, len(DNA_sequence), 3):
        if len(DNA_sequence[i:i+3]) < 3:
            break
        #Determine new amino acid with appropriate format
        if mode == 'COMPACT':
            new_codon = codontable[DNA_sequence[i:i+3]]
            out_seq += new_codon
        elif mode == 'VERBOSE':
            new_codon = codontable[DNA_sequence[i:i+3]]
            if new_codon == 'M':
                new_codon = 'Met'
            elif new_codon == '-':
                new_codon = 'Stop'
            out_seq += new_codon + ' '
        else:
            new_codon = codontable[DNA_sequence[i:i+3]]
            out_seq += ' ' + new_codon + ' '
    if out_seq:
        if out_seq[len(out_seq)-1] == ' ':
            out_seq = out_seq[:len(out_seq) - 1]
    else:
        out_seq = ''
    return out_seq

def reverse_complement(DNA_sequence):
    # Returns string containing reverse complement of DNA sequence
    complements_dic = {'A':'T','T':'A','C':'G','G':'C'}
    rev_DNA_sequence = ''
    for letter in DNA_sequence:
        rev_DNA_sequence += complements_dic[letter]
    return rev_DNA_sequence[::-1]

# Part I: Determine if the user has entered the appropriate number of argments when they called the script (one). Determine if the user entered one of three valid options for the mode. If there is an error in either of these, print out informative error messages indicating which error was made, what the three valid options are, and then quit the program.
if len(sys.argv) != 2:
    print("Invalid number of options")
    print("Mode can be one of the following options:")
    print("VERBOSE")
    print("COMPACT")
    print("DNA")
    sys.exit()

elif sys.argv[1].upper() not in ["VERBOSE", "COMPACT", "DNA"]:
    print(sys.argv[1] + " not a valid option")
    print("Mode can be one of the following options:")
    print("VERBOSE")
    print("COMPACT")
    print("DNA")
    sys.exit()

# Part II: Create loop to query user for DNA sequence
# Part III: Determine if user wants to exit program
# Part IV: Determine if user input is valid DNA sequence. If DNA sequence is not valid, print error message and allow user to enter new DNA sequence.
while True:
    user_input = input("Enter DNA sequence (or Exit to quit the program): ")
    if user_input.upper() == "EXIT":
        sys.exit()
    elif not valid_DNA_sequence(user_input):
        print("Invalid DNA sequence. Characters must be one of A, a, C, c, G, g, T, or t")
    else:
        # Part V: Print out 6 translated frames to the screen in appropriate format
        print_DNA_sequence(user_input.upper(), sys.argv[1].upper())

def reverse_complement(DNA_sequence):
    # Returns string containing reverse complement of DNA sequence
    complements_dic = {'A':'T','T':'A','C':'G','G':'C'}
    rev_DNA_sequence = ''
    for letter in DNA_sequence:
        rev_DNA_sequence = rev_DNA_sequence + complements_dic[letter]
    return rev_DNA_sequence[::-1]


print(reverse_complement('ATGACGGAG'))

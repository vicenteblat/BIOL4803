import pyfaidx, sys, pandas

def read_input_files(seq_file, ann_file):
    # This function should take in two files and open and perform error checking on them. It will return a pyfaidx, Fasta object and a pandas dataframe object
    good_flag = True # keep track of whether there was an error
    bad_nucleotides = set() # this object is similar to a list but can only hold unique values
    annotation_columns = ['GeneName', 'Chromosome', 'Strand', 'Start', 'Stop'] # Columns you expect in annotation file

    # First work on SEQUENCE FILE
    # Try to open sequence file with pyfaidx
    # #1. Is the file present?
    # #2. Can the file be opened by pyfaidx?
    # #3. Are the only nucleotides A,C,G,T?


    errors = 0
    #1. Is the file present?
    #2. Can the file be opened by pyfaidx?
    try:
        f = open(seq_file,'r')
        f.close()
        seq = pyfaidx.Fasta(seq_file)
    except FileNotFoundError:
        print("SequenceFileError: " + seq_file + " is not a valid filename")
        errors=+1
    except pyfaidx.FastaIndexingError:
        print("SequenceFileError: " + seq_file + " does not appear to be a fasta file")
        errors=+1
    else:
        #3. Are the only nucleotides A,C,G,T?
        # Return list of sub sequences in fasta file
        keys = list(seq.keys())
        # Iterates through sub sequences
        for sub_sequence_key in keys:
            # String representation of sub sequence
            sub_sequence = str(seq[sub_sequence_key])
            for char in sub_sequence:
                if char.upper() not in ['A', 'T', 'C', 'G']:
                    bad_nucleotides.add(char)

        if len(bad_nucleotides) != 0:
            print("SequenceFileError: The following bad nucleotides were found in your sequence file: " + str(list(bad_nucleotides)))
            errors=+1
        # Not sure if this is needed





    # Next work on the ANNOTATION FILE
    # Try to open sequence file with pandas
    # #1. Is the file present?
    # #2. Are there five columns?
    # #3. Are the headers what you expect?
    # #4. Is each gene listed once?
    # #5. Loop through all records. Determine if strand is +/-, start is less than stop, and length is divisible by three. Print out all error records
    # #6. Were there any errors? Exit. Otherwise return objects.


    #1. Is the file present?
    try:
        annotation = pandas.read_csv(ann_file, sep = "\t")

    except FileNotFoundError:
        print("AnnotationFileError: " + ann_file + " is not a valid filename")
        errors=+1
    else:
        #2. Are there five columns?
        #3. Are the headers what you expect?
        # n_columns = annotation.shape[1]
        headers = annotation.head(0)
        headers_list = headers.columns.tolist()
        # Checks number of columns and if headers_list contains all elements in annotation_columns
        # This next line only checks if they are in that specific order
        if headers_list != annotation_columns:
            print("AnnotationFileError: " + ann_file + " must contain five columns with the following headers: " + str(annotation_columns))
            errors=+1
        #4. Is each gene listed once?
        repeated_genes = []
        repeated_genes = annotation[annotation.duplicated()].GeneName.values.tolist()
        if repeated_genes != []:
            print("AnnotationFileError: " + ann_file + " can only contain each gene listed once. The following genes were listed more than once: " + str(repeated_genes))
            errors=+1

        #5. Loop through all records. Determine if strand is +/-, start is less than stop, and length is divisible by three.
        ### Print out all error records
        wrong_genes = annotation.loc[(annotation.Strand != '-') & (annotation.Strand != '+') | (annotation.Start > annotation.Stop) | ((annotation.Stop - annotation.Start + 1) % 3 != 0)]
        wrong_genes_T = wrong_genes.T
        for i, row in wrong_genes.iterrows():
            if (row.Strand != '-') & (row.Strand != '+'):
                print("AnnotationFileError: Strand must be + or -")
                print(wrong_genes_T[i])
                errors=+1

            if row.Start > row.Stop:
                print("AnnotationFileError: Stop must be greater than Start")
                print(wrong_genes_T[i])
                errors=+1

            if (row.Stop - row.Start + 1) % 3 != 0:
                print("AnnotationFileError: Start to stop must be divisible by three")
                print(wrong_genes_T[i])
                errors=+1


    if errors > 0:
        print("Exiting...")
        sys.exit()

    return seq, annotation

def chromosome_info(seqs, annotations):
    # This function takes in a pyfaidx.Fasta and pandas.DataFrame objects and prints out name, length, number of genes, and gc content of each chromosome
    print('\t'.join(['ContigName', 'Length', 'NumGenes', 'GC_content']))
    for contig in seqs:
        print(contig.name + ':\t' + str(len(contig)) + '\t' + str(len(annotations[annotations.Chromosome == contig.name])) + '\t' + str(contig[:].gc))

def codon_info(seqs, annotations):
    # This function prints out codon usage for all genes
    # Dictionary to hold information. One way to do this is to have a key 'All' that contains total number of amino acids, and a key/value pair 'aminoacid': dictionary, where the dictionary contains how often each codon is used. The loop below assumes use of that format (but you can do it other ways if you want) iterate through all rows of annotation file
    aa_usage = {
    'All':0,
    '-':{'All':0, 'TAA':0, 'TAG':0, 'TGA':0},
    'A':{'All':0, 'GCA':0, 'GCC':0, 'GCG':0, 'GCT':0},
    'C':{'All':0, 'TGC':0, 'TGT':0},
    'D':{'All':0, 'GAC':0, 'GAT':0},
    'E':{'All':0, 'GAA':0, 'GAG':0},
    'F':{'All':0, 'TTC':0, 'TTT':0},
    'G':{'All':0, 'GGA':0, 'GGC':0, 'GGG':0, 'GGT':0},
    'H':{'All':0, 'CAC':0, 'CAT':0},
    'I':{'All':0, 'ATA':0, 'ATC':0, 'ATT':0},
    'K':{'All':0, 'AAA':0, 'AAG':0},
    'L':{'All':0, 'CTA':0, 'CTC':0, 'CTG':0, 'CTT':0, 'TTA':0, 'TTG':0},
    'M':{'All':0, 'ATG':0},
    'N':{'All':0, 'AAC':0, 'AAT':0},
    'P':{'All':0, 'CCA':0, 'CCC':0, 'CCG':0, 'CCT':0},
    'Q':{'All':0, 'CAA':0, 'CAG':0},
    'R':{'All':0, 'AGA':0, 'AGG':0, 'CGA':0, 'CGC':0, 'CGG':0, 'CGT':0},
    'S':{'All':0, 'AGC':0, 'AGT':0, 'TCA':0, 'TCC':0, 'TCG':0, 'TCT':0},
    'T':{'All':0, 'ACA':0, 'ACC':0, 'ACG':0, 'ACT':0},
    'V':{'All':0, 'GTA':0, 'GTC':0, 'GTG':0, 'GTT':0},
    'W':{'All':0, 'TGG':0},
    'Y':{'All':0, 'TAC':0, 'TAT':0}}

    # iterate through all rows of annotation file
    for i, row in annotations.iterrows():
        # print(row)
        # print(row[0])
        gene_seq, protein_seq = ret_nucleotide_protein_sequence(seqs, annotations, row[0])
        # Populate dictionary info
        count = 0
        for aminoacid in protein_seq:
            aa_usage['All'] += 1
            aa_usage[aminoacid]['All'] += 1
            aa_usage[aminoacid][gene_seq[count:count+3]] += 1
            count += 3

    for aa in sorted(aa_usage.keys()):
        if aa != 'All':
            print(aa + '\t' + "{0:.2f}".format(100*aa_usage[aa]['All']/aa_usage['All']) + '% - ', end = '')
            for codon in sorted(aa_usage[aa].keys()):
                if codon != 'All':
                    print(codon + ': ' + "{0:.2f}".format(100*aa_usage[aa][codon]/aa_usage[aa]['All']) + '%. ', end = '')
            print()

def ret_nucleotide_protein_sequence(seqs, annotations, gene):
    # This function takes in a pyfaidx.Fasta and pandas.DataFrame objects and a string containing a gene name and returns the nucleotide sequence and amino acid sequence of the gene. It should raise a KeyError if the gene doesn't exist in the annotation file
    start = int(annotations.loc[annotations['GeneName'] == gene, 'Start']) - 1
    # print(start)
    stop = int(annotations.loc[annotations['GeneName'] == gene, 'Stop'])
    # print(stop)
    chromosome = annotations.loc[annotations['GeneName'] == gene, 'Chromosome']._get_values(0)
    # print(chromosome)

    strand = annotations.loc[annotations['GeneName'] == gene, 'Strand']._get_values(0)

    sequence = str(seqs[chromosome][start:stop])

    if strand == '-':
        sequence = reverse_complement(sequence)


    return sequence, translate(sequence)

def translate(sequence):
    # This function translates DNA sequence into amino acid sequence. Can be used on individual codons or cDNA sequences.
    out_seq = ''
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

    for i in range(0, len(sequence), 3):
        out_seq += codontable[sequence[i:i+3]]

    return out_seq


def reverse_complement(DNA_sequence):
    # Returns string containing reverse complement of DNA sequence
    complements_dic = {'A':'T','T':'A','C':'G','G':'C'}
    rev_DNA_sequence = ''
    for letter in DNA_sequence:
        rev_DNA_sequence += complements_dic[letter]
    return rev_DNA_sequence[::-1]

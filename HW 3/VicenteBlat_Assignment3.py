import argparse
import VicenteBlat_UsefulFunctions as UF


# Part I: Use argparse to define the required and optional arguments for this script. For the gene flag, you should use nargs to accept multiple arguments
parser = argparse.ArgumentParser(description = "This program is used to analyze bacterial genomes")
parser.add_argument("SequenceFile", help = "Fasta file containing the DNA sequence of the bacterial chromosomes")
parser.add_argument("AnnotationFile", help = "Tab-delimited file containing the location of all the genes. \
                                               File format should be: <GeneName> <Chromosome> <Strand> <Start> <Stop>")
parser.add_argument("-c", "--codons", help = "Run analysis of amino acid and codon usage", action = "store_true")
parser.add_argument("-g", "--GENES", help = "Return protein sequence of a specific gene or set of genes", nargs = '+')
args = parser.parse_args()

# Read in sequence and error file and parse it using pyfaidx and pandas
sequence_obj, annotation_obj = UF.read_input_files(args.SequenceFile, args.AnnotationFile)

# Create logic to determine if optional flags are present. If not run normal analysis
if args.codons is False and args.GENES is None:
    UF.chromosome_info(sequence_obj, annotation_obj)
else:
    # Create logic to determine if codon flag is present. If it is run codon analysis
    if args.codons:
        UF.codon_info(sequence_obj, annotation_obj)
    # Create logic to determine if it is loop through genes and print out sequence
    if args.GENES:
        for gene in args.GENES:
            try:
                coding_seq, protein_seq = UF.ret_nucleotide_protein_sequence(sequence_obj,annotation_obj,gene)
                print('>' + gene)
                print(protein_seq)
            except TypeError:
                print('Unable to find ' + gene + ' in the annotation file')

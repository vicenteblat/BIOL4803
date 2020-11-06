from Bio.PDB import *
from Bio import pairwise2
from Bio.SeqRecord import SeqRecord
from Bio.SubsMat import MatrixInfo as matlist
from Bio.Align import MultipleSeqAlignment, Seq

pdbp = PDBParser(QUIET = True)
ref_structure = pdbp.get_structure('Ref', 'B2R_Inactive.pdb')
mov_structure = pdbp.get_structure('Mov', 'B2R_Active.pdb')


#2: Get the sequence of two pdb files and save the sequences as SeqRecord objects
ppb = PPBuilder()
polypeptides_ref = ppb.build_peptides(ref_structure)
polypeptides_mov = ppb.build_peptides(mov_structure)

for pp_ref in polypeptides_ref:
    try:
        sequence_ref += SeqRecord(pp_ref.get_sequence())
    except:
        sequence_ref = SeqRecord(pp_ref.get_sequence())

for pp_mov in polypeptides_mov:
    try:
        sequence_mov += SeqRecord(pp_mov.get_sequence())
    except:
        sequence_mov = SeqRecord(pp_mov.get_sequence())


#3: Create a pairwise alignment of the two sequences using BioPython. Use the blosum62 matrix with a gap open penalty of -10 and a gap extend penalty of -0.5. Create a MultipleSeqAlignment object and save it as an object attribute
matrix = matlist.blosum62
gap_open = -10
gap_extend = -0.5

alns = pairwise2.align.globalds(sequence_ref.seq, sequence_mov.seq, matrix, gap_open, gap_extend)

top_aln = alns[0]
alignment = MultipleSeqAlignment([SeqRecord(Seq(top_aln[0])),SeqRecord(Seq(top_aln[1]))])

structure_alignment = StructureAlignment(alignment, ref_structure, mov_structure)

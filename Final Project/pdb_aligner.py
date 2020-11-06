import Bio.PDB
from Bio.SeqIO import SeqRecord
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
from Bio.Align import MultipleSeqAlignment, Seq
import pdb

# Read in structre
pdbp = Bio.PDB.PDBParser(QUIET = True)
structure = pdbp.get_structure('GluRec', '3RIF.pdb')

# Get protein sequence
ppb = Bio.PDB.PPBuilder()
polypeptides = ppb.build_peptides(structure)
seq1 = polypeptides[0].get_sequence()
seq2 = polypeptides[1].get_sequence()


matrix = matlist.blosum62
gap_open = -10
gap_extend = -0.5

alns = pairwise2.align.globalds(seq1,seq2,matrix, gap_open, gap_extend)
top_aln = alns[0]

alignment = MultipleSeqAlignment([SeqRecord(Seq(top_aln[0])),SeqRecord(Seq(top_aln[1]))])
structure_alignment = Bio.PDB.StructureAlignment(alignment,structure[0]['A'], structure[0]['B'])

sup = Bio.PDB.Superimposer()
ref_atoms = []
mov_atoms = []
for duo in structure_alignment.duos:
	res1 = duo[0]
	res2 = duo[1]
	if res1 and res2:
		ref_atoms.append(res1['CA'])
		mov_atoms.append(res2['CA'])

sup.set_atoms(ref_atoms, mov_atoms)
sup.apply(structure[0]['B'].get_atoms())
io = Bio.PDB.PDBIO()
io.set_structure(structure)
io.save('Temp.pdb')

from Bio.PDB import *
from Bio import pairwise2
from Bio.SeqRecord import SeqRecord
from Bio.SubsMat import MatrixInfo as matlist
from Bio.Align import MultipleSeqAlignment, Seq

class StructureComparer():
    def __init__(self, ref_pdb_file, mov_pdb_file):

        #1: Read in these pdbs using Bio.PDB and save them as object attributes
        pdbp = PDBParser(QUIET = True)
        self.ref_structure = pdbp.get_structure('Ref', ref_pdb_file)
        self.mov_structure = pdbp.get_structure('Mov', mov_pdb_file)


        #2: Get the sequence of two pdb files and save the sequences as SeqRecord objects
        ppb = PPBuilder()
        polypeptides_ref = ppb.build_peptides(self.ref_structure)
        polypeptides_mov = ppb.build_peptides(self.mov_structure)

        for pp_ref in polypeptides_ref:b
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
        self.alignment = MultipleSeqAlignment([SeqRecord(Seq(top_aln[0])),SeqRecord(Seq(top_aln[1]))])


        #4: Create a StructureAlignment object to pair the two pdb objects
        self.structure_alignment = StructureAlignment(self.alignment, self.ref_structure, self.mov_structure)

        #5: Create a list of atoms to use for the structural alignment. Use the duos that fall within one of the 7 transmembrane regions (based on the reference file). Use the 'N', 'CA', 'O', and 'C' atoms in the backbone.
        self.ref_residues = []
        ref_atoms = []
        self.mov_residues = []
        mov_atoms = []
        for duo in self.structure_alignment.duos:
            res1 = duo[0]
            res2 = duo[1]
            if res1 and res2:
                self.ref_residues.append(res1)
                self.mov_residues.append(res2)
                # ranges = [(31,56),(70,92),(104,129),(150,171),(198,220),(275,298),(306,327)]
                # if any(lower <= res1.get_id()[1] <= upper for (lower,upper) in ranges):
                if self.ret_region(res1.get_id()[1]):
                    ref_atoms.append(res1['N'])
                    mov_atoms.append(res2['N'])
                    ref_atoms.append(res1['CA'])
                    mov_atoms.append(res2['CA'])
                    ref_atoms.append(res1['O'])
                    mov_atoms.append(res2['O'])
                    ref_atoms.append(res1['C'])
                    mov_atoms.append(res2['C'])

        #6: Calculate the rotation and translation matrices using the Superimposer object using the ref_atoms and mov_atoms. Apply these to the second pdb file. You probably want to save the aligned pdb file to verify it is working.
        sup = Superimposer()
        sup.set_atoms(ref_atoms, mov_atoms)
        sup.apply(self.mov_structure.get_atoms())
        io = PDBIO()
        io.set_structure(self.mov_structure)
        io.save('Temp.pdb')

        for residue in self.ref_structure.get_residues():
            if residue.resname == 'LIG':
                self.ligand = residue

    def avg_residue_distance(self, res1, res2):
        # This method takes in two residues and returns the average distance between all the atoms in the two residues
        sum = 0
        count = 0
        for atom1 in res1.get_atoms():
            for atom2 in res2.get_atoms():
                sum += atom2 - atom1
                count += 1

        return sum/count



    def min_residue_distance(self, res1, res2):
        # This method takes in two residues and returns the minimum distance between the two residues
        min = list(res2.get_atoms())[0] - list(res1.get_atoms())[0]
        for atom1 in res1.get_atoms():
            for atom2 in res2.get_atoms():
                if atom2 - atom1 < min:
                    min = atom2 - atom1

        return min


    def ret_ligand(self, pdb, ligand_name):
        # This method returns the ligand residue for a pdb file (named ligand name)
        pdbp = PDBParser(QUIET = True)
        structure = pdbp.get_structure('pdb', pdb)
        for residue in structure.get_residues():
            if residue.resname == ligand_name:
                return residue


    def ret_region(self, loc):
        # This function returns the transmembrane domain of a residue found in the b2 inactive receptor. You do not need to modify it
        if loc >= 31 and loc <= 56:
            return 'TM1'
        if loc >= 70 and loc <= 92:
            return 'TM2'
        if loc >= 104 and loc <= 129:
            return 'TM3'
        if loc >= 150 and loc <= 171:
            return 'TM4'
        if loc >= 198 and loc <= 220:
            return 'TM5'
        if loc >= 275 and loc <= 298:
            return 'TM6'
        if loc >= 306 and loc <= 327:
            return 'TM7'
        return None

    def ret_residue_distances(self):
        # This function returns a dictionary of lists containing the reference location (key'd by the transmembrane domain) and a dictionary of lists containing the distances between two superimposed residues (key'd by the transmembrane domain)
        ref_loc = {'TM1':[], 'TM2':[], 'TM3':[], 'TM4':[], 'TM5':[], 'TM6':[], 'TM7':[]}
        ref_dis = {'TM1':[], 'TM2':[], 'TM3':[], 'TM4':[], 'TM5':[], 'TM6':[], 'TM7':[]}

        for duo in self.structure_alignment.duos:
            res1 = duo[0]
            res2 = duo[1]
            if res1 and res2:
                if self.ret_region(res1.get_id()[1]):
                    ref_loc[self.ret_region(res1.get_id()[1])].append(res1.get_id()[1])
                    ref_dis[self.ret_region(res1.get_id()[1])].append(self.avg_residue_distance(res1, res2))

        return ref_loc, ref_dis

    def ret_ligand_distances(self):
        # This function prints returns any residues closer than 4 angstroms to the ligand.
        residues_list = []
        ligand = self.ligand

        for residue in self.ref_structure.get_residues():
            # if residue.resname == 'LIG':
            #     import pdb; pdb.set_trace()
            #     print(residue)
            if self.avg_residue_distance(ligand, residue) <= 4 and residue.resname != 'LIG':
                residues_list.append(residue)

        return residues_list

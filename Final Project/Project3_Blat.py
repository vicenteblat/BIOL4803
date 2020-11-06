# Change the name of the module as appropriate to your last name
from HelperClass import StructureComparer

import argparse
import matplotlib.pyplot as plt

# Part A. Set up the options for the program (Don't change anything)
parser = argparse.ArgumentParser()
parser.add_argument('Inactivated_pdb', type = str, help = 'Name of a pdb file containing reference pdb')
parser.add_argument('Activated_pdb', type = str, help = 'Name of a pdb file containing activated pdb')
parser.add_argument('Paralogous_pdb', type = str, help = 'Name of a pdb file containing activated pdb')

args = parser.parse_args()

# Part B. Create StructureComparer objects using the files taken in by the argparse module (Don't change anything)
inact_act = StructureComparer(args.Inactivated_pdb, args.Activated_pdb)
inact_par = StructureComparer(args.Inactivated_pdb, args.Paralogous_pdb)

# Part C. Run the ret_residue_distances() method to return the distances between the aligned pdb files (Don't change anything)
ref_loc1, act_dis = inact_act.ret_residue_distances()
ref_loc2, par_dis = inact_par.ret_residue_distances()

inact_act_residues = inact_act.ret_ligand_distances()
print(inact_act_residues)


# Part D. Plot distances between the aligned pdb files (Don't change anything)

for key in ref_loc1:
    a, = plt.plot(ref_loc1[key], act_dis[key], label = 'Active/Inactive', color = 'black')
for key in ref_loc2:
    b, = plt.plot(ref_loc2[key], par_dis[key], label = 'Paralogous', color = 'red')

plt.legend(handles = [a,b])
plt.xlabel('Position (amino acids)')
plt.ylabel('Distance from reference (Angstroms)')
plt.ylim([0,20])
plt.text(35, 6, 'tm1')
plt.text(75, 6, 'tm2')
plt.text(110, 6, 'tm3')
plt.text(155, 6, 'tm4')
plt.text(205, 6, 'tm5')
plt.text(280, 6, 'tm6')
plt.text(310, 6, 'tm7')

plt.show()

# Part E. Create chimera command file from ref_loc1 and act_dis to color all of the residues that move more than 4 angstroms red
# (Change this section)
with open('Moving_residues.com', 'w') as f:
    for key in act_dis:
        distances = act_dis[key]
        x = 0
        for dis in distances:
            if dis > 4:
                location = ref_loc1[key][x]
                for residue in inact_act.ref_structure.get_residues():
                    if residue.id[1] == location:
                        f.write('color red :' + residue.resname)
            x += 1
f.close()


# Part F. Print out position and amino acid type for all residues within 4 angstroms of the ligand
# Also create a chimera command file from ligand_residues to color all of the resdiues within 4 angstroms of the ligand red
# Modify this section
inact_act_residues = inact_act.ret_ligand_distances()
for residue in inact_act_residues:
    print('Residue: ' + residue.resname + ' | Location: ' + str(residue.id[1]) + ' | Distance: ' + str(inact_act.min_residue_distance(residue,inact_act.ligand)))


with open('Ligand_residues.com', 'w') as f:
    for residue in inact_act_residues:
        f.write('color green :' + residue.resname)
f.close()

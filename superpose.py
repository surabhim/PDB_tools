import Bio.PDB
from math import acos, atan2, cos, pi, sin, asin
import sys

# ARGV[1] : pdb of the  reference protein
# ARGV[2] : pdb of the sample protein
# This script superposes the sample protein onto the reference protein using the pdb coordinates.

pdb_parser = Bio.PDB.PDBParser(QUIET = True)

ref_structure = pdb_parser.get_structure("reference",sys.argv[1])

sample_structure = pdb_parser.get_structure("samlpe",sys.argv[2])

ref_model    = ref_structure[0]
sample_model = sample_structure[0]

# Make a list of the atoms (in the structures) you wish to align.
# In this case we use CA atoms whose index is in the specified range
ref_atoms = []
sample_atoms = []

# Iterate of all chains in the model in order to find all residues
for ref_chain in ref_model:
  # Iterate of all residues in each model in order to find proper atoms
  for ref_res in ref_chain:
        ref_atoms.append(ref_res['CA'])

# Do the same for the sample structure
for sample_chain in sample_model:
  for sample_res in sample_chain:
    sample_atoms.append(sample_res['CA'])
    
# Now we initiate the superimposer:
super_imposer = Bio.PDB.Superimposer()
super_imposer.set_atoms(ref_atoms, sample_atoms)
super_imposer.apply(sample_model.get_atoms())

# Print RMSD:
#print super_imposer.rotran

translation=super_imposer.rotran[1]
#print translation
rot_matrix=super_imposer.rotran[0]

x1=atan2(rot_matrix[1,2], rot_matrix[2,2])
y1=asin(rot_matrix[0,2])
z1=atan2(rot_matrix[0,1],rot_matrix[0,0])

#print x1,y1,z1

x2=translation[0]
y2=translation[1]
z2=translation[2]
#print x2,y2,z2

print sys.argv[3],round(x1,5),round(y1,5),round(z1,5),round(x2,5),round(y2,5),round(z2,5)





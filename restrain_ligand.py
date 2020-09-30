#!/usr/bin/env python

"""
Calculates distances between neighboring residues of a ligand molecule and produces a set of
unambiguous distance restraints for HADDOCK to keep it in place during semi-flexible refinement.

Produces, at most, one restraint per ligand atom.

Creates a file containing the restraints, named after the original input PDB file and the ligand
residue name, with '.tbl' extension.
"""

from __future__ import print_function
import argparse
import os
import sys

try:
    import numpy as np
    from Bio.PDB import PDBParser
    from Bio.PDB import NeighborSearch
    from Bio.PDB.Polypeptide import is_aa
except ImportError as e:
    print('[!] Could not import module \'biopython\': {0}'.format(e), file=sys.stderr)
    sys.exit(1)

cline_parser = argparse.ArgumentParser(description=__doc__)
cline_parser.add_argument('pdbf', help='PDB file')
cline_parser.add_argument('-l', '--ligand', help='Ligand residue name', required=True)
cline_parser.add_argument('-p', '--pml', action='store_true', help='Write Pymol file with restraints')
args = cline_parser.parse_args()

if not os.path.isfile(args.pdbf):
    print('[!!] File not found: {0}'.format(args.pdbf))
    sys.exit(1)
else:
    root_name, _ = os.path.splitext(os.path.basename(args.pdbf))

# Read in structure
pdb_parser = PDBParser(QUIET=1)
structure = pdb_parser.get_structure(root_name, args.pdbf)

# Remove hydrogens
atom_lst = list(structure.get_atoms())
for atom in atom_lst:
    if atom.element == 'H':
        res = atom.parent
        res.detach_child(atom.name)

ligand = None
for residue in structure.get_residues():
    if residue.resname.strip() == args.ligand.strip():
        ligand = residue
        break

if not ligand:
    print('[!!] Ligand residue \'{0}\' not found in structure'.format(args.ligand), file=sys.stderr)
    sys.exit(1)

# Calculate center of mass of the ligand
ligand_com = map(lambda x: sum(x)/len(x), zip(*[at.coord for at in ligand]))
ligand_com = np.asarray(ligand_com, dtype=np.float32)

# Calculate neighbors considering only aminoacid/nucleotide atoms (excl. waters, other ligands, etc)
sel_atoms = [at for at in structure.get_atoms() if at.parent.id[0] == ' ']
ns = NeighborSearch(sel_atoms)
neighbors = ns.search(ligand_com, 10.0, level='R') # 10A radius, return residues

# Calculate residue closer to each ligand atom and the respective distance
ligand_atoms = ligand.child_list
min_dist_list, _seen = [], set()

for l_at in ligand_atoms:
    distances = []
    for residue in neighbors:
        for r_at in residue:
            distances.append((r_at, l_at, r_at - l_at))

    distances.sort(key=lambda x: x[-1])
    min_dist = distances[0]
    # One restraint per residue to keep the number of restraints small
    if min_dist[0].parent not in _seen:
        min_dist_list.append(min_dist)
        _seen.add(min_dist[0].parent)

# Output
# TBL: use chain ID as segid ID
with open('{}_{}.tbl'.format(root_name, args.ligand), 'w') as tbl_handle:
    _str = 'assi (segi {:4s} and resi {:4d} and name {:6s})\n' + \
           '     (segi {:4s} and resi {:4d} and name {:6s}) {:6.3f} 1.0 1.0'

    print('! Restraints to fix \'{}\' in its initial position'.format(args.ligand), file=tbl_handle)
    for dist in min_dist_list:
        r_at, l_at, d = dist
        print(_str.format(r_at.parent.parent.id, r_at.parent.id[1], '"' + r_at.name + '"',
                          l_at.parent.parent.id, l_at.parent.id[1], '"' + l_at.name + '"',
                          d), file=tbl_handle)

# PML (if requested)
if args.pml:
    with open('{}_{}.pml'.format(root_name, args.ligand), 'w') as pml_handle:
        _str = 'distance rst_{}, c. {} and i. {} and n. "{}", c. {} and i. {} and n. "{}"'
        for irest, dist in enumerate(min_dist_list, start=1):
            r_at, l_at, _ = dist
            print(_str.format(irest,
                              r_at.parent.parent.id, r_at.parent.id[1], r_at.name,
                              l_at.parent.parent.id, l_at.parent.id[1], l_at.name), file=pml_handle)


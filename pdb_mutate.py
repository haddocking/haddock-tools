#!/usr/bin/python
"""
(dummy-)Mutates single residue on a PDB-formatted structure.
HADDOCK will then reconstruct the residue according to its topology.

Usage1: python pdb_mutate.py <pdb file> <mutation chain> <mutation residue number> <wildtype residue name> <mutant residue name>
Example: python pdb_mutate.py 1A22.pdb A 7 SER ALA

Usage2: python pdb_mutate.py <mutation list file>
        The format of mutation list:
        pdbfile chain resi resn_wt resn_mut

Example: python pdb_mutate.py list_mutations
        for mutating residue 7 Serine of chain A to Alanine:
        1A22.pdb A 7 SER ALA

Author: {0}
Email: {1}
"""
from __future__ import print_function

import os
import sys

__author__ = "Cunliang Geng; Joao Rodrigues"
__email__ = "gengcunliang@gmail.com; j.p.g.l.m.rodrigues@gmail.com"

USAGE = __doc__.format(__author__, __email__)


def check_input(args):
    """Checks whether to read from stdin/file and validates user input/options."""

    if not len(args):
        # Read from pipe
        if not sys.stdin.isatty():
            mutfh = sys.stdin
        else:
            sys.stderr.write(USAGE)
            sys.exit(1)
    elif len(args) == 1 or len(args) == 5:
        if not os.path.isfile(args[0]):
            sys.stderr.write('File not found: ' + args[0] + '\n')
            sys.stderr.write(USAGE)
            sys.exit(1)
        mutfh = open(args[0], 'r')
    else:
        sys.stderr.write(USAGE)
        sys.exit(1)
    return mutfh


def mutate(structure_fhandle, chain, resi, resn_wt, resn_mut):

    mutated_structure = []
    for line in structure_fhandle:
        if line[0:4] == 'ATOM' or line[0:6] == 'HETATM':
            s_chain = line[21].strip()
            s_resi = line[22:26].strip()
            s_resn = line[17:20].strip()
            if s_chain == chain and s_resi == resi:
                # check whether wildtype residue name is same as given resn.
                if s_resn == resn_wt:
                    line = line[0:17]+resn_mut+line[20:]
                    mutated_structure.append(line)
                else:
                    sys.stderr.write('Error: Wildtype residue of chain {0} resi {1} is not {2} but {3}.\n'.
                                     format(chain, resi, resn_wt, s_resn))
                    return None
            else:
                mutated_structure.append(line)
        else:
            mutated_structure.append(line)
    return mutated_structure


def _print_mutant(mut_fhandle, chain, resi, resn_wt, resn_mut):
    mutant = mutate(mut_fhandle, chain, resi, resn_wt, resn_mut)
    mutpdb = ''.join(mutant).rstrip("\n")
    print(mutpdb)


def _print_mutants(mut_fhandle):

    mut_list = [l.split() for l in mut_fhandle if l.strip()]

    print("Generated mutant files:")
    for i in mut_list:
        if len(i) == 5:
            pdb_path, chain, resi, resn_wt, resn_mut = i
        else:
            sys.stderr.write('Error: Unrecognized mutation format in line "{0}"\n'.format(" ".join(i)))
            continue

        if not os.path.exists(pdb_path):
            sys.stderr.write('Error: PDB file {0} not found\n'.format(pdb_path))
            continue

        with open(pdb_path) as pdb_fhandle:
            structure = [l for l in pdb_fhandle]

        mutant = mutate(structure, chain, resi, resn_wt, resn_mut)
        if mutant:
            m_file = open('{0}_{1}_{2}{3}{4}.pdb'.format(os.path.basename(pdb_path).split('.')[0], chain, resn_wt, resi,
                                                         resn_mut), 'w')
            print("{0}_{1}_{2}{3}{4}.pdb".format(os.path.basename(pdb_path).split('.')[0], chain, resn_wt, resi,
                                                 resn_mut))
            m_file.write(''.join(mutant))
            m_file.close()
        else:
            continue


if __name__ == "__main__":
    mut_fhandle = check_input(sys.argv[1:])

    if len(sys.argv[1:]) == 1:
        _print_mutants(mut_fhandle)
    else:
        _print_mutant(mut_fhandle, *sys.argv[2:])

    mut_fhandle.close()

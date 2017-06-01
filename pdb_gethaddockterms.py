#!/usr/bin/python
"""
Extract HADDOCK terms (Eair, Evdw, Eelec, Edesolv and BSA) from a HADDOCK-format PDB file.

Usage: python pdb_gethaddockterms.py <HADDOCK PDB file>
Example: python pdb_gethaddockterms.py cluster1_1.pdb

Author: {0}
Email: {1}
"""
from __future__ import print_function
import os
import sys
import re

__author__ = "Cunliang Geng"
__email__ = "gengcunliang@gmail.com"

USAGE = __doc__.format(__author__, __email__)

def check_input(args):
    """Checks whether to read from stdin/file and validates user input/options."""

    if not len(args):
        # Read from pipe
        if not sys.stdin.isatty():
            fhandle = sys.stdin
        else:
            sys.stderr.write(USAGE)
            sys.exit(1)
    elif len(args) == 1:
        if not os.path.isfile(args[0]):
            sys.stderr.write('File not found: ' + args[0] + '\n')
            sys.stderr.write(USAGE)
            sys.exit(1)
        fhandle = open(args[0], 'r')
    else:
        sys.stderr.write(USAGE)
        sys.exit(1)
    return fhandle


def ExtractHaddockTerms(pdbfhandle):
    terms = {}
    for line in pdbfhandle:
        if line[0:6] == 'REMARK':
            if re.match(r'REMARK energies', line):
                ene = line.split(',')
                terms['Evdw'] = ene[5].strip()
                terms['Eelec'] = ene[6].strip()
                terms['Eair'] = ene[7].strip()
            elif re.match(r'REMARK Desolvation energy', line):
                terms['Edesolv'] = line.split()[3].strip()
            elif re.match(r'REMARK buried surface area', line):
                terms['BSA'] = line.split()[4].strip()
            else:
                continue

    return terms


def WriteHaddockTerms(pdbfhandle):
    terms = ExtractHaddockTerms(pdbfhandle)
    if terms:
        header = ['Eair', 'Evdw', 'Eelec', 'Edesolv', 'BSA']
        haddockterms = [terms[i] for i in header]
        print("\t".join(header))
        print("\t".join(haddockterms))
    else:
        sys.stderr.write("Error: Input is not a HADDOCK-format PDB file\n")
        sys.stderr.write(USAGE)
        sys.exit(1)


if __name__ == "__main__":

    pdbfhandle =  check_input(sys.argv[1:])
    WriteHaddockTerms(pdbfhandle)
    pdbfhandle.close()

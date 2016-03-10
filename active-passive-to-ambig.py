#!/usr/bin/env python

"""
Python script to convert a list of active and passive residues into 
ambiguous interaction restraints for HADDOCK
"""

import os, sys, time
import subprocess
import tempfile

def active_passive_to_ambig(active1, passive1, active2, passive2, 
        segid1='A', segid2='B', out='ambig.tbl'):
    """Convert active and passive residues to Ambiguous Interaction Restraints

    Parameters
    ----------
    active1 : list
        List of active residue numbers of the first segid

    passive1 : list
        List of passive residue numbers of the first segid

    active2 : list
        List of active residue numbers of the second segid
    
    active2 : list
        List of passive residue numbers of the second segid

    segid1 : string
        Segid to use for the first model

    segid2 : string
        Segid to use for the second model

    out : string
        Name of the output AIR file.
    """

    all1 = active1 + passive1
    all2 = active2 + passive2

    with open(out, 'w') as fh:
        for resi1 in active1:
            fh.write('assign (resi {:d} and segid {:s})\n'.format(resi1, segid1))
            fh.write('(\n')
            lines = []
            for resi2 in all2:
                lines.append('       (resi {:d} and segid {:s})\n'.format(resi2, segid2))
            lines = '        or\n'.join(lines)
            for line in lines:
                fh.write(line)

            fh.write(') 2.0 2.0 0.0\n\n')
                
        for resi2 in active2:
            fh.write('assign (resi {:d} and segid {:s})\n'.format(resi2, segid2))
            fh.write('(\n')
            lines = []
            for resi1 in all1:
                lines.append('       (resi {:d} and segid {:s})\n'.format(resi1, segid1))
            lines = '        or\n'.join(lines)
            for line in lines:
                fh.write(line)

            fh.write(') 2.0 2.0 0.0\n\n')

def main():
    import sys
    if len(sys.argv) != 3:
        print '\nUsage:\n     python active-passive_to_ambig.py <active-passive-file1> <active-passive-file2>\n\n' +\
              'where <active-passive-file> is a file consisting of two space-delimited lines with\n' +\
              'the first line active residues numbers and the second line passive residue numbers\n'
        sys.exit()

    active1, passive1 = [[int(x) for x in line.split()] for line in open(sys.argv[1])]
    active2, passive2 = [[int(x) for x in line.split()] for line in open(sys.argv[2])]
    active_passive_to_ambig(active1, passive1, active2, passive2)


if __name__ == '__main__':
    main()

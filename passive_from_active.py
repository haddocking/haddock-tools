#!/usr/bin/env python

from __future__ import print_function

"""
Output a list of active residues surrounding active residues given as input

usage: passive_from_active.py <pdb> <active-residues> [-c CHAIN_ID] [-s SURFACE_RES_LIST]
"""

import argparse
import os
import sys

import contextlib

try:
    from Bio.PDB import *
except ImportError as e:
    print("Could not import BioPython\n{0}".format(e))
    sys.exit(1)


# Scaling factors for relative ASA
# Calculated using extended ALA-X-ALA peptides
# Taken from NACCESS
rel_asa = {
    'total':
        {
            'ALA': 107.95,
            'CYS': 134.28,
            'ASP': 140.39,
            'GLU': 172.25,
            'PHE': 199.48,
            'GLY': 80.10,
            'HIS': 182.88,
            'ILE': 175.12,
            'LYS': 200.81,
            'LEU': 178.63,
            'MET': 194.15,
            'ASN': 143.94,
            'PRO': 136.13,
            'GLN': 178.50,
            'ARG': 238.76,
            'SER': 116.50,
            'THR': 139.27,
            'VAL': 151.44,
            'TRP': 249.36,
            'TYR': 212.76,
        },
    'bb':
        {
            'ALA': 38.54,
            'CYS': 37.53,
            'ASP': 37.70,
            'GLU': 37.51,
            'PHE': 35.37,
            'GLY': 47.77,
            'HIS': 35.80,
            'ILE': 37.16,
            'LYS': 37.51,
            'LEU': 37.51,
            'MET': 37.51,
            'ASN': 37.70,
            'PRO': 16.23,
            'GLN': 37.51,
            'ARG': 37.51,
            'SER': 38.40,
            'THR': 37.57,
            'VAL': 37.16,
            'TRP': 38.10,
            'TYR': 35.38,
        },
    'sc':
        {
            'ALA': 69.41,
            'CYS': 96.75,
            'ASP': 102.69,
            'GLU': 134.74,
            'PHE': 164.11,
            'GLY': 32.33,
            'HIS': 147.08,
            'ILE': 137.96,
            'LYS': 163.30,
            'LEU': 141.12,
            'MET': 156.64,
            'ASN': 106.24,
            'PRO': 119.90,
            'GLN': 140.99,
            'ARG': 201.25,
            'SER': 78.11,
            'THR': 101.70,
            'VAL': 114.28,
            'TRP': 211.26,
            'TYR': 177.38,
        }
}


@contextlib.contextmanager
def stdchannel_redirected(stdchannel, dest_filename):
    """
    A context manager to temporarily redirect stdout or stderr
    https://stackoverflow.com/questions/977840/redirecting-fortran-called-via-f2py-output-in-python/978264#978264
    e.g.:
    with stdchannel_redirected(sys.stderr, os.devnull):
        if compiler.has_function('clock_gettime', libraries=['rt']):
            libraries.append('rt')
    """
    oldstdchannel, dest_file = None, None
    try:
        oldstdchannel = os.dup(stdchannel.fileno())
        dest_file = open(dest_filename, 'w')
        os.dup2(dest_file.fileno(), stdchannel.fileno())

        yield
    finally:
        if oldstdchannel is not None:
            os.dup2(oldstdchannel, stdchannel.fileno())
        if dest_file is not None:
            dest_file.close()


def get_surface_resids(structure, cutoff=15, config_path=os.environ.get('FREESASA_CONFIG')):
    """
    Calls freesasa using its Python API and returns
    per-residue accessibilities.
    """
    try:
        from freesasa import Classifier, structureFromBioPDB, calc
    except ImportError as err:
        print('[!] The binding affinity prediction tools require the \'freesasa\' Python API', file=sys.stderr)
        raise ImportError(err)
    import pkg_resources

    asa_data, rsa_data, rel_main_chain, rel_side_chain = {}, {}, {}, {}
    _rsa = rel_asa['total']
    _rsa_bb = rel_asa['bb']
    _rsa_sc = rel_asa['sc']

    classifier = Classifier(config_path)
    pkg_resources.cleanup_resources()

    with stdchannel_redirected(sys.stderr, os.devnull):
        struct = structureFromBioPDB(structure, classifier, )
        result = calc(struct)

    # iterate over all atoms to get SASA and residue name
    for idx in range(struct.nAtoms()):
        atname = struct.atomName(idx).strip()
        resname = struct.residueName(idx)
        resid = int(struct.residueNumber(idx))
        chain = struct.chainLabel(idx)
        at_uid = (chain, resname, resid, atname)
        res_uid = (chain, resname, resid)

        asa = result.atomArea(idx)
        asa_data[at_uid] = asa
        # add asa to residue
        rsa_data[res_uid] = rsa_data.get(res_uid, 0) + asa

        if atname in ('C', 'N', 'O'):
            rel_main_chain[res_uid] = rel_main_chain.get(res_uid, 0) + asa
        else:
            rel_side_chain[res_uid] = rel_side_chain.get(res_uid, 0) + asa

    # convert total asa ro relative asa
    rsa_data.update((res_uid, asa / _rsa[res_uid[1]]) for res_uid, asa in rsa_data.items())
    rel_main_chain.update((res_uid, asa / _rsa_bb[res_uid[1]] * 100) for res_uid, asa in rel_main_chain.items())
    rel_side_chain.update((res_uid, asa / _rsa_sc[res_uid[1]] * 100) for res_uid, asa in rel_side_chain.items())

    # We format to fit the pipeline
    resid_access = {}
    for res_uid, access in rel_main_chain.items():
        resid_access[res_uid[2]] = {'side_chain_rel': rel_side_chain.get(res_uid), 'main_chain_rel': access}
    surface_resids = [r for r, v in resid_access.items() if v['side_chain_rel'] >= cutoff or
                      v['main_chain_rel'] >= cutoff]
    return surface_resids


if __name__ == "__main__":

    args_parser = argparse.ArgumentParser(description=__doc__)
    args_parser.add_argument('pdb_file', type=str, help='PDB file')
    args_parser.add_argument('active_list', type=str, help='List of active residues IDs (int) separated by commas')
    args_parser.add_argument('-c', '--chain-id', type=str, help='Chain id to be used in the PDB file (default: All)')
    args_parser.add_argument('-s', '--surface-list', type=str,
                             help='List of surface residues IDs (int) separated by commas')
    args = args_parser.parse_args()

    # Parse the PDB file
    if os.path.isfile(args.pdb_file):
        try:
            p = PDBParser(QUIET=True)
            s = p.get_structure('pdb', args.pdb_file)
        except Exception as e:
            print('Error while parsing the PDB file: {0}'.format(e))
            sys.exit(1)
    else:
        print('File not found: {0}'.format(args.pdb_file))
        sys.exit(1)

    try:
        if args.chain_id:
            atom_list = [a for a in s[0][args.chain_id].get_atoms()]
        else:
            atom_list = [a for a in s[0].get_atoms()]
    except KeyError as e:
        print('Chain {0} does not exist in the PDB file {1}, please enter a proper chain id'.
              format(args.chain_id, args.pdb_file))
        sys.exit(1)

    try:
        active_list = [int(res) for res in args.active_list.split(',')]
        act_atoms = [a.get_coord() for a in atom_list if a.parent.id[1] in active_list]
    except:
        print('The list of active residues must be provided as a comma-separated list of integers')
        sys.exit(1)

    try:
        if args.surface_list:
            surface_list = [int(res) for res in args.active_list.split(',')]
        else:
            surface_list = get_surface_resids(s)
    except Exception as e:
        print("There was an error while calculating surface residues: {}".format(e))
        sys.exit(1)

    ns = NeighborSearch(atom_list)
    neighbors = []
    for a in act_atoms:
        neighbors.append(ns.search(a, 6.5, "R"))  # HADDOCK used 6.5A as default

    passive_list = set()
    for n in neighbors:
        for r in n:
            passive_list.add(r.id[1])
    tmp = passive_list & set(surface_list)
    passive_list = tmp - set(active_list)
    print(' '.join([str(r) for r in sorted(passive_list)]))

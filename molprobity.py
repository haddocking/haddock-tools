#!/usr/bin/env python
# Copyright 2013, Joao Rodrigues

"""
Python script to assign histidine protonation states using
Molprobity / Reduce.

Joao Rodrigues @ 2013
Adapted from Sjoerd's WHATIF code

syntax: molprobity.py <PDB-file>        
"""

import io
import os
import sys
import subprocess
import tempfile

def _check_molprobity_path(custom_path=None):
    """
    Tries to find 'reduce' executable in the system path.
    """
    
    if custom_path:
        if os.path.isfile(custom_path) and os.access(custom_path, os.X_OK):
            return custom_path
        else:
            raise Exception("Could not find path to 'reduce' executable: {0} does not exist or is not executable\n".format(custom_path))
    else:
        try:
            path = os.getenv('PATH').split(os.pathsep)
        except KeyError:
            raise Exception("Could not find path to 'reduce' executable: environment variable PATH not defined\n")
        for directory in path:
            if not os.path.exists(directory):
                continue
            if 'reduce' in os.listdir(directory):
                reduce_path = os.path.join(directory, 'reduce')
                if os.path.isfile(reduce_path) and os.access(reduce_path, os.X_OK):
                    return reduce_path
                else:
                    raise Exception("Found 'reduce' but it is either not executable or a directory.. ({0})\n".format(reduce_path))
    raise Exception("Could not find path to 'reduce' executable: Are you sure it is installed?\n")

def run_molprobity(pdbdata, molprobity_executable=None):
    """
    Reads a PDB file and outputs the corrected structure and a dictionary with protonation states.
    Expects either an open file handle or a string with a PDB formatted structure.
    
    Option strip_header removes all lines not starting with ATOM, TER, END, etc.. (check PDB format)
    """

    reduce_exec = _check_molprobity_path(molprobity_executable)
    cmd_string = [reduce_exec, '-build', '-Xplor', '-quiet']
    
    # File Handle vs Data String
    if isinstance(pdbdata, io.TextIOBase):
        cmd_stdin = pdbdata.read()
    else:
        cmd_stdin = pdbdata

    # Temporary File for Reduce
    tmp_file = tempfile.NamedTemporaryFile(mode="w+t")
    tmp_file.write(cmd_stdin)
    tmp_file.flush() # Force write to file otherwise might be incomplete

    cmd_string.append(tmp_file.name)

    try:
        process_handle = subprocess.Popen(cmd_string, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, close_fds=False)
    except Exception as e:
        raise Exception("There was an error running the 'reduce': {0}\n({1})\n".format(' '.join(cmd_string), e))    

    p_stdout, p_stderr = process_handle.communicate()

    return p_stdout, p_stderr

def analyze_protonation_state(pdbdata,pdbname="the PDB file"):      
    
    his_db = {}
    hisprotonatoms = {" HD1": 'd1', " HD2": 'd2', " HE1": 'e1', " HE2": 'e2'}

    # Build Histidine 'Database'
    # Assign protonation states based on presence of atoms in the PDB file
    for l in pdbdata.splitlines():
        if not l.startswith('ATOM'):
            continue

        aname = l[12:16]
        resn = l[17:20]
        resi = int(l[22:26])
        if resn != "HIS":
            continue
        elif resn == "HIS" and aname in hisprotonatoms:
            if resi not in his_db:
                his_db[resi] = {}
            currhis = his_db[resi]
            histidine_state = hisprotonatoms.get(aname)
            currhis[histidine_state] = True
    
    # Decide on Protonation State for CNS/HADDOCK
    ret = []
    for resi in his_db:
        his = his_db[resi]
        dcount = his.get('d1', 0) + his.get('d2', 0)
        ecount = his.get('e1', 0) + his.get('e2', 0)
        total_count = dcount + ecount
        if total_count == 4:
            ret.append(dict(resid=resi,state="HIS+"))
        elif total_count == 3:
            if dcount ==2:
                ret.append(dict(resid=resi,state="HISD"))
            else:
                ret.append(dict(resid=resi,state="HISE"))
        else:
            raise Exception("Molprobity could not guess the protonation state of histidine {0:n} in {1}: dcount={2}, ecount={3}".format(resi, pdbname, dcount, ecount))
    return ret

if __name__ == "__main__":
    # Quick and dirty usage example
    if not sys.argv[1:]:
        print("usage: {0} <pdb file> <pdb file 2> ...".format(sys.argv[0]))
        sys.exit(1)
    for ppath in sys.argv[1:]:
        print("## Executing Reduce to assign histidine protonation states")
        print("## Input PDB: {0} ".format(ppath))
        open_fhandle = open(ppath)
        hadded, process_error = run_molprobity(open_fhandle)

        # Write Optimized PDB, without USER nor Hs
        pname = os.path.basename(ppath)[:-4]
        fout = open("{0}_optimized.pdb".format(pname), 'w')
        for line in hadded.splitlines():
            fout.write(line+'\n')
        fout.close()


        ret = analyze_protonation_state(hadded)
        open_fhandle.close()
        #print '\n'.join([pdbx_line[:78] for pdbx_line in hadded.splitlines() if not pdbx_line.startswith('USER')])
        #sorted_histidines = sorted(ret, key=lambda x: x['resid'])
        sorted_states = sorted(ret, key=lambda x: (x['state'], x['resid']))
        for his in sorted_states:
            print("HIS ( {0:n} )\t-->\t{1}".format(his['resid'], his['state']))


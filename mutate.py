#!/usr/bin/python

"""
(dummy-)Mutates residues on a PDB-formatted structure. 
Reads mutations from a list with the following format:

	resi resn chain

where resn is the mutated residue.

Example, for mutating residue 18 of chain C to an Alanine:
	18 ALA C

HADDOCK will then reconstruct the residue according to its topology.
"""

def mutate(structure, resi, chain, mutation):

	mutated_structure = []
	ori_res = None
	atomlist = ["CA", "N", "O", "C"]

	for line in structure:
		if line[0:4] == 'ATOM':
			s_resi = line[22:26].strip()
			s_chain = line[21].strip()
			s_atom = line[12:16].strip()
			if s_resi == resi and s_chain == chain and s_atom in atomlist:
				ori_res = line[17:20].strip()
				line = line[0:17]+mutation+line[20:]
				mutated_structure.append(line)
			elif s_resi == resi and s_chain == chain and s_atom not in atomlist: 
				continue
			else:
				mutated_structure.append(line)
		else: 
			mutated_structure.append(line)
	return (mutated_structure, ori_res)

if __name__ == "__main__":
	
	import sys
	import os

	USAGE = "\nusage: python {0} <pdb file> <mutation list>\n\n".format(sys.argv[0])
	USAGE+= "(dummy-)Mutates residues on a PDB-formatted structure.\n"
	USAGE+= "Reads mutations from a list with the following format:\n"
	USAGE+= "\tresi resn chain\n"
	USAGE+= "Example, for mutating residue 18 of chain C to an Alanine:\n"
	USAGE+= "\t18 ALA C\n"

	if len(sys.argv[1:]) != 2:
		sys.stderr.write(USAGE)
		sys.exit(1)
		
	s_path = sys.argv[1]
	ml_path = sys.argv[2]


	if not os.path.exists(s_path) or not os.path.exists(ml_path):
		sys.stderr.write('Input files not found\n')
		sys.exit(1)

	with open(s_path) as pdb_fhandle:
			structure = [l for l in pdb_fhandle]
	with open(ml_path) as mut_fhandle:
			mut_list = [l.split() for l in mut_fhandle if l.strip()]

	for mutation in mut_list:
		if len(mutation) == 3:
			resi, resn, chain = mutation
		elif len(mutation) == 2:
			resi, resn = mutation
			chain = ""
		else:
			sys.stderr.write('Unrecognized mutation format in line {0}'.format(lineno))
			continue

		#print "Mutating residue {0}{1} to {2}".format(resi, chain, resn)
		mutant, original_residue = mutate(structure, resi, chain, resn)
		if not original_residue:
			sys.stderr.write('Residue not found: {0}{1}\n'.format(resi, chain))
			continue

		if chain.strip():
			m_file = open('{0}_{1}_{2}{3}{4}.pdb'.format(os.path.basename(s_path).split('.')[0], chain, original_residue, resi, resn), 'w')
			print '{0}_{1}_{2}{3}{4}.pdb'.format(os.path.basename(s_path).split('.')[0], chain, original_residue, resi, resn)
		else:
			m_file = open('{0}_{1}{2}{3}{4}.pdb'.format(os.path.basename(s_path).split('.')[0], chain, original_residue, resi, resn), 'w')
			print '{0}_{1}{2}{3}{4}.pdb'.format(os.path.basename(s_path).split('.')[0], chain, original_residue, resi, resn)
		m_file.write(''.join(mutant))
		m_file.close()

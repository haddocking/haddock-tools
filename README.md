
 haddock-tools 
================================================

Set of useful HADDOCK utility scripts 

About 
--------- 
This is a collection of scripts useful for pre- and post-processing and analysis for HADDOCK runs.  
Requests for new scripts will be taken into consideration, depending on the effort and general usability of the script.  


Installation
------------
Download the zip archive or clone the repository with git. This last is the recommended option as it
is then extremely simple to get updates.

```bash
# To download
git clone https://github.com/haddocking/haddock-tools

# To compile the executables
cd haddock-tools
make

# To update
cd haddock-tools && git pull origin master
```

Scripts
------------


## Restraints-related

#### active-passive_to_ambig.py
A python script to create ambiguous interaction restraints for use in HADDOCK based on list of active and passive residues (refer to the [HADDOCK software page](http://www.bonvinlab.org/software/haddock2.2/haddock.html) for more information)

Usage:
```bash
     ./active-passive_to_ambig.py <active-passive-file1> <active-passive-file2>
```

where <active-passive-file> is a file consisting of two space-delimited lines with
the first line active residues numbers and the second line passive residue numbers. One file per input structure should thus be provided.


#### restrain_bodies.py
A python script to creates distance restraints to lock several chains together. 
Useful to avoid unnatural flexibility or movement due to 
sequence/numbering gaps during the refinement stage of HADDOCK.

Usage:
```bash
./restrain_bodies.py [-h] [--exclude EXCLUDE [EXCLUDE ...]] [--verbose] structures [structures ...]

  positional arguments:
    structures            PDB structures to restraint

  optional arguments:
    -h, --help            show this help message and exit
    --exclude EXCLUDE [EXCLUDE ...], -e EXCLUDE [EXCLUDE ...] Chains to exclude from the calculation
    --verbose, -v
```

#### restrain_ligand.py
Calculates distances between neighboring residues of a ligand molecule and produces a set of
unambiguous distance restraints for HADDOCK to keep it in place during semi-flexible refinement.
Produces, at most, one restraint per ligand atom.

Usage:
```bash
./restrain_ligand.py [-h] -l LIGAND [-p] pdbf

positional arguments:
  pdbf                  PDB file

optional arguments:
  -h, --help            show this help message and exit
  -l LIGAND, --ligand LIGAND
                        Ligand residue name
  -p, --pml             Write Pymol file with restraints
```

#### haddock_tbl_validation
The validate_tbl.py script in that directoy will check the correctness of your restraints (CNS format) for HADDOCK.

Usage:
```bash
usage: python validate_tbl.py [-h] [--pcs] file

This script validates a restraint file (*.tbl).

positional arguments:
  file        TBL file to be validated

  optional arguments:
    -h, --help  show this help message and exit
    --pcs       PCS mode
```


------------
## PDB-related

#### contact-segid
A c++ program to calculate all heavy atom interchain contacts (where the chain identification is taken from the segid) within a given distance cutoff in Angstrom.

Usage:
```bash
   contact-segid <pdb file> <cutoff>
```

#### contact-chainID
A c++ program to calculate all heavy atom interchain contacts (where the chain identification is taken from the chainID) within a given distance cutoff in Angstrom.

Usage:
```bash
   contact-chainID <pdb file> <cutoff>
```

#### molprobity.py
A python script to predict the protonation state of Histidine residues for HADDOCK. It uses molprobity for this, calling the [Reduce](http://kinemage.biochem.duke.edu/software/reduce.php) software which should in the path.

Usage:
```bash
    ./molprobity.py <PDBfile>
```

Example:
```bash
./molprobity.py 1F3G.pdb
## Executing Reduce to assign histidine protonation states
## Input PDB: 1F3G.pdb
HIS ( 90 )	-->	HISD
HIS ( 75 )	-->	HISE
```
An optimized file is also written to disk, in this example it would be called ```1F3G_optimized.pdb```. 


#### pdb_blank_chain
Simple perl script to remove the chainID from a PDB file

Usage:
```bash
    pdb_blank_chain inputfile > outputfile
```

#### pdb_blank_segid
Simple perl script to remove the segid from a PDB file

Usage:
```bash
    pdb_blank_segid inputfile > outputfile
```

#### pdb_blank_chain-segid
Simple perl script to remove both the chainID and segid from a PDB file

Usage:
```bash
    pdb_blank_chain-segid inputfile > outputfile
```

#### pdb_chain-to-segid
Simple perl script to copy the chainID to the segid in a PDB file

Usage:
```bash
    pdb_chain-to-segid inputfile > outputfile
```

#### pdb_segid-to-chain
Simple perl script to copy the segid to the chainID in a PDB file

Usage:
```bash
    pdb_segid-to-chain inputfile > outputfile
```

#### pdb_chain-segid
Simple perl script to copy the chainID to segid in case the latter is empty (or vice-verse) in a PDB file

Usage:
```bash
    pdb_chain-segid inputfile > outputfile
```

#### pdb_setchain
Simple perl script to set the chainID in a PDB file

Usage:
```bash
     pdb_setchain -v CHAIN=chainID inputfile > outputfile
```

#### joinpdb
Simple perl script to concatenate separate single structure PDB files into a multi-model PDB file.
Usage:
```bash
     joinpdb  -o outputfile  [inputfiles]

    where inputfiles are a list of PDB files to be concatenated
```

#### pdb_mutate.py
A python script to mutate residues for HADDOCK. A mutation list file is used as input, and the output is/are corresponding PDB file(s) of mutant(s). The format of mutation in the mutation list file is "PDBid ChainID ResidueID ResidueNameWT ResidueNameMut".


Usage:
```bash
    ./pdb_mutate.py <mutation list file>
```

Example:
```bash
./pdb_mutate.py mut_1A22.list

## In  mut_1A22.list, the residue 14, 18 and 21 in chain A will be mutated to ALA:
## 1A22.pdb A 14 MET ALA
## 1A22.pdb A 18 HIS ALA
## 1A22.pdb A 21 HIS ALA
```

#### pdb_strict_format.py
A python script to check format of PDB files with respect to HADDOCK format rules. A PDB file is used as input, and the output is a console message if an error or a warning is triggered by a bad formmated line. The script uses wwPDB format guidelines [wwwPDB guidelines](http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html) and check resid against a list of known ligands and amino-acids recognized by HADDOCK.


Usage:
```bash
./pdb_strict_format.py [-h] [-nc] pdb

This script validates a PDB file (*.pdb).

positional arguments:
  pdb                PDB file

optional arguments:
  -h, --help         show this help message and exit
  -nc, --no_chainid  Ignore empty chain ids
```

License
---------

Apache Licence 2.0

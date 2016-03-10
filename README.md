================================================
============== haddock-tools ===================
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

# To update
cd haddock-tools && git pull origin master
```

Scripts
------------


## Restraints-related

#### active-passive_to_ambig.py
A python script to create ambiguous interaction restraints for use in HADDOCK based on list of active and passive residues (refer to the [HADDOCK software page](http://www.bonvinlab.org/software/haddock2.2/haddock.html) for more infmation)

Usage:
```bash
     python active-passive_to_ambig.py <active-passive-file1> <active-passive-file2>
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
A python script to predict the protonation state of Histidine residues for HADDOCK. It uses molprobity for this, calling the reduce software which should in the path.

Usage:
```bash
    molprobity.py <PDBfile>
```

Example:
```bash
molprobity.py 1F3G.pdb
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



License
---------

Apache Licence 2.0

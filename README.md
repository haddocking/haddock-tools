
 haddock-tools 
================================================

Set of useful HADDOCK utility scripts, which requires Python 3.7+. 

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

#### passive_from_active.py
A python script to obtain a list of passive residues providing a PDB file and a list of active residues.
This will automatically calculate a list of surface residues from the PDB to filter out buried residues except if
a surface list is provided.
By default, neighbors of the active residues are searched within 6.5 Angstroms and surface residues are residues whose
relative side chain accessibility or main chain accessibility is above 15%.

Requirements:
* [FreeSASA](https://freesasa.github.io/)

`pip install freesasa`
* [BioPython](http://biopython.org/)

`pip install biopython`

Usage:
```bash
./passive_from_active.py [-h] [-c CHAIN_ID] [-s SURFACE_LIST]
                              pdb_file active_list

positional arguments:
  pdb_file              PDB file
  active_list           List of active residues IDs (int) separated by commas

optional arguments:
  -h, --help            show this help message and exit
  -c CHAIN_ID, --chain-id CHAIN_ID
                        Chain id to be used in the PDB file (default: All)
  -s SURFACE_LIST, --surface-list SURFACE_LIST
                        List of surface residues IDs (int) separated by commas
```

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

### calc-accessibility.py

```bash
$ python3 haddock-CSB-tools/calc-accessibility.py -h                                                                                                                                                                                                               [17:06:52]
usage: calc-accessibility.py [-h] [--log-level {DEBUG,INFO,WARNING,ERROR,CRITICAL}] [--cutoff CUTOFF] pdb_input

positional arguments:
  pdb_input             PDB structure

optional arguments:
  -h, --help            show this help message and exit
  --log-level {DEBUG,INFO,WARNING,ERROR,CRITICAL}
  --cutoff CUTOFF       Relative cutoff for sidechain accessibility

$ python3 haddock-CSB-tools/calc-accessibility.py complex_1w.pdb --cutoff 0.4                                                                                                                                                                                      [17:10:51]
02/11/2020 17:10:57 L157 INFO - Calculate accessibility...
02/11/2020 17:10:57 L228 INFO - Chain: A - 115 residues
02/11/2020 17:10:57 L228 INFO - Chain: B - 81 residues
02/11/2020 17:10:57 L234 INFO - Applying cutoff to side_chain_rel - 0.4
02/11/2020 17:10:57 L244 INFO - Chain A - 82,83,84,85,86,87,88,90,91,94,95,98,99,102,104,106,109,113,116,117,118,122,128,129,130,132,139,141,144,145,148,149,150,151,153,156,158,160,162,163,167,168,169,170,171,173,174,175,176,178,179,180,181,183,184,186,188,194,196
02/11/2020 17:10:57 L244 INFO - Chain B - 1,2,4,5,8,11,12,15,18,21,23,24,25,26,27,30,31,33,34,37,38,41,43,44,45,46,47,50,63,64,67,69,70,73,74,76,77,78,79,80,81
```

### create_cif.py

Converts the `cluster*.pdb` files in a run directory to [IHM mmCIF](https://mmcif.wwpdb.org/dictionaries/mmcif_ihm.dic/Index/) format

Warning: Limited functionally, still work in progress! Tested for hetero-complexes with ambig restraints.

Needs `ihm` and `biopython`, install it with
```bash
$ pip install ihm --install-option="--without-ext"
$ pip install biopython
```

```
$ python3 create_cif.py -h
usage: create_cif.py [-h] run_directory

positional arguments:
  run_directory  Location of the uncompressed run, ex:
                 /home/rodrigo/runs/47498-protein-protein

optional arguments:
  -h, --help     show this help message and exit


$ python3 haddock-CSB-tools/create_cif.py ~/projects/cif_parser/47518-cif
[23/02/2021 13:45:25] INFO Converting the cluster*.pdb structures to .cif
[23/02/2021 13:45:25] INFO Looking for models in /Users/rodrigo/projects/cif_parser/47518-cif
[23/02/2021 13:45:25] INFO Found 4 structures
[23/02/2021 13:45:25] INFO Looking for the tblfile field in /Users/rodrigo/projects/cif_parser/47518-cif/job_params.json
[23/02/2021 13:45:25] INFO tblfile field found, extracting information
[23/02/2021 13:45:25] INFO Converting /Users/rodrigo/projects/cif_parser/47518-cif/cluster1_1.pdb
[23/02/2021 13:45:26] INFO Saving as cluster1_1.cif
[23/02/2021 13:45:26] INFO Converting /Users/rodrigo/projects/cif_parser/47518-cif/cluster1_2.pdb
[23/02/2021 13:45:26] INFO Saving as cluster1_2.cif
[23/02/2021 13:45:27] INFO Converting /Users/rodrigo/projects/cif_parser/47518-cif/cluster1_3.pdb
[23/02/2021 13:45:27] INFO Saving as cluster1_3.cif
[23/02/2021 13:45:28] INFO Converting /Users/rodrigo/projects/cif_parser/47518-cif/cluster1_4.pdb
[23/02/2021 13:45:29] INFO Saving as cluster1_4.cif
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

#### param_to_json.py
A python script to transform a haddockparam.web file into a JSON structure. It is possible to use it as a class and then access
extra functions like: `change_value(key, value)` ; `update(subdict_to_replace)` ; `dump_keys()` ; `get_value(key)` ; `write_json()`

Usage:
```bash
./param_to_json.py [-h] [-o OUTPUT] [-g GET] [-e [EXAMPLE]] web

This script parses a HADDOCK parameter file (*.web) and transforms it to JSON
format. It also allows to change a parameter of the haddockparam.web

positional arguments:
  web                   HADDOCK parameter file

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT
                        Path of JSON output file
  -g GET, --get GET     Get value of a particular parameter
  -e [EXAMPLE], --example [EXAMPLE]
                        Print an example
```

License
---------

Apache Licence 2.0

haddock-tools ================================================ Set of useful HADDOCK utility scripts About --------- This is a collection of scripts useful for pre- and post-processing and analysis for HADDOCK runs.  Requests for new scripts will be taken into consideration, depending on the effort and general usability of the script.  
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

### active-passive_to_ambig.py
A python script to create ambiguous interaction restraints for use in HADDOCK based on list of active and passive residues (refer to the [HADDOCK software page](http://www.bonvinlab.org/software/haddock2.2/haddock.html) for more infmation)

Usage:
```bash
     python active-passive_to_ambig.py <active-passive-file1> <active-passive-file2>
```

where <active-passive-file> is a file consisting of two space-delimited lines with
the first line active residues numbers and the second line passive residue numbers. One file per input structure should thus be provided.

### molprobity.py
A python script to predict the protonation state of Histidine residues for HADDOCK. It uses molprobity for this, calling the reduce software which should in the path.

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


### restrain_bodies.py
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

License
---------

Apache Licence 2.0

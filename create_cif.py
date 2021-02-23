import os
import logging
import pathlib
import argparse
import json
import re
from itertools import groupby, count
import ihm
import ihm.location
import ihm.dataset
import ihm.representation
import ihm.restraint
import ihm.protocol
import ihm.model
import ihm.dumper
from Bio.PDB import PDBParser as BioParser, Selection, NeighborSearch

VALID_BASES = ['DA', 'DC', 'DG', 'DT', 'A', 'C', 'G', 'U', 'DJ']
VALID_AMINOACIDS = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE',
                    'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL', 'ACE', 'ALY', 'ASH', 'CFE', 'CHX', 'CSP', 'CTN', 'CYC',
                    'CYF', 'CYM', 'DDZ', 'DUM', 'GLH', 'HLY', 'HYP', 'M3L', 'MLY', 'MLZ', 'MSE', 'NEP', 'PNS', 'PTR',
                    'QSR', 'SEP', 'TOP', 'TYP', 'TYS']


AMINOACID_CODES = [
    # 20 canonical amino acids
    ('CYS', 'C'), ('ASP', 'D'), ('SER', 'S'), ('GLN', 'Q'),
    ('LYS', 'K'), ('ILE', 'I'), ('PRO', 'P'), ('THR', 'T'),
    ('PHE', 'F'), ('ASN', 'N'), ('GLY', 'G'), ('HIS', 'H'),
    ('LEU', 'L'), ('ARG', 'R'), ('TRP', 'W'), ('ALA', 'A'),
    ('VAL', 'V'), ('GLU', 'E'), ('TYR', 'Y'), ('MET', 'M'),
    # Non-canonical amino acids
    ('MSE', 'M'), ('SOC', 'C')]

BASE_CODES = [
    # 5 canonical nucleic acids
    ('ADE', 'A'), ('GUA', 'G'), ('THY', 'T'), ('CYT', 'C'), ('URI', 'U'),
    ('DA', 'A'), ('DG', 'G'), ('DT', 'T'), ('DC', 'C'),
    ('A', 'A'), ('G', 'G'), ('T', 'T'), ('C', 'C'), ('U', 'U')]

NARES3 = ['---', 'CYT', 'THY', 'GUA', 'ADE', 'URI']
NARES2 = ['--', 'DC', 'DT', 'DG', 'DA', 'DU']

THREE_TO_ONE = dict(AMINOACID_CODES + BASE_CODES)
THREE_TO_TWO = dict([(three, two) for three, two in zip(NARES3, NARES2)])


class Restraint:
    """Represent a CNS restraint between two residues as a tuple"""

    def __init__(self, partner1_chainid, partner1_resid, partner2_chainid, partner2_resid, distance, groupid):
        #: Chain ID of 1st residue
        self.partner1_chainid = partner1_chainid
        #: Residue ID of 1st residue
        self.partner1_resid = partner1_resid
        #: Chain ID of 2nd residue
        self.partner2_chainid = partner2_chainid
        #: Residue ID of 2nd residue
        self.partner2_resid = partner2_resid
        #: Distance restraint as a list [upper_d, lower_d, target]
        self.distance = distance
        self.groupid = groupid

    def __repr__(self):
        return f'Restraint ({self.groupid}): {self.partner1_chainid}.{self.partner1_resid} <-> '\
               f'{self.partner2_chainid}.{self.partner2_resid} : {self.distance}'


def _parse_block(block, group=1):
    """
    Parse an `assign` block containing one or multiple selection and a distance value

    :param str block: Block of assign statement starting with assign and finishing with distance value
    :param int group: ID of the first group to be found
    :return: List of restraints based on the selections
    :rtype: list of :class:`~tools.TBLHandler.Restraint`

    """
    block = "".join(block)
    distance = re.findall(r'(\d+\.\d+) (\d+\.\d+) (\d+\.\d+)', block)
    distance = tuple([float(d) for d in distance[0]])
    selections = [s.strip('() ') for s in re.findall(r'\([^(]*\)', block)]
    partner1 = selections[0]
    partners = selections[1:]
    resid = r'resid +(\d+)'
    chainid = r'segid +([A-Z])'
    partner1_chainid = re.findall(chainid, partner1)[0]
    partner1_resid = int(re.findall(resid, partner1)[0])
    restraints = []
    for p in partners:
        partner2_chainid = re.findall(chainid, p)[0]
        partner2_resid = int(re.findall(resid, p)[0])
        restraints.append(Restraint(partner1_chainid, partner1_resid, partner2_chainid, partner2_resid,
                                    distance, group))
    return restraints


def extract_restraints(tbl_content):
    """
    Extract restraints from a TBL file content

    :param list tbl_content: List of TBL file lines to be parsed
    :return: List of all restraints present in the file
    :rtype: list of :class:`~tools.TBLHandler.Restraint`

    """
    restraints = []
    current_assign = []
    num_open = 0
    num_close = 0
    assign_block_found = False
    for line in tbl_content.split(os.linesep):
        line = line.strip()
        if line and line[0] != "!":
            if assign_block_found:
                if num_open != num_close or len(current_assign) == 1:
                    current_assign.append(line)
                    num_open += line.count("(")
                    num_close += line.count(")")
                else:
                    restraints.extend(_parse_block(current_assign))
                    current_assign = []
                    num_open = 0
                    num_close = 0
                    assign_block_found = False
            if line.startswith("assign"):
                assign_block_found = True
                current_assign.append(line)
                num_open += line.count("(")
                num_close += line.count(")")
    # Process the last line
    restraints.extend(_parse_block(current_assign))
    return restraints


def get_interface_ranges(interface):
    """
    Convert lists of interface residues to lists of interface residue ranges

    :param interface: Dictionary of interface residues per chain IDs (e.g. {'A':[1,2,3], 'B': [10,11,12], ...}
    :return: interface_ranges: Dictionary of interface residues as ranges (e.g. {'A': [(1,3), (8,10)], 'B':[(10,12)]}
    :rtype: dict

    """
    def as_range(iterable):  # not sure how to do this part elegantly
        l = list(iterable)
        if len(l) > 1:
            return (l[0], l[-1])
        else:
            return (l[0], l[0])

    interface_ranges = {}
    for p in interface:
        interface_ranges[p] = [as_range(g) for _, g in groupby(interface[p], lambda n, c=count(): n - next(c))]
    return interface_ranges


def get_interface_residues(pdb_file, radius=10.0):
    """Return a list of interacting residues based on accessibility and distance to partner

    :param str pdb_file: PDB file path
    :param int,float radius: Maximum distance to be considered to define a residue as interacting with another chain
    :return: interface: Dictionary of interface residues per chain IDs (e.g. {'A':[1,2,3], 'B': [10,11,12], ...}
    :rtype: dict

    """
    p = BioParser(QUIET=True)
    s = p.get_structure('pdb', pdb_file)
    if sum(1 for _ in s.get_chains()) < 2:
        logging.error("Less than 2 chains have been detected in your PDB, no interface can be extracted.")
        # raise PDBParsingError("Only one chain has been found, no interface can be extracted.")
    m = s[0]
    all_atoms = Selection.unfold_entities(m, 'A')
    # Unfold atoms for NeighborSearch algorithm to work
    ns = NeighborSearch(all_atoms)
    interface = ns.search_all(radius, "R")
    # Filter redundant residues
    buffer = dict([(ch_id.id, []) for ch_id in m.get_chains()])
    for r in interface:
        if r[0].parent.id != r[1].parent.id:
            if r[0].id[1] not in buffer[r[0].parent.id]:
                buffer[r[0].parent.id].append(r[0].id[1])
            if r[1].id[1] not in buffer[r[1].parent.id]:
                buffer[r[1].parent.id].append(r[1].id[1])
    interface = buffer
    for ch_id in m.get_chains():
        interface[ch_id.id] = sorted(interface[ch_id.id])
    return interface


def _update_peptide_alphabet():
    """Update the IHM mmCIF alphabet with modified amino-acids and ligands supported by HADDOCK

    :return dict_extension: A dictionary extension containing the information for new amino-acids and ligands
    :rtype: dict

    """
    dict_extension = {}
    dict_extension['ACE'] = ihm.LPeptideChemComp(id='ACE', code='ACE', code_canonical='n', name='Nterm acetyl group',
                                                 formula='C2 H3 O2')
    dict_extension['ALY'] = ihm.LPeptideChemComp(id='ALY', code='ALY', code_canonical='K', name='Acetylated LYS',
                                                 formula='C8 H14 O2 N2')
    dict_extension['ASH'] = ihm.LPeptideChemComp(id='ASH', code='ASH', code_canonical='D', name='Protonated ASP',
                                                 formula='C4 H5 O3 N1')
    dict_extension['CFE'] = ihm.LPeptideChemComp(id='CFE', code='CFE', code_canonical='C',
                                                 name='CYS with iron sulfur cluster', formula='C3 H4 O1 N1 S3 F2')
    dict_extension['CSP'] = ihm.LPeptideChemComp(id='CSP', code='CSP', code_canonical='C', name='Phosphorylated CYS',
                                                 formula='C3 H4 O4 N1 S1 P1')
    dict_extension['CTN'] = ihm.LPeptideChemComp(id='CTN', code='CTN', code_canonical='c', name='Cterm amide group',
                                                 formula='H2 N')
    dict_extension['CYF'] = ihm.LPeptideChemComp(id='CYF', code='CYF', code_canonical='C',
                                                 name='CYS without the sulfur H',
                                                 formula='C3 H4 O1 N1 S1')
    dict_extension['CYM'] = ihm.LPeptideChemComp(id='CYM', code='CYM', code_canonical='C', name='CYS with MTSL',
                                                 formula='C12 H5 O2 N2 S2')
    dict_extension['DDZ'] = ihm.LPeptideChemComp(id='DDZ', code='DDZ', code_canonical='A', name='3,3,-dihydroxy ALA',
                                                 formula='C2 H5 O2 N1')
    dict_extension['GLH'] = ihm.LPeptideChemComp(id='GLH', code='GLH', code_canonical='Q', name='Protonated GLU',
                                                 formula='C5 H7 O3 N1')
    dict_extension['HEB'] = ihm.LPeptideChemComp(id='HEB', code='HEB', code_canonical='h', name='Heme B',
                                                 formula='C34 H30 O4 N4 F1')
    dict_extension['HEC'] = ihm.LPeptideChemComp(id='HEC', code='HEC', code_canonical='h', name='Heme C',
                                                 formula='C34 H32 O4 N4 F1')
    dict_extension['HYP'] = ihm.LPeptideChemComp(id='HYP', code='HYP', code_canonical='P', name='4R-hydroxyproline',
                                                 formula='C5 H7 O2 N1')
    dict_extension['M3L'] = ihm.LPeptideChemComp(id='M3L', code='M3L', code_canonical='K', name='Trimethyl LYS',
                                                 formula='C9 H21 O1 N2')
    dict_extension['MLY'] = ihm.LPeptideChemComp(id='MLY', code='MLY', code_canonical='K', name='Dimethyl LYS',
                                                 formula='C8 H18 O1 N2')
    dict_extension['MLZ'] = ihm.LPeptideChemComp(id='MLZ', code='MLZ', code_canonical='K', name='Monomethyl LYS',
                                                 formula='C7 H15 O1 N2')
    dict_extension['NEP'] = ihm.LPeptideChemComp(id='NEP', code='NEP', code_canonical='H', name='NE phosphorylated HIS',
                                                 formula='C6 H7 O4 N3 P1')
    dict_extension['PTR'] = ihm.LPeptideChemComp(id='PTR', code='PTR', code_canonical='Y', name='O-Phosphotyrosine',
                                                 formula='C9 H10 O6 N1 P3')
    dict_extension['SEP'] = ihm.LPeptideChemComp(id='SEP', code='SEP', code_canonical='S', name='Phosphorylated SER',
                                                 formula='C3 H4 O5 N1 P1')
    dict_extension['TIP'] = ihm.LPeptideChemComp(id='TIP', code='TIP', code_canonical='w', name='Water model',
                                                 formula='H2O')
    dict_extension['TOP'] = ihm.LPeptideChemComp(id='TOP', code='TOP', code_canonical='T', name='Phosphorylated THR',
                                                 formula='C4 H7 O5 N1 P1')
    dict_extension['TYP'] = ihm.LPeptideChemComp(id='TYP', code='TYP', code_canonical='Y', name='Phosphorylated TYR',
                                                 formula='C9 H8 O5 N1 P1')
    dict_extension['TYS'] = ihm.LPeptideChemComp(id='TYS', code='TYS', code_canonical='Y', name='Sulfonated TYR',
                                                 formula='C9 H8 O5 N1 S1')
    dict_extension['PNS'] = ihm.LPeptideChemComp(id='PNS', code='PNS', code_canonical='S',
                                                 name='Phosphopantetheine SER',
                                                 formula='C11 H23 O7 N2 S1 P1')
    return dict_extension


def create_ihm_mmcif(pdb_file, params, interface, restraints=None):
    """Create a IHM mmCIF compatible file from a PDB file and a list of HADDOCK parameters and interface residue IDs

    :param str pdb_file: PDB file path
    :param dict params: HADDOCK parameter file as a dictionary
    :param list interface: List of interface residue IDs
    :param restraints: List of Restraints used during the docking process
    :type: list of :class:`~tools.TBLHandler.Restraint`

    """
    # Load run parameters
    # with open(param_path) as o:
    #     params = json.load(o)
    # Stores main information - chain IDs, active residues, passives residues
    partners = {}
    for i, p in params['partners'].items():
        partners[p['segid']] = {}
        partners[p['segid']]['activeres'] = p['activereslist']
        partners[p['segid']]['passiveres'] = p['passivereslist']
        partners[p['segid']]['moleculetype'] = p['moleculetype']

    # First, we create a system, which contains everything we know about the
    # modeling. A single mmCIF file can contain multiple Systems, but in most
    # cases we use just one:
    system = ihm.System()

    p = BioParser(QUIET=True)
    try:
        s = p.get_structure('HADDOCK model', pdb_file)
    except Exception as e:
        logging.error("Error while loading the PDB: {}".format(e))
        raise

    # Sequences are concatenated string, we ignore HOH
    seqs_id = {}
    seqs = {}

    # Update LPeptideAlphabet to add modified aa supported by HADDOCK
    dict_extension = _update_peptide_alphabet()

    # We store sequences and match between author resids and mmCIF consecutive resids (1->n)
    for r in s.get_residues():
        if r.get_parent().id in partners and r.resname != "HOH":
            segid = r.get_parent().id
            if segid not in seqs_id:
                seqs_id[segid] = 1
            else:
                seqs_id[segid] += 1
            if segid not in seqs:
                seqs[segid] = {'seq': [], 'id_from_pdb': {}, 'id_from_seq': {}, 'hetero': False}

            seqs[segid]['id_from_pdb'][seqs_id[segid]] = r.id[1]
            seqs[segid]['id_from_seq'][r.id[1]] = seqs_id[segid]
            # Detect whether the residue is a non-standard aa (but supported by HADDOCK)
            if r.resname not in VALID_BASES and r.resname not in VALID_AMINOACIDS:
                pass
            # We only store one-letter code for standard aa
            try:
                if partners[segid]['moleculetype'] == "DNA":
                    seqs[segid]['seq'].append(THREE_TO_ONE[r.resname.strip()])
                else:
                    seqs[segid]['seq'].append(THREE_TO_ONE[r.resname.strip()])
            except KeyError as e:
                # Detect modified aa and supported ligands
                if partners[segid]['moleculetype'] not in ("DNA", "RNA") and r.resname.strip() in dict_extension:
                    seqs[segid]['seq'].append(dict_extension[r.resname.strip()])
                elif partners[segid]['moleculetype'] not in ("DNA", "RNA") and not seqs[segid]['hetero']:
                    seqs[segid]['hetero'] = True
                    seqs[segid]['seq'].append(r.resname.strip())
                logging.error(e)
            except Exception as e:
                logging.error(f"An issue occurred when trying to add residue {r}", e)
                raise Exception
            # seqs[segid].append(r.resname.strip())

    # Create different entities with specific dictionary depending on the unit type (Protein, RNA, DNA)
    entities = []
    for i, ch_id in enumerate(partners):
        if partners[ch_id]['moleculetype'] == 'DNA':
            alphabet = ihm.DNAAlphabet
        elif partners[ch_id]['moleculetype'] == 'RNA':
            alphabet = ihm.RNAAlphabet
        else:
            alphabet = ihm.LPeptideAlphabet
        if not seqs[ch_id]['hetero']:
            entities.append(ihm.Entity(seqs[ch_id]['seq'], alphabet=alphabet, description=f'Partner {i+1}'))
        else:
            entities.append(ihm.Entity([ihm.NonPolymerChemComp(seqs[ch_id]['seq'][0])], description=f'Partner {i+1}'))
    system.entities.extend(tuple(entities))

    # Create different asymetric units based on author sequence ids
    asym_units = {}
    for i, ch_id in enumerate(partners):
        asym_units[ch_id] = (ihm.AsymUnit(entities[i], auth_seq_id_map=seqs[ch_id]['id_from_pdb'], id=ch_id))
    system.asym_units.extend(tuple([unit for k, unit in asym_units.items()]))

    # Create asymetric unit ranges from interface residues
    asym_unit_ranges = []
    for ch_id, unit in asym_units.items():
        if not seqs[ch_id]['hetero']:
            for r in interface[unit.id]:
                asym_unit_ranges.append(ihm.AsymUnitRange(unit, seqs[unit.id]['id_from_seq'][r[0]],
                                                          seqs[unit.id]['id_from_seq'][r[1]]))

    # Create Model assembly
    modeled_assembly = ihm.Assembly([unit for k, unit in asym_units.items()], name='Modeled assembly')
    assemblies = []
    for ch_id, unit in asym_units.items():
        assemblies.append(ihm.Assembly((unit,), name=f'Subunit {unit.id}'))

    # In HADDOCK interface is flexible whereas the rest of the molecules are rigid
    # We need to create new asymmetric units to store the interface
    segm = [ihm.representation.AtomicSegment(unit, rigid=True) for ch_id, unit in asym_units.items()]
    segm += [ihm.representation.AtomicSegment(asym, rigid=False) for asym in asym_unit_ranges]
    rep = ihm.representation.Representation(segm)

    # Define dataset location (in auto mode we defined them as external file)
    pdb_datasets = []
    em_dataset = []
    for i, p in enumerate(partners):
        loc = ihm.location.InputFileLocation(pdb_file.parent / f'data/sequence/protein{i+1}.pdb')
        dataset = ihm.dataset.PDBDataset(loc)
        pdb_datasets.append(dataset)

    if params['em_rest']:
        loc = ihm.location.InputFileLocation(pdb_file.parent / f'data/sequence/protein{i+1}.pdb')
        dataset = ihm.dataset.EMDensityDataset(loc)
        em_dataset.append(dataset)

    all_datasets = ihm.dataset.DatasetGroup(pdb_datasets + em_dataset)
    # Extract restraints
    # TODO (un)ambiguous restraints
    poly_res_features = []
    if restraints:
        for r in restraints:
            a1 = asym_units[r.partner1_chainid]
            f = ihm.restraint.ResidueFeature([a1(seqs[r.partner1_chainid]['id_from_seq'][r.partner1_resid],
                                                 seqs[r.partner1_chainid]['id_from_seq'][r.partner1_resid])])
            if f.ranges[0] not in poly_res_features:
                system.orphan_features.append(f)
                poly_res_features.append(f.ranges[0])
            a2 = asym_units[r.partner2_chainid]
            if seqs[r.partner2_chainid]['hetero']:
                f2 = ihm.restraint.ResidueFeature([a2])
            else:
                f2 = ihm.restraint.ResidueFeature([a2(seqs[r.partner2_chainid]['id_from_seq'][r.partner2_resid],
                                                      seqs[r.partner2_chainid]['id_from_seq'][r.partner2_resid])])
            if f2.ranges[0] not in poly_res_features:
                system.orphan_features.append(f2)
                poly_res_features.append(f2.ranges[0])
    # We do not know from where the ambiguous restraints come from so we create a generic dataset that will be used
    # as reference for the ihm_derived_distance_restraint list
    loc = ihm.location.InputFileLocation(pdb_file.parent / 'data/distances/ambig.tbl')
    derived_dist_dataset = ihm.dataset.Dataset(loc)
    dist = ihm.restraint.UpperBoundDistanceRestraint(2.0)
    probability = 1 if not params["noecv"] else 1.0 / params["ncvpart"]
    restraints = []
    for feature in system.orphan_features:
        for feature2 in system.orphan_features:
            # We create restraints only between residues of different subunits (chain IDs)
            if feature.ranges[0].asym.id != feature2.ranges[0].asym.id and (feature, feature2) not in restraints \
                    and (feature2, feature) not in restraints:
                r = ihm.restraint.DerivedDistanceRestraint(dataset=derived_dist_dataset, feature1=feature,
                                                           feature2=feature2, distance=dist, probability=probability)
                system.restraints.append(r)
                restraints.append((feature, feature2))
    if params['em_rest']:
        em_rsr = ihm.restraint.EM3DRestraint(dataset=em_dataset[0], assembly=modeled_assembly)
        system.restraints.append(em_rsr)

    # Add protocol used
    if params['em_rest']:
        protocol = ihm.protocol.Protocol(name='HADDOCK-EM')
    else:
        protocol = ihm.protocol.Protocol(name='HADDOCK')

    protocol.steps.append(ihm.protocol.Step(
        assembly=modeled_assembly,
        dataset_group=all_datasets,
        method='Rigid-body minimization in HADDOCK (it0)',
        name='Rigid-body minimization',
        num_models_begin=0,
        num_models_end=params['structures_0'], multi_scale=True))
    protocol.steps.append(ihm.protocol.Step(
        assembly=modeled_assembly,
        dataset_group=all_datasets,
        method='Semi-flexible SA in HADDOCK (it1)',
        name='Simulated annealing',
        num_models_begin=0,
        num_models_end=params['structures_1'], multi_scale=True))
    protocol.steps.append(ihm.protocol.Step(
        assembly=modeled_assembly,
        dataset_group=all_datasets,
        method='Water refinement in HADDOCK (itw)',
        name='Refinement',
        num_models_begin=0,
        num_models_end=params['waterrefine'], multi_scale=True))

    atoms = dict([(ch_id.id, []) for ch_id in s.get_chains()])
    for c in s.get_chains():
        for r in c:
            for a in r:
                atoms[c.id].append((c.id, seqs[c.id]['id_from_seq'][r.id[1]], a.name[0], a.name, a.coord[0], a.coord[1],
                                    a.coord[2]))

    # Subclassing the IHM Model class and overriding the get_atoms method:
    class MyModel(ihm.model.Model):
        # Map our asym unit names A and B to IHM asym_unit objects:
        asym_unit_map = dict([(unit.id, unit) for ch_id, unit in asym_units.items()])

        def get_atoms(self):
            for ch_id in self.asym_unit_map:
                for asym, seq_id, type_symbol, atom_id, x, y, z in atoms[ch_id]:
                    yield ihm.model.Atom(asym_unit=self.asym_unit_map[asym], type_symbol=type_symbol, seq_id=seq_id,
                                         atom_id=atom_id, x=x, y=y, z=z)
    clus = [int(s) for s in os.path.basename(pdb_file) if s.isdigit()][0]
    best = [int(s) for s in os.path.basename(pdb_file) if s.isdigit()][1]
    if clus == best == 1:
        model = MyModel(assembly=modeled_assembly, protocol=protocol, representation=rep, name='Best scoring model')
    else:
        model = MyModel(assembly=modeled_assembly, protocol=protocol, representation=rep,
                        name=f'Best {best} of cluster {clus}')

    # Add cross-correlation value if EM map was used
    if params['em_rest']:
        with open(pdb_file) as o:
            for l in o:
                if 'cross-correlation' in l:
                    cc = float(l.split()[3])
                    break
        em_rsr.fits[model] = ihm.restraint.EM3DRestraintFit(cross_correlation_coefficient=cc)

    # Similar models can be grouped together. Here we only have a single model
    # in the group
    model_group = ihm.model.ModelGroup([model], name='All models')

    # Groups are then placed into states, which can in turn be grouped. In this
    # case we have only a single state:
    state = ihm.model.State([model_group])
    system.state_groups.append(ihm.model.StateGroup([state]))

    # Dump mmCIF file
    cif_file = pdb_file.parent / pdb_file.name.replace('.pdb', '.cif')
    try:
        logging.info(f"Saving as {pdb_file.name.replace('.pdb', '.cif')}")
        with open(cif_file, 'w') as cif_fh:
            ihm.dumper.write(cif_fh, [system])
    except Exception as e:
        logging.error("Error saving mmCIF archive", e)


def load_run_parameters(json_file):
    job_params = json.load(open(json_file))
    return job_params


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("run_directory",
                        help="Location of the uncompressed run, ex: /home/rodrigo/runs/47498-protein-protein")

    args = parser.parse_args()

    logging.basicConfig(level='DEBUG',
                        format='[%(asctime)s] %(levelname)s %(message)s',
                        datefmt='%d/%m/%Y %H:%M:%S')

    logging.info('Converting the cluster*.pdb structures to .cif')

    run_path = pathlib.Path(args.run_directory)
    logging.info(f'Looking for models in {run_path}')

    cluster_pdbs = list(run_path.glob('cluster*.pdb'))
    cluster_pdbs.sort()
    logging.info(f'Found {len(cluster_pdbs)} structures')

    parameter_file = load_run_parameters(run_path / 'job_params.json')

    # this can be none
    logging.info(f"Looking for the tblfile field in {run_path / 'job_params.json'}")
    tblfile = parameter_file["tblfile"]
    if tblfile:
        logging.info('tblfile field found, extracting information')
        restraints = extract_restraints(tblfile)
    else:
        logging.warning('tblfile field not found')
        restraints = False

    for cluster_strct in cluster_pdbs:
        logging.info(f'Converting {cluster_strct}')
        interface_res = get_interface_residues(cluster_strct)
        interface_ranges = get_interface_ranges(interface_res)
        create_ihm_mmcif(cluster_strct, parameter_file, interface_ranges, restraints)

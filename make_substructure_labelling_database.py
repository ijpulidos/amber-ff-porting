import os
import rdkit
from rdkit import Chem
import glob
import json
from collections import defaultdict
import openff
from openff.toolkit.topology import Molecule
import copy
#import MDAnalysis as mda
import mdtraj as md
from simtk.openmm import unit


def remove_charge_and_bond_order_from_guanidinium(rdmol):
    """
    To correct for chemical perception issues with possible resonance states of arginine, 
    remove all charge from the guanidinium group, and set all bond orders to 4. This will
    mark the resonant bonds with a unique "$" character in the SMARTS, which we can later 
    replace. 
    """
    for atom in rdmol.GetAtoms():
        #print(dir(atom))
        if atom.GetAtomicNum() != 6: # element.symbol != "C":
            continue
        nitrogen_neighbors = 0
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 7:
                nitrogen_neighbors += 1
        if nitrogen_neighbors != 3:
            continue
        atom.SetFormalCharge(0)
        for neighbor in atom.GetNeighbors():
            neighbor.SetFormalCharge(0)
        for bond in atom.GetBonds():
        # Set bond order 4, which will produce a "$" character. We later replace this with "~".
            #print(dir(bond))
            bond.SetBondType(Chem.BondType.QUADRUPLE)
            #bond.SetBondType(Chem.BondType.AROMATIC)
            #bond.SetBondType(Chem.BondType.OTHER)
            
def remove_charge_and_bond_order_from_imidazole(rdmol):
    """
    To correct for chemical perception issues with possible resonance states of histidine, 
    remove all charge from the imidazole group, and set all bond orders to 4. This will
    mark the resonant bonds with a unique "$" character in the SMARTS, which we can later 
    replace. 
    """
    #matches = offmol.chemical_environment_matches('[C:1]1~[C:2]~[N:3]~[C:4]~[N:5]1')
    imidazole_substructure = Chem.MolFromSmarts('[C:1]1~[C:2]~[N:3]~[C:4]~[N:5]1')
    matches = rdmol.GetSubstructMatches(imidazole_substructure)
    all_imidazole_atoms = set()
    for match in matches:
        for idx in match:
            all_imidazole_atoms.add(idx)

    for atom in rdmol.GetAtoms():
        if atom.GetIdx() in all_imidazole_atoms:
            atom.SetFormalCharge(0)
        
    for bond in rdmol.GetBonds():
        #print(dir(bond))
        if ((bond.GetBeginAtomIdx() in all_imidazole_atoms) and
            (bond.GetEndAtomIdx() in all_imidazole_atoms)):
            bond.SetBondType(Chem.BondType.QUADRUPLE)

def fix_carboxylate_bond_orders(offmol):
    """Fix problem where leap-produced mol2 files have carboxylates defined with all single bonds"""
    # First, find carbanions
    for atom1 in offmol.atoms:
        if atom1.atomic_number == 6 and atom1.formal_charge.value_in_unit(unit.elementary_charge) == -1:
            # Then, see if they're bound to TWO oxyanions
            oxyanion_seen = False
            for bond in atom1.bonds:
                atom2 = [atom for atom in bond.atoms if not atom == atom1][0]
                if atom2.element.atomic_number == 8 and atom2.formal_charge.value_in_unit(unit.elementary_charge) == -1:
                    # If we find a bond to a SECOND oxyanion, then zero both 
                    # the carbon and the second oxygen's formal charges, and 
                    # set the bond order to 2
                    if oxyanion_seen:
                        atom1._formal_charge = 0 * unit.elementary_charge
                        atom2._formal_charge = 0 * unit.elementary_charge
                        bond._bond_order = 2
                    oxyanion_seen = True


def get_subrdmol(rdmol: rdkit.Chem.rdchem.Mol, indices: list=[],
                 sanitize: bool=False):
    """Create new sub-molecule from selected atom indices
    Parameters
    ----------
    rdmol: rdkit.Chem.rdchem.Mol
        Input molecule
    indices: iterable of ints
        atom indices to include from input molecule, indexed from 0
    sanitize: bool
        whether to sanitize the molecule (recommend: no)
    Returns
    -------
    rdkit.Chem.rdchem.Mol: subset of molecule
    """
    submol = Chem.RWMol(rdmol)
    ix = sorted([at.GetIdx() for at in rdmol.GetAtoms()
                 if at.GetIdx() not in indices])
    for i in ix[::-1]:
        submol.RemoveAtom(int(i))
    if sanitize:
        Chem.SanitizeMol(submol)
        
    for atom in submol.GetAtoms():
        #print(dir(atom))
        atom.SetNoImplicit(True)
        
    remove_charge_and_bond_order_from_imidazole(submol)
    remove_charge_and_bond_order_from_guanidinium(submol)
    return submol


def get_smarts(mol: openff.toolkit.topology.molecule.Molecule,
               indices: list=[], label_indices: list=[]):
    """Get SMARTS of selected atoms in molecule
    Parameters
    ----------
    mol: openff.toolkit.topology.molecule.Molecule
        Input molecule
    indices: iterable of ints
        atom indices to include in SMARTS, indexed from 0
    label_indices: iterable of ints
        atom indices to label, indexed from 0. The atoms
        will be labelled in the order specified. Labels
        begin from 1.
    Returns
    -------
    str: SMARTS string
    """
    rdmol = mol.to_rdkit()
    for i, lix in enumerate(label_indices, 1):
        at = rdmol.GetAtomWithIdx(int(lix))
        at.SetAtomMapNum(i)
    indices = sorted(set(indices) | set(label_indices))
    submol =  get_subrdmol(rdmol, indices)
    #smarts = Chem.MolToSmarts(submol, isomericSmiles=False)
    smarts = Chem.MolToSmarts(submol, isomericSmiles=True)
    smarts = smarts.replace('$', '~')
    return smarts


def get_smiles(mol: openff.toolkit.topology.molecule.Molecule,
               indices: list=[], label_indices: list=[]):
    """Get SMARTS of selected atoms in molecule
    Parameters
    ----------
    mol: openff.toolkit.topology.molecule.Molecule
        Input molecule
    indices: iterable of ints
        atom indices to include in SMARTS, indexed from 0
    label_indices: iterable of ints
        atom indices to label, indexed from 0. The atoms
        will be labelled in the order specified. Labels
        begin from 1.
    Returns
    -------
    str: SMARTS string
    """
    rdmol = mol.to_rdkit()
    for i, lix in enumerate(label_indices, 1):
        at = rdmol.GetAtomWithIdx(int(lix))
        at.SetAtomMapNum(i)
    indices = sorted(set(indices) | set(label_indices))
    submol =  get_subrdmol(rdmol, indices)
    smiles = Chem.MolToSmiles(submol, allHsExplicit=True)
    smiles = smiles.replace('$','~')
    return smiles

def build_residue_library(output_file):
    """Build residue library file from MainChain SMARTS representations."""
    data_dict = defaultdict(list)  # dictionary where to store the residues substructures
    
    for representation in ['MainChain', 'NTerminal', 'CTerminal']:
        single_res_files = glob.glob(f'{representation}/???/???.sdf')
        for filename in single_res_files:
            offmol = Molecule.from_file(filename)
            pdb_name = filename[:-3] + 'pdb'
            md_top = md.load(pdb_name).topology
            # Get the name of the residue from the filename
            base_filename = os.path.basename(filename)
            residue_filename = os.path.splitext(base_filename)[0]
            for residue in md_top.residues:
                atom_names = []
                atom_indices = []
                residue_name = residue.name
                for atom in residue.atoms:
                    atom_indices.append(atom.index)
                    atom_names.append(atom.name)
                smarts = get_smiles(offmol, 
                                    indices=atom_indices, 
                                    label_indices=atom_indices)
                #smarts = get_smarts(offmol, 
                #                    indices=atom_indices, 
                #                    label_indices=atom_indices)
                #print([residue_filename, smarts])
                # Change the residue name to the one in file for unrecognized protonation states
                if residue_filename in ['CYX', 'HID', 'HIE', 'HIP'] and residue_name not in ['ACE', 'NME']:
                    residue_name = residue_filename
                data_dict[residue_name].append(smarts)

    # Remove duplicated smarts
    for residue, smarts_list in data_dict.items():
        data_dict[residue] = list(set(smarts_list))
    # Write dictionary into library text file
    with open(output_file, 'w') as library_file:
        json.dump(data_dict, library_file, indent=4)

build_residue_library('single_res_library.json')
        
    


#single_res_mols = glob.glob('*/???/???.sdf')
#for filename in single_res_mols:
    #print()
    #print(filename)
    #offmol = Molecule.from_file(filename)
    #bond_orderless_offmol = copy.deepcopy(offmol)
    #for bond in bond_orderless_offmol.bonds:
        #bond.bond_order = 1
    #for atom in bond_orderless_offmol.atoms:
        #atom.formal_charge = 0
    #pdb_name = filename[:-3] + 'pdb'
    #md_top = md.load(pdb_name).topology
    #for residue in md_top.residues:
        #atom_names = []
        #atom_indices = []
        #residue_name = residue.name
        #for atom in residue.atoms:
            #atom_indices.append(atom.index)
            #atom_names.append(atom.name)
        #smarts = get_smarts(offmol, 
                            #indices = atom_indices, 
                            #label_indices=atom_indices)
        #smiles = get_smiles(offmol,
                            #indices=atom_indices,
                            #label_indices=atom_indices)
        #bo_less_smarts = get_smarts(bond_orderless_offmol, 
                            #indices = atom_indices, 
                            #label_indices=atom_indices)
        #print()
        #print(residue_name)
        #print(smarts)
        #print(smiles)
        #print(bo_less_smarts)
        #print(atom_names)
        
# How to load PDB?
#    * Use MDTraj, 
#    * OpenMM, 
#    * MDAnalysis, or
#    * RDKit
# How to do substructure matching to PDB?
#    * Pretend it's a cheminformatics mol with all bond orders=1 and all formal charges = 0, then use substructure matching
#    * Make it into a networkX graph and do substructurre matching using that

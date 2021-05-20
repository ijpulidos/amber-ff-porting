import os
import glob
from openff.toolkit.topology import Molecule
from openeye.oechem import OEGraphMol, OEParseSmiles, OECreateCanSmiString


def oe_canonicalize(smiles):
    """Canonicalize smiles using OpenEye."""
    mol = OEGraphMol()
    OEParseSmiles(mol, smiles)
    return OECreateCanSmiString(mol)

sdfs = glob.glob('./*/*/*sdf')
for sdf in sdfs:
    print()
    #print(sdf)
    #print(offmol.to_smiles(mapped=True))
    offmol = Molecule.from_file(sdf)
    smiles = offmol.to_smiles(mapped=True)
    #smiles = offmol.to_smiles(mapped=False, explicit_hydrogens=False)
    smiles_no_nums = "".join([element for element in smiles if not element.isdigit()])
    smiles_clean = smiles_no_nums.replace(':', '')
    print(smiles_clean)
    #print(smiles_clean)
    sdf_split = sdf.split('/')
    sequence = sdf_split[-2]
    sequence = sequence.replace('_', ' ')
    #print(sequence)
    #if 'MainChain' in sdf:
        #print(offmol.to_smiles(), '-->', 'ACE', sequence[:3], sequence[-3:], 'NME')
    #if 'CTerminal' in sdf:
        #print(offmol.to_smiles(), '-->', 'ACE', sequence[:3], sequence[-3:])
    #if 'NTerminal' in sdf:
        #print(offmol.to_smiles(), '-->', sequence[:3], sequence[-3:], 'NME')
    if 'MainChain' in sdf:
        print(sdf_split[1], oe_canonicalize(offmol.to_smiles()), '-->', 'ACE', sequence, 'NME')
    if 'CTerminal' in sdf:
        print(sdf_split[1], oe_canonicalize(offmol.to_smiles()), '-->', 'ACE', sequence)
    if 'NTerminal' in sdf:
        print(sdf_split[1], oe_canonicalize(offmol.to_smiles()), '-->', sequence, 'NME')

from openff.toolkit.topology import Molecule
from utils import fix_carboxylate_bond_orders
import glob

mol2s = glob.glob('./*/*/*.mol2')
for mol2 in mol2s:
    offmol = Molecule.from_file(mol2)
    outfile_name = mol2.replace('.mol2','.sdf')
    fix_carboxylate_bond_orders(offmol)
    offmol.to_file(outfile_name, file_format='sdf')
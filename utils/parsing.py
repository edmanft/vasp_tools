import numpy as np
import matplotlib.pyplot as plt
import ase
from ase.build import fcc111
from ase.visualize.plot import plot_atoms
from ase.constraints import FixAtoms
from ase.io import write
from ase.build import molecule
from ase.io import read
import json
from ase import Atoms
import pubchempy as pcp

def parse_pubchem_json(path):
    """Helper function for parsing the Pubchem Conformer JSON of the molecule to
        an Atoms object from ase """
    with open(path) as f:
        conformer_json = json.load(f)
    elements = np.asarray(conformer_json['PC_Compounds'][0]['atoms']['element'])
    positions = np.zeros([len(elements), 3])
    positions[:, 0] = np.asarray(conformer_json['PC_Compounds'][0]['coords'][0]['conformers'][0]['x'])
    positions[:, 1] = np.asarray(conformer_json['PC_Compounds'][0]['coords'][0]['conformers'][0]['y'])
    positions[:, 2] = np.asarray(conformer_json['PC_Compounds'][0]['coords'][0]['conformers'][0]['z'])
    
    return Atoms(elements, positions=positions)



def parse_from_CID(cid):
    """Helper function for accesing the Pubchem API to obtain the Conformer JSON and parse to
        an Atoms object from ase """
    conformer_json = pcp.get_json(cid, record_type='3d')
    elements = np.asarray(conformer_json['PC_Compounds'][0]['atoms']['element'])
    positions = np.zeros([len(elements), 3])
    positions[:, 0] = np.asarray(conformer_json['PC_Compounds'][0]['coords'][0]['conformers'][0]['x'])
    positions[:, 1] = np.asarray(conformer_json['PC_Compounds'][0]['coords'][0]['conformers'][0]['y'])
    positions[:, 2] = np.asarray(conformer_json['PC_Compounds'][0]['coords'][0]['conformers'][0]['z'])
    
    return Atoms(elements, positions=positions)

def atom_sites(atoms):
    '''Automatically obtain the atom and bond sites in the correct format for the AFM simulation 
        input file'''
    # Obtain atom sites
    atoms_dict = parse_formula(str(atoms.symbols))
    atom_list = list()
    for symbol, count in atoms_dict.items():
        for i in range(1, 1+count):
            atom_list.append(symbol + str(i))
   
    positions = [list(atom.position) for atom in iodo_phenil]
    atom_pos_dict = dict(zip(atom_list, positions))
    
    # Obtain bond sites
    nb_list = build_neighbor_list(atoms)
    connectivity_m = nb_list.get_connectivity_matrix()
    connectivity_m.setdiag(0)
    
    bond_dict = dict()
    for (index_1, index_2), _ in connectivity_m.items():
        bond_name = str(atom_list[index_1])+'-' + str(atom_list[index_2])
        bond_dict[bond_name] = list((iodo_phenil[index_1].position + iodo_phenil[index_1].position)/2)
        
    for symbol, position in atom_pos_dict.items():
        print(symbol, *position)
    for bond_name, position in bond_dict.items():
        print(bond_name, *position)
    return None

class Sites():
    def __init__(self, names, positions, numbers):
        self.names = names
        self.positions = positions
        self.numbers = numbers
        self.colors = jmol_colors[self.numbers]
        
    def describe_sites(self):
        for names, position in zip(self.names, self.positions):
            print(names,'  ' , *position)

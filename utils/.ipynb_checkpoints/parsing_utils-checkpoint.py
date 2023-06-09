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
    """Helper function for parsing the Pubchem JSON of the molecule to
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
    # Download the compound data
    conformer_json = pcp.get_json(cid, record_type='3d')
    elements = np.asarray(conformer_json['PC_Compounds'][0]['atoms']['element'])
    positions = np.zeros([len(elements), 3])
    positions[:, 0] = np.asarray(conformer_json['PC_Compounds'][0]['coords'][0]['conformers'][0]['x'])
    positions[:, 1] = np.asarray(conformer_json['PC_Compounds'][0]['coords'][0]['conformers'][0]['y'])
    positions[:, 2] = np.asarray(conformer_json['PC_Compounds'][0]['coords'][0]['conformers'][0]['z'])
    
    return Atoms(elements, positions=positions)
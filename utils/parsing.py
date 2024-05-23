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

def extract_sites(file_path, exclude_commented=True):
    sites_dict = {}
    with open(file_path, 'r') as file:
        lines = file.readlines()
        
    in_sites_section = False
    for line in lines:
        line = line.strip()
        if line == "&SITES":
            in_sites_section = True
            continue
        if line.startswith("&") and in_sites_section:
            break
        if in_sites_section:
            if line and (not exclude_commented or not line.startswith('#')):
                parts = line.split()
                site_name = parts[0]
                if not site_name.startswith('#'):
                    coordinates = np.array([float(coord) for coord in parts[1:]])
                    sites_dict[site_name] = coordinates
                
    return sites_dict
    
    def parse_dft_to_dataframe(file_path):
    """
    Parse the DFT calculation results from the given file into a pandas DataFrame.
    
    Args:
    - file_path (str): Path to the file containing the DFT calculation results.
    
    Returns:
    - pandas.DataFrame: A DataFrame with 'z' as index and site names as columns.
    """
    # Dictionary to hold our data, with sites as keys and lists of tuples as values
    data = {}

    with open(file_path, 'r') as file:
        for line in file:
            # New site calculation
            if line.startswith("Calculating site..."):
                current_site = line.strip().split("... ")[1]
                if current_site not in data:
                    data[current_site] = []
            elif "eV" in line:
                z, energy = line.strip().split("\t")
                data[current_site].append((float(z), float(energy.split(" ")[0])))
                
    # Convert to DataFrame
    # First, transform data to a format suitable for DataFrame construction
    df_data = {}
    for site, values in data.items():
        for z, energy in values:
            if z not in df_data:
                df_data[z] = {}
            df_data[z][site] = energy
            
    # Create DataFrame
    df = pd.DataFrame.from_dict(df_data, orient='index').sort_index()
    df = df.rename_axis('z').reset_index()
    
    return df

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

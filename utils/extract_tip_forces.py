"""
Python function to extract the forces on the tip after calculating the DFT spectrocopy.
Here we assume that the system is written as sample+tip without posterior sorting, so the last atoms
in the POSCAR are directly the tip. If this is not the case, please change the script to point to the 
specific indexs of the tip in the system.

Example python extract_tip_forces.py --folder_list top_* hollow_* --tip_poscar tip/POSCAR --out2xyz_path $PUBLIC/out2xyz_forces_repeat.pl --forces_csv_path spectroscopy_data/oh_h_clean_dE.csv
"""

import argparse
import os
import shutil
import glob
import subprocess
import numpy as np
import pandas as pd
from ase import Atoms
from ase.io import read


def parse_arguments():
    """Parses arguments from the CLI."""
    parser = argparse.ArgumentParser(description="Extract DFT forces on the tip from spectroscopy")
    
    parser.add_argument('--folder_list', nargs='+', required=True, help='List of folder patterns to include, e.g., top_*, hollow_*')
    parser.add_argument('--tip_poscar', default='tip/POSCAR', type=str, help='Path to the tip POSCAR file.')
    parser.add_argument('--out2xyz_path', required=True, help='Path to the out2xyz_forces_repeat.pl script')
    parser.add_argument('--forces_csv_path', required=True, help='Path of final csv')
    parser.add_argument('--debug', action='store_true', help='True if we just want to print the folders')
    
    args = parser.parse_args()

    args.tip_poscar = os.path.abspath(args.tip_poscar)
    args.forces_csv_path = os.path.abspath(args.forces_csv_path)
    
    return args


def _is_float(path: str) -> bool:
    """Check if the path basename is a float indicating the folder represents a height."""
    try:
        float(path)
        return True
    except ValueError:
        return False


def _filter_heights(paths: list[str]) -> list[str]:
    """Filter for directories that represent a physical height."""
    filtered_paths = []
    for path in paths:
        if os.path.isdir(path) and _is_float(os.path.basename(path)):
            filtered_paths.append(path)
    return filtered_paths


def read_last_xyz(xyz_path: str):
    """Read the atoms and forces from a .xyz file."""
    with open(xyz_path, 'r') as file:
        lines = file.readlines()

    natoms = int(lines[0].strip())
    symbols = []
    positions = []
    forces = []
    
    for i in range(2, 2 + natoms):
        line = lines[i].strip().split()
        symbols.append(line[0])
        positions.append([float(x) for x in line[1:4]])
        forces.append([float(x) for x in line[4:7]])

    atoms = Atoms(symbols=symbols, positions=positions)
    forces = np.array(forces)
    
    return atoms, forces


def _extract_tip_force(folder: str, args, tip_len: int):
    """Extract the force on the tip from spectroscopy data."""
    heights = sorted(glob.glob(os.path.join(folder, '*')))
    heights = _filter_heights(heights)
    
    z_list = []
    force_list = []
    
    for height in heights:
        z_list.append(os.path.basename(height))
        try:
            shutil.copy(args.out2xyz_path, os.path.join(height, 'out2xyz_forces_repeat.pl'))
            os.chdir(height)
            subprocess.run('./out2xyz_forces_repeat.pl 1, 1, 1', shell=True)
            
            last_xyz = os.path.join(height, 'last.xyz')
            if os.path.exists(last_xyz):
                system, forces = read_last_xyz(last_xyz)
                tip_sum_forces = forces[-tip_len:, 2].sum()
                force_list.append(tip_sum_forces)
        except Exception as e:
            print(f"Error processing {height}: {e}")
    
    return z_list, force_list

def store_forces_dict(forces_dict, args):
    """
    Checks if all z values are the same across sites, then converts the dictionary containing 
    force data into a DataFrame.
    
    Parameters:
    forces_dict (dict): Dictionary with force data, structured with site keys and 'z' and 'force' subkeys.
    args (Namespace): Arguments from argparse, including the path for saving the CSV file.
    """
    z_values = [data['z'] for site, data in forces_dict.items()]
    
    if not all(z == z_values[0] for z in z_values):
        raise ValueError("Not all sites have the same z values")
    
    df = pd.DataFrame({
        'z': z_values[0]
    })

    for site, data in forces_dict.items():
        df[site] = data['force']

    df.reset_index(drop=True, inplace=True)

    csv_dir = os.path.dirname(args.forces_csv_path)
    if not os.path.exists(csv_dir):
        os.makedirs(csv_dir)
        print(f"Created directory {csv_dir}")

    df.to_csv(args.forces_csv_path, index=False)
    print(f"Data saved to {args.forces_csv_path}")

    

def main():
    args = parse_arguments()
    
    if args.debug:
        for pattern in args.folder_list:
            folders = glob.glob(pattern)
            for folder in folders:
                print(folder)
                heights = _filter_heights(glob.glob(os.path.join(os.path.abspath(folder), '*')))
                for height in heights:
                    print(height)
    else:
        work_dir = os.getcwd()
        tip = read(args.tip_poscar)
        tip_len = len(tip)
        sites_force_dict = {}

        for pattern in args.folder_list:
            folders = glob.glob(pattern)
            for folder in folders:
                print(f"Processing {folder}:")
                z_list, force_list = _extract_tip_force(os.path.abspath(folder), args, tip_len)
                sites_force_dict[os.path.basename(folder)] = {'z': z_list, 'force': force_list}
                os.chdir(work_dir)
                print(f"Finished processing {folder}.")

        # Convert sites_force_dict to pandas Dataframe to store as csv.
        # If the heights of the spectroscopy don't match among sites, this step won't work.
        # In that case, simply store as json directly.
        store_forces_dict(sites_force_dict, args)



if __name__ == "__main__":
    main()

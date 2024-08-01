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
from ase import Atoms
from ase.io import read 

def parse_arguments():
    """Parses arguments from the CLI."""
    parser = argparse.ArgumentParser(description="Extract DFT forces on the tip from spectroscopy")
    
    # Adding arguments
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
    if path is None:
        return False
    try:
        float(path)
        return True
    except ValueError:
        return False

def _filter_heights(paths: list[str]) -> list[str]:
    filtered_paths = []
    for path in paths:
        if os.path.isdir(path) and _is_float(os.path.basename(path)):
            filtered_paths.append(path)
    return filtered_paths

def read_last_xyz(xyz_path):

    with open(xyz_path, 'r') as file:
        lines = file.readlines()

    natoms = int(lines[0].strip())
    comment = lines[1].strip()
    
    symbols = []
    positions = []
    forces = []
    
    for i in range(2, 2 + natoms):
        line = lines[i].strip().split()
        symbols.append(line[0])
        position = [float(x) for x in line[1:4]]
        positions.append(position)
        force = [float(x) for x in line[4:7]]
        forces.append(force)

    atoms = Atoms(symbols=symbols, positions=positions)
    forces = np.array(forces)
    return atoms, forces
    


def _extract_tip_force(folder, args, tip_len):
    """Generates the last.xyz on each tip position and gets the forces"""
    heights = sorted(glob.glob(os.path.join(folder, '*')))
    heights = _filter_heights(heights)
    
    z_list = []
    force_list = []
    for height in heights:
        z_list.append(os.path.basename(height))
        try:
            
            shutil.copy(args.out2xyz_path, os.path.join(height, 'out2xyz_forces_repeat.pl'))
            
            # Run the script using subprocess
            os.chdir(height)
            subprocess.run('./out2xyz_forces_repeat.pl 1, 1, 1', shell=True)
            
            # Assume the script generates a file named 'last.xyz' in the height directory
            last_xyz = os.path.join(height, 'last.xyz')
            if os.path.exists(last_xyz):
                # read forces of tip
                system, forces = read_last_xyz(last_xyz)
                tip_sum_forces = forces[-tip_len:, 2].sum()
                force_list.append(tip_sum_forces)

        except Exception as e:
            print(f"Error processing {height}: {e}")
    
    return z_list, force_list
    
    

def main():
    args = parse_arguments()
    if args.debug:
        for pattern in args.folder_list:
            folders = glob.glob(pattern)
            for folder in folders:
                print(folder)
                folder = os.path.abspath(folder)
                heights = glob.glob(os.path.join(folder, '*'))
                heights = _filter_heights(heights)
                for height in heights:
                    print(height)
    else:
        work_dir = os.getcwd()
        
        tip = read(args.tip_poscar)
        tip_len = len(tip)

        sites_force_dict = dict()
        # Iterate over all folders matching the patterns in folder_list
        for pattern in args.folder_list:
            folders = glob.glob(pattern)
            print(folders)
            for folder in folders:
                print(f"Processing {folder}:")
                folder = os.path.abspath(folder)
                z_list, force_list = _extract_tip_force(folder, args, tip_len)

                sites_force_dict[os.path.basename(folder)] = {'z': z_list, 'force': force_list}
                os.chdir(work_dir)
                print(f"Finished processing {folder}.")
        print(sites_force_dict)




if __name__ == "__main__":
    main()
    



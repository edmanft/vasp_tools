import re
import argparse
import os
import sys
import shutil
import json
import subprocess
from pathlib import Path
import numpy as np
import ase
from ase.io import write, read



def parse_arguments():
    """Parses argument from the CLI. All paths are converted to absolute paths."""
    parser = argparse.ArgumentParser(description="Create DFT spectroscopies of a tip moving in the z direction on top of a sample.")
    
    # Adding default values to the arguments
    parser.add_argument('--label', required=True, type=str, help='Name of the folder for multiple calculations.')
    parser.add_argument('--incar', default='INCAR', type=str, help='Path to the INCAR file.')
    parser.add_argument('--kpoints', default='KPOINTS', type=str, help='Path to the KPOINTS file.')
    parser.add_argument('--potcar', default='POTCAR', type=str, help='Path to the POTCAR file.')
    parser.add_argument('--sample_poscar', default='sample/POSCAR', type=str, help='Path to the sample POSCAR file.')
    parser.add_argument('--tip_poscar', default='tip/POSCAR', type=str, help='Path to the tip POSCAR file.')
    parser.add_argument('--x', required=True, type=float, help='X position of the sample.')
    parser.add_argument('--y', required=True, type=float, help='Y position of the sample.')
    parser.add_argument('--zref', required=True, type=float, help='Reference height of the sample')
    parser.add_argument('--z0', required=True, type=float, help='Starting height.')
    parser.add_argument('--zf', required=True, type=float, help='Finishing height.')
    parser.add_argument('--dz', required=True, type=float, help='Step height.')
    parser.add_argument('--vasp_instruct', required=True, type=str, help='Path to file with commands to execute VASP and postprocessing.')
    parser.add_argument('--restart', action='store_true', help='Restart a calculation from an intermediate position.')
    args = parser.parse_args()
    
    # Check that dz is positive
    if args.dz <= 0:
        print("Error: dz must be a positive value.")
        sys.exit(1)
    
    # Get the absolute paths
    args.label = os.path.abspath(args.label)
    args.incar = os.path.abspath(args.incar)
    args.kpoints = os.path.abspath(args.kpoints)
    args.potcar = os.path.abspath(args.potcar)
    args.sample_poscar = os.path.abspath(args.sample_poscar)
    args.tip_poscar = os.path.abspath(args.tip_poscar)
    args.vasp_instruct = os.path.abspath(args.vasp_instruct)
    
    return args

def create_calculation_folders(args):
    """Create the folders and copy the INCAR, KPOINTS and POTCAR for the calculation"""
    
    if not os.path.exists(args.label):
        os.makedirs(args.label)
    
    args_dict = vars(args)
    if args.restart:
        with open(os.path.join(args.label, 'metadata_restart.json'), 'w') as f:
            json.dump(args_dict, f, indent=4)
    
    else:
        with open(os.path.join(args.label, 'metadata.json'), 'w') as f:
            json.dump(args_dict, f, indent=4)
    
    # Create all the calculation folders and copy the INCAR, POTCAR, KPOINTS and execute_relaxation_sh 

    z0, zf, dz = args.z0, args.zf, args.dz 

    # Compute the heights of the spectroscopy.
    # The tip can approach to the sample (dz<0) or get away from it (dz>0)
    npoints = int(np.abs(zf-z0)/dz) +1
    sign_dz = np.sign(zf-z0) 
    
    zrange = np.arange(0, npoints)*dz*sign_dz + z0

    for z in zrange:
        calc_dir = os.path.join(args.label,f"{z:.2f}")  
        if not os.path.exists(calc_dir):
            os.makedirs(calc_dir)
    
        try:
            shutil.copy(args.incar, os.path.join(calc_dir, "INCAR"))
            shutil.copy(args.kpoints, os.path.join(calc_dir, "KPOINTS"))
            shutil.copy(args.potcar, os.path.join(calc_dir, "POTCAR"))
            shutil.copy(args.vasp_instruct, os.path.join(calc_dir, "execute_relaxation_sh"))
            
            
            
        except Exception as e:
            raise Exception(f"Failed to copy files to {calc_dir}: {e}")

def prepare_sample_tip(args):
    sample = read(args.sample_poscar) 
    tip = read(args.tip_poscar)
    if sample.constraints == []:
        print(f"Couldn't find any constraints in sample: {args.sample_poscar}")
        print("Defaulting to static sample")
        sample_constraint = ase.constraints.FixAtoms([atom.index for atom in sample])
        sample.set_constraint(sample_constraint)        

    if tip.constraints == []:
        print(f"Couldn't find any constraints in tip: {args.tip_poscar}")
        print("Defaulting to static tip")
        tip_constraint = ase.constraints.FixAtoms([atom.index for atom in tip])
        tip.set_constraint(tip_constraint)
    tip.cell = sample.cell
    
    return sample, tip
def remove_slashes_from_files(directory):
    # List of files to check and modify
    files_to_check = ['CONTCAR', 'CHGCAR', 'LOCPOT']
    
    # Check and modify each file directly in the directory
    for filename in files_to_check:
        filepath = Path(directory) / filename
        if filepath.exists():
            subprocess.run(f"sed -i '6s|/||g' {filepath}", shell=True)

def check_calculation_success(calc_dir):
    """
    Checks if the VASP calculation in the specified directory was successful by searching for a specific
    pattern in the OUTCAR file using regular expressions.
    """
    success_pattern = re.compile(r"aborting loop because EDIFF is reached")
    try:
        with open(f"{calc_dir}/OUTCAR", "r") as file:
            content = file.read()
            if re.search(success_pattern, content):
                return True
    except FileNotFoundError:
        print(f"OUTCAR file not found in {calc_dir}")
    except Exception as e:
        print(f"Error reading from {calc_dir}/OUTCAR: {e}")
    return False


def main():
    args = parse_arguments()
    x_pos, y_pos = args.x, args.y
    zref, z0, zf, dz = args.zref, args.z0, args.zf, args.dz 

    
    

    # Compute the heights of the spectroscopy.
    # The tip can approach to the sample (dz<0) or get away from it (dz>0)
    npoints = int(np.abs(zf-z0)/dz) +1
    sign_dz = np.sign(zf-z0) 
    zrange = np.arange(0, npoints)*dz*sign_dz + z0
    
    create_calculation_folders(args)

    # Read the sample and tip POSCAR and set constraints if there aren't any.
    sample, tip = prepare_sample_tip(args)

    
    disp = np.array([x_pos, y_pos, z0 + zref])

    tip_disp = tip.copy()
    tip_disp.positions += disp
    
    system = sample + tip_disp
    
    # Add the constraints for the first configuration of the system
    sample_constraint = sample.constraints[0].get_indices()
    tip_constraint = tip.constraints[0].get_indices() + len(sample) 
    system_constraint = list(sample_constraint) + list(tip_constraint)
    system_constraint = ase.constraints.FixAtoms(system_constraint)
    system.set_constraint(system_constraint)

    tip_indexs = len(sample) + np.arange(len(tip))


    

    calc_dir, prev_calc_dir = None, None
    if args.restart: 
        prev_calc_dir = os.path.join(args.label,f"{z0-dz*sign_dz:.2f}")
        print(f"Restarting calculation from {prev_calc_dir}")

    for z in zrange:
        calc_dir = os.path.join(args.label,f"{z:.2f}")

        if prev_calc_dir is None: ## First calculation
            write(os.path.join(calc_dir, "POSCAR"), system, direct=False)

            
            os.chdir(calc_dir)
            os.system('./execute_relaxation_sh' )
            
                
            remove_slashes_from_files(calc_dir) 

            
        else:  # Posterior calculation, we have CHGCAR and WAVECAR
            try:
                for file in ["WAVECAR", "CHGCAR"]:
                    shutil.move(os.path.join(prev_calc_dir, file), os.path.join(calc_dir, file))
            except Exception as e:
                print(f"Warning: Failed to move WAVECAR and CHGCAR to {calc_dir}: {e}")
                print("""We recommend initializing the calculation with the 
                        WAVECAR and CHGCAR from the previous step""")
            
            
            system = read(os.path.join(prev_calc_dir, "CONTCAR"))
            # God knows why system[tip_indexs].positions += np.array([0, 0, dz*sign_dz]) doesn't work
            # Displace the atoms one at a time
            for index in tip_indexs: 
                system[index].position += np.array([0, 0, dz*sign_dz])
            
            
            write(os.path.join(calc_dir, "POSCAR"), system, direct=False)

            retry_count = 0
            max_retries = 3 # Maximum number of calculations we try

            while retry_count < max_retries:
                os.chdir(calc_dir)
                os.system('./execute_relaxation_sh')

                if check_calculation_success(calc_dir):
                    print(f"Calculation at {calc_dir} successful.")
                    break
                else:
                    print(f"Calculation at {calc_dir} failed.")
                    retry_count += 1

            if retry_count == max_retries:
                print("Maximum retries reached. Exiting job.")
                sys.exit()
            
            remove_slashes_from_files(calc_dir)   
              

        prev_calc_dir = calc_dir




if __name__ == "__main__":
    main()

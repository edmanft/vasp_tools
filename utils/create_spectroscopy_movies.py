"""Example python create_spectroscopy_movies.py --folder_list top_* hollow_* --out2xyz_path $PUBLIC/out2xyz_forces_repeat.pl --cell_repetitions 2 2 1 --spec_movie_path spectroscopy_movies"""

import argparse
import os
import shutil
import glob
import subprocess

def parse_arguments():
    """Parses arguments from the CLI."""
    parser = argparse.ArgumentParser(description="Create spectroscopy movies based on the DFT calculations.")
    
    # Adding arguments
    parser.add_argument('--folder_list', nargs='+', required=True, help='List of folder patterns to include, e.g., top_*, hollow_*')
    parser.add_argument('--ascending', action='store_true', help='True if the tip distancing from sample (dz>0)')
    parser.add_argument('--out2xyz_path', required=True, help='Path to the out2xyz_forces_repeat.pl script')
    parser.add_argument('--cell_repetitions', nargs=3, type=int, default=[2, 2, 1], help='Repetitions for the unit cell in x, y, z directions')
    parser.add_argument('--spec_movie_path', required=True, help='Path where the movie files will be stored')
    parser.add_argument('--debug', action='store_true', help='True if we just want to print the folders')
    args = parser.parse_args()

    args.spec_movie_path = os.path.abspath(args.spec_movie_path)
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


def _generate_movie(folder, args):
    """Generates the spectroscopy movie for a specific site (folder)."""
    heights = sorted(glob.glob(os.path.join(folder, '*')), reverse=not args.ascending)
    heights = _filter_heights(heights)
    
    movie_frames = []
    for height in heights:
        try:
            
            shutil.copy(args.out2xyz_path, os.path.join(height, 'out2xyz_forces_repeat.pl'))
            
            # Run the script using subprocess
            cell_str = ' '.join(str(x) for x in args.cell_repetitions)
            submit_command = "./out2xyz_forces_repeat.pl " + cell_str
            os.chdir(height)
            subprocess.run(submit_command, shell=True)
            
            # Assume the script generates a file named 'last.xyz' in the height directory
            last_xyz = os.path.join(height, 'last.xyz')
            if os.path.exists(last_xyz):
                movie_frames.append(last_xyz)
        except Exception as e:
            print(f"Error processing {height}: {e}")
    
    # Concatenate the .xyz files into a single movie file
    movie_file = os.path.join(args.spec_movie_path, f"movie_{os.path.basename(folder)}.xyz")
    with open(movie_file, 'w') as outfile:
        for frame in movie_frames:
            with open(frame, 'r') as infile:
                outfile.write(infile.read())

def main():
    args = parse_arguments()
    if args.debug:
        for pattern in args.folder_list:
            folders = glob.glob(pattern)
            for folder in folders:
                print(folder)
                folder = os.path.abspath(folder)
                heights = sorted(glob.glob(os.path.join(folder, '*')), reverse=not args.ascending)
                heights = _filter_heights(heights)
                for height in heights:
                    print(height)
    else:
        work_dir = os.getcwd()
        # Create the directory for the movie files if it does not exist
        os.makedirs(args.spec_movie_path, exist_ok=True)
        
        # Iterate over all folders matching the patterns in folder_list
        for pattern in args.folder_list:
            folders = glob.glob(pattern)
            print(folders)
            for folder in folders:
                print(f"Processing {folder}:")
                folder = os.path.abspath(folder)
                _generate_movie(folder, args)
                os.chdir(work_dir)
                print(f"Finished processing {folder}.")




if __name__ == "__main__":
    main()
    

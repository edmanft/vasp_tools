def save_povray(trilayer, category, base_path, file_name, rotations, canvas_width=2000, show_unit_cell=False, show_bonds=False):
    """
    Save .pov, .ini, and rendered .png images for the given `trilayer` structure using POV-Ray.
    
    Args:
        trilayer (Atoms): The ASE Atoms object to render.
        category (str): The category ("Clean" or "Ov").
        base_path (str): Base directory where files are saved.
        file_name (str): Name of the input structure file.
        rotations (dict): Dictionary of rotation labels and their settings.
        canvas_width (int): Resolution width for rendering (default: 2000px).
    """
    category_path = os.path.join(base_path, "3x3", category, file_name)
    os.makedirs(category_path, exist_ok=True)
    
    for view_name, rotation in rotations.items():
        pov_file = os.path.join(category_path, f"{view_name}.pov")
        ini_file = f'{view_name}.ini'
        png_file = ini_file.replace('.pov', '.png')

        if show_bonds:
            bondatoms = get_bondpairs(trilayer)
            povray_settings=dict(canvas_width=canvas_width, bondatoms = bondatoms)
        else:
            povray_settings=dict(canvas_width=canvas_width)
        
        # Write the .pov file
        write(
            pov_file,
            trilayer,
            format='pov',
            rotation=rotation,
            povray_settings=povray_settings, 
            show_unit_cell=show_unit_cell, 
        )
        
        # Run POV-Ray to render the PNG
        os.chdir(category_path)
        print(f'category_path: {category_path}')
        print(f'pov_file: {ini_file}')
        os.system(f"povray {ini_file} +O{png_file} +A")

def process_files(file_paths, sketches_path, canvas_width=2000, show_unit_cell=False, show_bonds=False):
    """
    Process a list of file paths and generate rendered images with POV-Ray.
    
    Args:
        file_paths (list): List of file paths to process.
        sketches_path (str): Path to save rendered images and files.
        canvas_width (int): Resolution width for rendering (default: 2000px).
    """
    rotations = {
        "top": "",
        "front1": "-90x",
        "front2": "-90x,-45y"
    }
    
    for path in file_paths:
        sample = read(path)
        trilayer = Atoms([atom for atom in sample if atom.position[2] > 9])
        trilayer.cell = sample.cell
        trilayer.pbc = sample.pbc
        trilayer.positions[:, 2] -= trilayer.positions[:, 2].min()
        
        # Modify atoms based on position
        for atom in trilayer:
            if atom.position[2] < 0.4:
                atom.symbol = 'S'
        
        # Determine category
        if 'clean' in path.lower():
            category = "Clean"
        elif 'ov' in path.lower():
            category = "Ov"
            trilayer[0].symbol = 'Cl'
            trilayer[1].symbol = 'Cl'
        else:
            continue  # Skip files that don't match Clean or Ov
        
        file_name = os.path.basename(path)
        save_povray(trilayer, category, sketches_path, file_name, rotations, canvas_width, show_unit_cell=show_unit_cell, show_bonds=show_bonds)

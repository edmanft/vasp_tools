def atoms_description(atoms, print_flag=False):
    '''Describes both atom types, coordinates and cell.
    If print_flag = True, plots the sample from different perspectives.'''
    print('*Printing atomic positions:*')
    for atom in atoms:
        print(f'{atom}')
    print('**----------------- \nCell:**')
    for lat_vec in atoms.cell:
        print(lat_vec)
    if print_flag:
        fig, ax = plt.subplots(1,3)
        plot_atoms(atoms, ax[0], radii=0.3, rotation=('90x,45y,0z'))
        plot_atoms(atoms, ax[1], radii=0.3, rotation=('45x,45y,0z'))
        plot_atoms(atoms, ax[2], radii=0.3)
        plt.show()
    return None

def center_atoms(atoms):
    '''Translates the geometrical center of the molecule to the center of the cell'''
    cell = atoms.cell
    cell_center = np.sum(cell, axis = 0)/2
    
    geom_center = np.mean(atoms.positions, axis = 0)
    
    atoms.positions = atoms.positions - geom_center + cell_center
    return atoms


def parse_formula(formula):
    # Regular expression to match element symbols and counts
    pattern = r'([A-Z][a-z]*)(\d*)'

    # Use re.findall to get all matches in the formula
    matches = re.findall(pattern, formula)

    # Create a dictionary to store the element symbols and counts
    counts = {}
    for symbol, count in matches:
        if count == '':
            count = 1
        else:
            count = int(count)
        if symbol in counts:
            counts[symbol] += count
        else:
            counts[symbol] = count

    return counts

def plot_with_positions(atoms, position_dict, color_dict, offset=None):
    # Plot the atoms using ASE's plot_atoms function
    fig, ax = plt.subplots()
    plot_atoms(atoms, ax=ax, rotation=('0x,0y,0z'))
    if offset is None:
        offset = [1.7, 0] # the bottom left corner of the cell is not drawn in (0,0)
    
    # Iterate through the position_dict to draw crosses and legends
    for label, position in position_dict.items():
        # Draw a cross at the specified position
        color = color_dict.get(label, 'black')  # Default color is black if not specified
        ax.scatter(position[0]+offset[0], position[1]+offset[1], color=color, label=label, s=100, marker='x', linewidths=2)
    
    # Add a legend outside of the plot
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    
    plt.show()

    

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
def plot_with_positions(atoms, position_dict, color_dict):
    # Plot the atoms using ASE's plot_atoms function
    fig, ax = plt.subplots()
    plot_atoms(atoms, ax=ax, rotation=('0x,0y,0z'))
    offset = [1.7, 0] # the bottom left corner of the cell is not drawn in (0,0)
    
    # Iterate through the position_dict to draw crosses and legends
    for label, position in position_dict.items():
        # Draw a cross at the specified position
        color = color_dict.get(label, 'black')  # Default color is black if not specified
        ax.scatter(position[0]+offset[0], position[1]+offset[0], color=color, label=label, s=100, marker='x', linewidths=2)
    
    # Add a legend outside of the plot
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    
    plt.show()
    

def center_atoms(atoms, final_height=None):
    '''Translates the geometrical center of the molecule to the center of the cell.'
    final_height: final height of the plane containing the molecule (mean z position).'''
    cell = atoms.cell
    cell_center = np.sum(cell, axis = 0)/2
    
    geom_center = np.mean(atoms.positions, axis = 0)
    
    atoms.positions = atoms.positions - geom_center + cell_center
    
    if final_height is not None:
        min_height = np.min(atoms.positions[:,2])
        atoms.positions[:,2] = atoms.positions[:,2] - min_height + final_height
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

def plot_from_outcar(path, ylim_E=None, ylim_F=None):

    outcar_full = read(path, format='vasp-out', index=':')

    energies = [step.get_total_energy() for step in outcar_full]
    forces = [(step.get_forces()**2).sum() for step in outcar_full]

    fig, ax = plt.subplots(2,1)
    ax[0].plot(np.arange(len(energies))+1, energies, 'g-',label='Energy $(eV)$')
    ax[1].plot(np.arange(len(forces))+1, forces, 'b-', label='Forces $(eV/\AA³)$')
    if ylim_E is not None:
        ax[0].set_ylim(ylim_E)
    if ylim_F is not None:
        ax[1].set_ylim(ylim_F)
    fig.legend()

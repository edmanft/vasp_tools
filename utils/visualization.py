def atoms_description(atoms, print_flag=False):
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

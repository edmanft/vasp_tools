def parse_ddec6_bond_orders(file_path, min_bond_order=0.9):
    bond_orders = {}
    current_atom = None
    
    with open(file_path, 'r') as file:
        lines = file.readlines()
        
        for line in lines:
            line = line.strip()
            if not line:
                continue

            # Identify the atom number
            atom_match = re.match(r'Printing BOs for ATOM #\s+(\d+)', line)
            if atom_match:
                current_atom = int(atom_match.group(1))
                bond_orders[current_atom] = []
                continue
            
            # Parse bond orders
            if current_atom is not None:
                bond_order_match = re.search(r'atom number\s+(\d+).*bond order\s+=\s+([\d.]+)', line, re.IGNORECASE)
                if bond_order_match:
                    bonded_atom = int(bond_order_match.group(1))
                    bond_order = float(bond_order_match.group(2))
                    if bond_order >= min_bond_order:
                        bond_orders[current_atom].append((bonded_atom, bond_order))
    bond_order_dict = {}
    
    for atom, bonds in bond_orders.items():
        for bonded_atom, bond_order in bonds:
            bond = tuple(sorted((atom, bonded_atom)))
            bond_order_dict[bond] = bond_order
    
    return bond_order_dict

def plot_molecular_bond_order(molecule, bond_orders, title, 
                              plot_size=(10, 10), node_size=30, text_size=12, 
                              cmap_name='Blues', zoom_len=2, x0=[2, 2], xf=[15, 15], 
                              bond_order_precision=1):
    pos_dict = {atom.index+1:atom.position[:2] for atom in molecule if atom.symbol == 'C'}
    lx, ly = molecule.cell.diagonal()[:2]

    # Create a graph
    G = nx.Graph()

    # Add nodes with positions
    for atom, coords in pos_dict.items():
        G.add_node(atom, pos=coords)

    # Add edges with bond orders as weights
    for (atom, bonded_atom), bond_order in bond_orders.items():
        if atom in pos_dict and bonded_atom in pos_dict:
            G.add_edge(atom, bonded_atom, weight=bond_order)

    # Get positions for nodes
    pos = nx.get_node_attributes(G, 'pos')

    # Normalize bond order values for the color map
    if bond_orders:
        norm = mcolors.Normalize(vmin=min(bond_orders.values()), vmax=max(bond_orders.values()))
        cmap = cm.get_cmap(cmap_name)

        # Get edge colors based on bond orders
        edge_colors = [cmap(norm(G[u][v]['weight'])) for u, v in G.edges]

        # Set up the plot
        plt.figure(figsize=plot_size)
        ax = plt.gca()

        # Draw nodes
        nx.draw_networkx_nodes(G, pos, node_size=node_size, node_color='black')

        # Draw edges with colors
        nx.draw_networkx_edges(G, pos, edgelist=G.edges, edge_color=edge_colors, width=2)

        # Draw edge labels (bond orders)
        edge_labels = {(u, v): f"{d['weight']:.{bond_order_precision}f}" for u, v, d in G.edges(data=True)}
        nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_color='green', font_size=text_size)

        # Create a colorbar
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = plt.colorbar(sm)
        cbar.set_label('Bond Order')

        # Set the plot title and axis limits
        plt.title(title)
        ax.set_xlim([x0[0], xf[0]])
        ax.set_ylim([x0[1], xf[1]])
        ax.set_aspect('equal')

    # Show the plot
    plt.show()

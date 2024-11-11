#!/usr/bin/env python
# coding: utf-8

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import argparse
from ase.io import read
from ase.visualize.plot import plot_atoms

# Add the path to your custom modules (adjust as necessary)
sys.path.append('/home/eventura/codes/DBSPM')
from pydbspm.grid import Grid
from pydbspm.SPMGrid import SPMGrid



# Function to parse command line arguments
def parse_arguments():
    parser = argparse.ArgumentParser(description='Script for analyzing AFM data of a given molecule')
    parser.add_argument('-m', '--molecule', type=str, required=True, help='Name of the molecule')
    parser.add_argument('-hl', '--height_list', type=float, nargs='+', required=True,
                        help='List of heights (in Angstroms) at which to analyze the data')
    parser.add_argument('-zref', '--zref', type=float, default=3.0, help='Reference height (default: 3.0 Å)')
    parser.add_argument('-A', '--amplitude', type=float, default=0.6, help='Amplitude of the oscillation of the tip (default: 0.6 Å)')
    parser.add_argument('-k0', '--k0', type=float, default=1800.0, help='Spring constant of the cantilever (default: 1800 N/m)')
    parser.add_argument('-f0', '--f0', type=float, default=30300.0, help='Resonance frequency of the cantilever (default: 30300 Hz)')
    return parser.parse_args()

# Main function
def main():
    # Parse arguments
    args = parse_arguments()
    molec_name = args.molecule
    height_list = args.height_list

    # Function to calculate the index along the z-axis
    def nz(z, zref=args.zref, dz=None):
        return int(np.rint((z + zref) / dz))

    # Create the output directory if it doesn't exist
    output_dir = "pp_dbspm"
    os.makedirs(output_dir, exist_ok=True)

    # Read the sample data
    sample = read('sample/POSCAR')
    plot_atoms(sample)
    plt.savefig(os.path.join(output_dir, f'{molec_name}_sample_atoms.png'), bbox_inches='tight')
    plt.close()

    # Load sample data from .npz file
    samp_data = np.load('sample/sample.npz')
    print(samp_data.files)

    # Create a Grid object for the sample data
    samp_grid = Grid(samp_data['shape'], span=samp_data['cell'], origin=samp_data['origin'], zref=3.0)

    # Set densities for the grid
    samp_grid.set_density('rhoS', samp_data['rhoS'], safe=True)
    samp_grid.set_density('potS', samp_data['potS'], safe=True)

    # Calculate electric field and charge density
    dz = samp_grid.dr[2]  # Units: Angstroms

    # Determine subplot configuration for charge density and E_z (always 2 columns)
    n_cols_cd_ez = 2
    n_rows_cd_ez = (len(height_list) + n_cols_cd_ez - 1) // n_cols_cd_ez

    # Plot charge density at different heights and save the figures
    fig, axs = plt.subplots(n_rows_cd_ez, n_cols_cd_ez, figsize=(8, 3 * n_rows_cd_ez))
    axs = np.array(axs).reshape(-1)

    for i, height in enumerate(height_list):
        z_idx = nz(height, dz=dz)
        z_data = samp_data['rhoS'][..., z_idx] * 1000  # Convert to e/nm^3
        ax = axs[i]
        cont = ax.contourf(samp_grid.x, samp_grid.y, z_data, 64, cmap='hot')
        cbar = plt.colorbar(cont, ax=ax)
        cbar.set_label('Charge density $(e/nm^3)$')
        ax.set_aspect('equal')
        ax.axis('off')
        ax.set_title(rf'{molec_name} charge density at $z_{{ts}}={height:.2f}\ \AA$')

    # Hide any unused subplots
    for j in range(len(height_list), len(axs)):
        axs[j].axis('off')

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'{molec_name}_charge_density.png'), bbox_inches='tight')
    plt.close()

    # Plot electric field (E_z) at different heights and save the figures
    potS = samp_data['potS']
    fig, axs = plt.subplots(n_rows_cd_ez, n_cols_cd_ez, figsize=(8, 3 * n_rows_cd_ez))
    axs = np.array(axs).reshape(-1)

    for i, height in enumerate(height_list):
        z_idx = nz(height, dz=dz)
        Ez = -np.gradient(potS, dz, axis=2)  # Calculate E_z
        max_scale = np.max(np.abs(Ez[..., z_idx]))  # Scale based on max value at this height
        ax = axs[i]
        cont = ax.contourf(
            samp_grid.x, samp_grid.y,
            Ez[..., z_idx],
            levels=64,
            vmin=-max_scale,
            vmax=max_scale,
            cmap='RdBu_r'
        )
        cbar = plt.colorbar(cont, ax=ax)
        cbar.set_label(rf'$E_z$ at $z={height:.2f}\ \AA$ (V/\AA)')
        ax.set_title(rf'{molec_name} $E_z$ at $z_{{ts}}={height:.2f}\ \AA$')
        ax.set(aspect='equal')
        ax.axis('off')

    # Hide any unused subplots
    for j in range(len(height_list), len(axs)):
        axs[j].axis('off')

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'{molec_name}_Ez.png'), bbox_inches='tight')
    plt.close()

    # Load interaction grids and compute gradients
    sr = SPMGrid(f'{molec_name}_sr_a1.08.npz', 'sr').get_gradient()
    es = SPMGrid(f'{molec_name}_es.npz', 'es').get_gradient()
    vdw = SPMGrid(f'{molec_name}_vdw.npz', 'vdw').get_gradient()
    relax = SPMGrid(f'relax_{molec_name}_k0.2000_a1.08_V42.91.npz', 'E').get_gradient()

    # Multiply short-range interaction by its potential
    sr *= 42.91

    # Create the static image by summing the interactions
    static = sr + es + vdw

    # Determine subplot configuration for interaction plots (5 columns)
    n_cols_interactions = 5
    n_rows_interactions = len(height_list)

    # Plot individual interactions and save the figures
    fig, ax = plt.subplots(n_rows_interactions, n_cols_interactions, figsize=(20, 3 * n_rows_interactions), subplot_kw={'aspect': 'equal'})

    for i, height in enumerate(height_list):
        sr.plot(height, ax=ax[i, 0], color='w', cmap='gray', conv=1602, units='pN', text=True)
        es.plot(height, ax=ax[i, 1], color='k', cmap='gray', conv=1602, units='pN', text=True)
        vdw.plot(height, ax=ax[i, 2], color='k', cmap='gray', conv=1602, units='pN', text=True)
        static.plot(height, ax=ax[i, 3], color='w', cmap='gray', conv=1602, units='pN', text=True)
        relax.plot(height, ax=ax[i, 4], color='w', cmap='gray', conv=1602, units='pN', text=True)
        ax[i, 0].set_ylabel(fr'''$z_{{ts}}={height:.2f}\ \AA$''', fontsize=12)
        for j in range(n_cols_interactions):
            ax[i, j].set_xticks([])
            ax[i, j].set_yticks([])

    ax[0, 0].set_title('Short Range')
    ax[0, 1].set_title('Electrostatic')
    ax[0, 2].set_title('Van der Waals')
    ax[0, 3].set_title('Static')
    ax[0, 4].set_title('Relaxed')

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'{molec_name}_interactions.png'), bbox_inches='tight')
    plt.close()

    # Associate atomic positions with the grids for plotting
    static.atoms = sample
    vdw.atoms = sample
    es.atoms = sample
    sr.atoms = sample
    relax.atoms = sample

    # Plot interactions with atomic positions included
    fig, ax = plt.subplots(n_rows_interactions, n_cols_interactions, figsize=(20, 3 * n_rows_interactions), subplot_kw={'aspect': 'equal'})

    for i, height in enumerate(height_list):
        sr.plot(height, ax=ax[i, 0], color='w', cmap='gray', conv=1602, units='pN', text=True, show_atoms=True)
        es.plot(height, ax=ax[i, 1], color='k', cmap='gray', conv=1602, units='pN', text=True, show_atoms=True)
        vdw.plot(height, ax=ax[i, 2], color='k', cmap='gray', conv=1602, units='pN', text=True, show_atoms=True)
        static.plot(height, ax=ax[i, 3], color='w', cmap='gray', conv=1602, units='pN', text=True, show_atoms=True)
        relax.plot(height, ax=ax[i, 4], color='w', cmap='gray', conv=1602, units='pN', text=True, show_atoms=True)
        ax[i, 0].set_ylabel(fr'''$z_{{ts}}={height:.2f}\ \AA$''', fontsize=12)
        for j in range(n_cols_interactions):
            ax[i, j].set_xticks([])
            ax[i, j].set_yticks([])

    ax[0, 0].set_title('Short Range')
    ax[0, 1].set_title('Electrostatic')
    ax[0, 2].set_title('Van der Waals')
    ax[0, 3].set_title('Static')
    ax[0, 4].set_title('Relaxed')

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'{molec_name}_interactions_with_atoms.png'), bbox_inches='tight')
    plt.close()

    # Calculate frequency shift from relaxed grid
    fs = relax.get_fs(0.6, verbose=True)

    # Plot static, relaxed, and frequency shift images side by side
    fig, ax = plt.subplots(len(height_list), 3, figsize=(14, 3 * len(height_list)), subplot_kw={'aspect': 'equal'})

    for i, height in enumerate(height_list):
        static.plot(height, ax=ax[i, 0], cmap='gray', conv=1602, units='pN', text=True, color='w')
        relax.plot(height, ax=ax[i, 1], cmap='gray', conv=1602, units='pN', text=True, color='w')
        fs.plot(height, ax=ax[i, 2], cmap='gray', units='Hz', text=True, color='w')
        ax[i, 0].set_ylabel(fr'''$z_{{ts}}={height:.2f}\ \AA$''', fontsize=12)
        for j in range(3):
            ax[i, j].set_xticks([])
            ax[i, j].set_yticks([])

    ax[0, 0].set_title('Static Image')
    ax[0, 1].set_title('Relaxed')
    ax[0, 2].set_title('Frequency Shift')

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'{molec_name}_static_relaxed_fs.png'), bbox_inches='tight')
    plt.close()

    # Plot static, relaxed, and frequency shift images side by side
    fig, ax = plt.subplots(len(height_list), 3, figsize=(14, 3 * len(height_list)), subplot_kw={'aspect': 'equal'})

    for i, height in enumerate(height_list):
        static.plot(height, ax=ax[i, 0], cmap='gray', conv=1602, units='pN', text=True, show_atoms=True, color='w')
        relax.plot(height, ax=ax[i, 1], cmap='gray', conv=1602, units='pN', text=True, show_atoms=True, color='w')
        fs.plot(height, ax=ax[i, 2], cmap='gray', units='Hz', text=True, show_atoms=True, color='w')
        ax[i, 0].set_ylabel(fr'''$z_{{ts}}={height:.2f}\ \AA$''', fontsize=12)
        for j in range(3):
            ax[i, j].set_xticks([])
            ax[i, j].set_yticks([])

    ax[0, 0].set_title('Static Image')
    ax[0, 1].set_title('Relaxed')
    ax[0, 2].set_title('Frequency Shift')

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'{molec_name}_static_relaxed_fs_with_atoms.png'), bbox_inches='tight')
    plt.close()

if __name__ == '__main__':
    main()


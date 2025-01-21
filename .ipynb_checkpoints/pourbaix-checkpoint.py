#!/usr/bin/env python3

import os
import sys
import warnings
import pandas as pd
import matplotlib.pyplot as plt
from pymatgen.ext.matproj import MPRester
from pymatgen.core.ion import Ion
from pymatgen.analysis.pourbaix_diagram import PourbaixEntry, PourbaixDiagram, PourbaixPlotter
from pymatgen.entries.computed_entries import ComputedEntry

# Ignore warning messages
warnings.filterwarnings('ignore')

# Set plot file name
PLOT_NAME = 'Pourbaix_PBE+U_U2p75_allLayered'

# Load Materials Project API key from environment variable
API_KEY = os.getenv('MAPI_KEY')
if not API_KEY:
    sys.exit("Error: MAPI_KEY environment variable not set.")

mpr = MPRester(API_KEY)

# Define constants
kJmol = 96.485

df = pd.read_csv('/pscratch/sd/j/jiuy97/6_MNC/figures/pourbaix/1Fe_energies.tsv', delimiter='\t', index_col=0)
df.index = df.index.str.upper()
print(df.head())

# def get_pourbaix_energy(comp):
    
#     SAC_formation_energy=0.1+2*(-0.44)  -2*237/kJmol
#     print('SAC_formation_energy : ', SAC_formation_energy)
#     return pourbaix_energy

def get_solid_entries():
    """Generate solid entries."""
    ion_dict_solids_mnc = {
    # 'FeNC-OH': SAC_formation_energy,
    # 'FeNC-OH': SAC_formation_energy,
    # 'FeO3H2': SAC_formation_energy
    }

    solid_entries = []
    for key, energy in ion_dict_solids_mnc.items():
        comp = Ion.from_formula(key)
        entry = PourbaixEntry(ComputedEntry(comp, energy))
        solid_entries.append(entry)

    # for entry in solid_entries:
    #     entry.nPhi -= 2
    #     entry.nH2O -= 

    return solid_entries

def get_ion_entries():
    """Fetch ion entries from the Materials Project API."""
    try:
        ion_entries = mpr.get_pourbaix_entries(["Fe"])
    except Exception as e:
        sys.exit(f"Error fetching ion data from Materials Project API: {e}")

    return ion_entries

def plot_pourbaix(entries):
    """Plot and save Pourbaix diagram."""
    pourbaix = PourbaixDiagram(entries, filter_solids=False)
    plotter = PourbaixPlotter(pourbaix)

    # Generate the Pourbaix plot and get the Axes object
    ax = plotter.get_pourbaix_plot(limits=[[-2, 16], [-2, 4]])
    
    # Customize the plot
    for line in ax.lines:
        line.set_linewidth(1.0)  # Adjust line thickness
    
    for text in ax.texts:
        text.set_fontsize(14)  # Adjust phase label font size
    
    # Set axis labels and tick font size
    ax.set_xlabel("pH", fontsize=14)
    ax.set_ylabel("Potential (V vs SHE)", fontsize=14)
    ax.tick_params(axis='both', labelsize=14)
    
    fig = ax.figure
    fig.set_size_inches((8, 7))
    plt.tight_layout()

    plt.savefig(PLOT_NAME + '.pdf')
    plt.savefig(PLOT_NAME + '.png')
    plt.show()

def main():
    """Main execution function."""
    print('\n################## Solid Entries ##########################################\n')
    solid_entries = get_solid_entries()
    for entry in solid_entries:
        print(entry)

    print('\n################## Ion Entries ##########################################\n')
    ion_entries = get_ion_entries()
    for entry in ion_entries:
        print(entry)

    all_entries = solid_entries + ion_entries
    print("\nTotal Entries:", len(all_entries))

    for entry in all_entries:
        if entry.npH != entry.nPhi:
            print(entry)
    
    plot_pourbaix(all_entries)

if __name__ == "__main__":
    main()
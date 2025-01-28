#!/usr/bin/env python3

import os
import sys
import warnings
import pandas as pd
import matplotlib.pyplot as plt
from pymatgen.ext.matproj import MPRester
from pymatgen.core.ion import Ion
from pymatgen.core.composition import Composition
from pymatgen.analysis.pourbaix_diagram import PourbaixEntry, PourbaixDiagram, PourbaixPlotter
from pymatgen.analysis.pourbaix_diagram import IonEntry, PDEntry, ComputedEntry

warnings.filterwarnings('ignore')

png_name = '3Mo_pourbaix'
tsv_name = '3Mo_energies.tsv'

API_KEY = os.getenv('MAPI_KEY')
if not API_KEY:
    sys.exit("Error: MAPI_KEY environment variable not set.")
mpr = MPRester(API_KEY)
mpr_entries = mpr.get_pourbaix_entries(['Mo'])

mpr1_entries = []
mpr2_entries = []

for entry in mpr_entries:
    if 'ion' in entry.entry_id and entry.npH - entry.nPhi > 0:
        mpr1_entries.append(entry)
    else:
        mpr2_entries.append(entry)
        
kJmol = 96.485
calmol = 23.061
water = 2.4583 # the standard Gibbs free energy of formation of water

# gas
h2 = -6.77149190
h2o = -14.23091949

zpeh2o = 0.560
cvh2o = 0.103
tsh2o = 0.675

zpeh2 = 0.268
cvh2 = 0.0905
tsh2 = 0.408

gh2o = h2o + zpeh2o - tsh2o + cvh2o
gh2 = h2 + zpeh2 - tsh2 + cvh2

n2 = -16.64503942
zpen2 = 0.098
tsn2 = 0.592
gn2 = n2 + zpen2 - tsn2

# solid
gc = -.37429454E+02/4

N4C26 = -.27195317E+03
H2N4C26 = -.28183609E+03

# ads
zpeoh = 0.376
cvoh = 0.042
tsoh = 0.066

zpeo = 0.064
cvo = 0.034
tso = 0.060

zpeooh = 0.471
cvooh = 0.077
tsooh = 0.134

dgo = zpeo + cvo - tso
dgoh = zpeoh + cvoh - tsoh
dgooh = zpeooh + cvooh - tsooh
dgh = dgoh - dgo

elements_data = {
    "Ti": {"electrode_potential": -1.63, "cation_charge": 2},  # Ti^2+ + 2e- → Ti (Ti^2+ prioritized over Ti^4+)
    "V":  {"electrode_potential": -1.18, "cation_charge": 2},  # V^2+ + 2e- → V (V^2+ prioritized over V^3+)
    "Cr": {"electrode_potential": -0.91, "cation_charge": 2},  # Cr^2+ + 2e- → Cr (Cr^2+ prioritized over Cr^3+)
    "Mn": {"electrode_potential": -1.18, "cation_charge": 2},  # Mn^2+ + 2e- → Mn
    "Fe": {"electrode_potential": -0.44, "cation_charge": 2},  # Fe^2+ + 2e- → Fe (Fe^2+ prioritized over Fe^3+)
    "Co": {"electrode_potential": -0.28, "cation_charge": 2},  # Co^2+ + 2e- → Co (Co^2+ prioritized over Co^3+)
    "Ni": {"electrode_potential": -0.25, "cation_charge": 2},  # Ni^2+ + 2e- → Ni
    "Cu": {"electrode_potential": 0.34,  "cation_charge": 2},  # Cu^2+ + 2e- → Cu

    "Zr": {"electrode_potential": -1.45, "cation_charge": 2},  # Zr^2+ + 2e- → Zr (Zr^2+ prioritized over Zr^4+)
    "Nb": {"electrode_potential": -1.00, "cation_charge": 3},  # Nb^3+ + 3e- → Nb (Nb^3+ prioritized over Nb^5+)
    "Mo": {"electrode_potential": -0.20, "cation_charge": 3},  # Mo^3+ + 3e- → Mo (Mo^3+ prioritized over Mo^6+)
    "Tc": {"electrode_potential": 0.74,  "cation_charge": 7},  # Tc^7+ + 7e- → Tc (No lower oxidation state available)
    "Ru": {"electrode_potential": 0.45,  "cation_charge": 2},  # Ru^2+ + 2e- → Ru (Ru^2+ prioritized over Ru^3+)
    "Rh": {"electrode_potential": 0.76,  "cation_charge": 2},  # Rh^2+ + 2e- → Rh (Rh^2+ prioritized over Rh^3+)
    "Pd": {"electrode_potential": 0.95,  "cation_charge": 2},  # Pd^2+ + 2e- → Pd

    "Hf": {"electrode_potential": -1.55, "cation_charge": 4},  # Hf^4+ + 4e- → Hf (No lower oxidation state available)
    "Ta": {"electrode_potential": -1.10, "cation_charge": 3},  # Ta^3+ + 3e- → Ta (Ta^3+ prioritized over Ta^5+)
    "W":  {"electrode_potential": -0.11, "cation_charge": 4},  # W^4+ + 4e- → W (W^4+ prioritized over W^6+)
    "Re": {"electrode_potential": 0.30,  "cation_charge": 4},  # Re^4+ + 4e- → Re (Re^4+ prioritized over Re^7+)
    "Os": {"electrode_potential": 0.40,  "cation_charge": 2},  # Os^2+ + 2e- → Os (Os^2+ prioritized over Os^4+)
    "Ir": {"electrode_potential": 1.16,  "cation_charge": 2},  # Ir^2+ + 2e- → Ir (Ir^2+ prioritized over Ir^3+)
    "Pt": {"electrode_potential": 1.20,  "cation_charge": 2},  # Pt^2+ + 2e- → Pt
}

potential = elements_data['Mo']['electrode_potential']
charge = elements_data['Mo']['cation_charge']

metal_path = './metals.tsv'
metal_df = pd.read_csv(metal_path, delimiter='\t', index_col=0)
gm = metal_df.loc['Mo', 'energy']

df = pd.read_csv(tsv_name, delimiter='\t', index_col=0)

df['name'] = 'MoNC(' + df.index.str.upper() + ')'
df['comp'] = 'MoXO' + df.index.str.upper().str.replace("-", "")
df['comp'] = df['comp'].str.replace('MoXOVAC', 'H2XO')
df['name'] = df['name'].str.replace('MoNC(VAC)', 'Mo⁺²+H₂NC', regex=False)
df['comp'] = df['comp'].str.replace('MoXOCLEAN', 'MoXO')
df['name'] = df['name'].str.replace('MoNC(CLEAN)', 'MoNC(clean)')
df['comp'] = df['comp'].str.replace('MoXOMH', 'MoXOH')
df['comp'] = df['comp'].str.replace('MoXONH', 'MoXOH')

# df['energy'] = df['dG'] + df.loc['clean', 'G'] - gm - 2 * gn2 - 26 * gc - water * (df['#O'] + df['#OH'] + df['#OOH']*2)
# df['energy'] = df['dG'] + df.loc['clean', 'G'] - gm - N4C26 - water * (df['#O'] + df['#OH'] + df['#OOH']*2)
df['energy'] = df['dG'] + df.loc['clean', 'G'] + gh2 - gm - H2N4C26 - 2 * dgh - water * (df['#O'] + df['#OH'] + df['#OOH']*2)

df = df.drop(index='vac')
df = df.drop(index='o-oh')
df = df.drop(index='o-ooh')
df = df.drop(index='ooh-o')
df = df.drop(index='oh-ooh')
df = df.drop(index='ooh-oh')
df = df.drop(index='ooh-ooh')
df = df.dropna()
print(df)

def get_ref_entries():
    ref_entries = []
    
    refs={
        'Mo': 'Mo(s)',
        # 'N2': 'N2(g)',
        # 'C': 'C(s)',
        # 'X': 'N4C26',
        # 'H2X': 'H2NC(vac)',
        }
    
    for comp, name in refs.items():
        entry = PourbaixEntry(PDEntry(comp, 0.0, name=name))
        ref_entries.append(entry)
        
    entry = PourbaixEntry(PDEntry('H2XO', -water, name='H2NC(vac)'))
    ref_entries.append(entry)
    
    return ref_entries
    
def get_sac_entries():
    sac_entries = []
    
    for index, row in df.iterrows(): 
        entry = PourbaixEntry(ComputedEntry(row['comp'], row['energy']-water, entry_id=row['name']))
        sac_entries.append(entry)
        
    energy = H2N4C26 + 2 * dgh - gh2 - N4C26
    entry = PourbaixEntry(ComputedEntry('X', -energy-water, entry_id='NC(vac)'))
    sac_entries.append(entry)
    
    return sac_entries

def get_solid_entries():
    solid_entries = []
    solids={
        'Mo': 0,
        'MoO2' : -120.000/calmol,
        # 'MoO3': -227.000/calmol,
        'MoO3': -161.950/calmol,
        }
    
    for solid, energy in solids.items():
        entry = PourbaixEntry(PDEntry(solid, energy))
        solid_entries.append(entry)

    return solid_entries

def get_ion_entries():
    ion_entries = []
    ions={
        'Mo+++': -13.800/calmol,
        # 'HMoO4-': -213.600/calmol,
        'MoO4--': -205.420/calmol,
        }
    
    for ion, energy in ions.items():
        comp = Ion.from_formula(ion)
        entry = PourbaixEntry(IonEntry(comp, energy), concentration=1e-6)
        ion_entries.append(entry)

    return ion_entries
    
def plot_pourbaix(entries, png_name):
    pourbaix = PourbaixDiagram(entries, filter_solids=False)
    plotter = PourbaixPlotter(pourbaix)

    ax = plotter.get_pourbaix_plot(limits=[[-2, 16], [-2, 4]])
    
    for line in ax.lines:
        line.set_linewidth(1.0)
    for text in ax.texts:
        text.set_fontsize(14)
    
    ax.set_xlabel("pH", fontsize=14)
    ax.set_ylabel("Potential (V vs SHE)", fontsize=14)
    ax.tick_params(axis='both', labelsize=14)
    
    fig = ax.figure
    fig.set_size_inches((8, 8))

    plt.savefig(png_name, bbox_inches='tight')
    # plt.close()
    plt.show()

def main():
    print('\n################## ReMorence Entries ##########################################\n')
    ref_entries = get_ref_entries()
    for entry in ref_entries:
        print(entry)
    
    print('\n################## SAC Entries ##########################################\n')
    sac_entries = get_sac_entries()
    for entry in sac_entries:
        print(entry)

    print('\n################## Solid Entries ##########################################\n')
    solid_entries = get_solid_entries()
    for entry in solid_entries:
        print(entry)
        
    print('\n################## Ion Entries ##########################################\n')
    ion_entries = get_ion_entries()
    for entry in ion_entries:
        print(entry)

    all_entries = ref_entries + sac_entries + solid_entries + ion_entries
    print("\nTotal Entries:", len(all_entries))
    
    all_entries = ref_entries + sac_entries
    plot_pourbaix(all_entries, f'{png_name}_sac4.png')
    
    # plot_pourbaix(solid_entries, f'{png_name}_solid.png')
    # plot_pourbaix(ion_entries, f'{png_name}_ion.png')
    # all_entries = solid_entries + ion_entries
    # plot_pourbaix(all_entries, f'{png_name}_exp.png')

    # plot_pourbaix(mpr1_entries, f'{png_name}_mpr1.png')
    # plot_pourbaix(mpr2_entries, f'{png_name}_mpr2.png')
    
    all_entries = ref_entries + sac_entries + solid_entries + ion_entries
    plot_pourbaix(all_entries, f'{png_name}_bulk4.png')


if __name__ == "__main__":
    main()
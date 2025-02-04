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
png_name = 'iron_sac_vacancy_reference'

# API_KEY = os.getenv('MAPI_KEY')
# if not API_KEY:
#     sys.exit("Error: MAPI_KEY environment variable not set.")
# mpr = MPRester(API_KEY)
# mpr_entries = mpr.get_pourbaix_entries(['Fe'])

# mpr1_entries = []
# mpr2_entries = []

# for entry in mpr_entries:
#     if 'ion' in entry.entry_id and entry.npH - entry.nPhi > 0:
#         mpr1_entries.append(entry)
#     else:
#         mpr2_entries.append(entry)
        
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

gh = gh2 / 2
go = gh2o - gh2
goh = gh2o - gh2 / 2
gooh = 2 * gh2o - 1.5 * gh2

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

# metal_path = './metals.tsv'
# metal_df = pd.read_csv(metal_path, delimiter='\t', index_col=0)
# gm = metal_df.loc['Fe', 'energy']
gm = -5.041720865 # Fe

df = pd.DataFrame()
df.loc['clean', ['E', '#H', '#O', '#OH', '#OOH']] = [-280.18, 0, 0, 0, 0]
df.loc['mh', ['E', '#H', '#O', '#OH', '#OOH']] = [-282.56, 1, 0, 0, 0]
df.loc['nh', ['E', '#H', '#O', '#OH', '#OOH']] = [-282.64, 1, 0, 0, 0]
df.loc['o', ['E', '#H', '#O', '#OH', '#OOH']] = [-285.95, 0, 1, 0, 0]
df.loc['oh', ['E', '#H', '#O', '#OH', '#OOH']] = [-290.98, 0, 0, 1, 0]
df.loc['ohoh', ['E', '#H', '#O', '#OH', '#OOH']] = [-300.74, 0, 0, 2, 0]
df.loc['oh-oh', ['E', '#H', '#O', '#OH', '#OOH']] = [-300.63, 0, 0, 2, 0]
df.loc['ooh', ['E', '#H', '#O', '#OH', '#OOH']] = [-295.14, 0, 0, 0, 1]
df.loc['oho', ['E', '#H', '#O', '#OH', '#OOH']] = [-295.23, 0, 1, 1, 0]
df.loc['oh-o', ['E', '#H', '#O', '#OH', '#OOH']] = [-295.53, 0, 1, 1, 0]
df.loc['o-o', ['E', '#H', '#O', '#OH', '#OOH']] = [-289.68, 0, 2, 0, 0]

df['comp'] = 'FeX' + df.index.str.upper().str.replace("-", "")
df['comp'] = df['comp'].str.replace('FeXCLEAN', 'FeX')
df['comp'] = df['comp'].str.replace('FeXMH', 'FeXH')
df['comp'] = df['comp'].str.replace('FeXNH', 'FeXH')

df['G'] = df['E'] + dgh * df['#H'] + dgo * df['#O'] + dgoh * df['#OH'] + dgooh * df['#OOH']
df['dG'] = df['G'] - df.loc['clean', 'E'] - gh * df['#H'] - go * df['#O'] - goh * df['#OH'] - gooh * df['#OOH']
df['energy'] = df['dG'] + df.loc['clean', 'G'] + 2 * gh - gm - H2N4C26 - 2 * dgh - water * (df['#O'] + df['#OH'] + df['#OOH']*2)

df_print = pd.DataFrame()
df_print = df
df_print['E'] = df_print['E'].astype(float).round(2)
df_print['G'] = df_print['G'].astype(float).round(2)
df_print['dG'] = df_print['dG'].astype(float).round(2)
df_print['energy'] = df_print['energy'].astype(float).round(2)
print(df_print)

def get_ref_entries():
    ref_entries = []
    
    refs={
        'Fe': 'Fe(s)',
        # 'N2': 'N2(g)',
        # 'C': 'C(s)',
        # 'X': 'N4C26',
        'H2X': 'H2NC(vac)',
        }
    
    for comp, name in refs.items():
        entry = PourbaixEntry(PDEntry(comp, 0.0, name=name))
        ref_entries.append(entry)

    return ref_entries
    
def get_sac_entries():
    sac_entries = []
    
    for index, row in df.iterrows():    
        entry = PourbaixEntry(ComputedEntry(row['comp'], row['energy']))
        sac_entries.append(entry)
        
    energy = H2N4C26 + 2 * dgh - gh2 - N4C26
    entry = PourbaixEntry(ComputedEntry('X', -energy, entry_id='NC(vac)'))
    sac_entries.append(entry)
    
    return sac_entries  

def get_solid_entries():
    solid_entries = []
    solids={
        'Fe': 0,
        'FeO' : -58.880/calmol,
        'Fe3O4': -242.400/calmol,
        'Fe2O3': -177.100/calmol,
        'Fe2O3': -161.930/calmol,
        'Fe(OH)2': -115.570/calmol,
        'Fe(OH)3': -166.000/calmol,
        }
    
    for solid, energy in solids.items():
        entry = PourbaixEntry(PDEntry(solid, energy))
        solid_entries.append(entry)

    return solid_entries

def get_ion_entries():
    ion_entries = []
    ions={
        'Fe++': -20.300/calmol,
        'HFeO2-': -90.627/calmol,
        'Fe+++': -2.530/calmol,
        'FeOH++': -55.910/calmol,
        'Fe(OH)2+': -106.200/calmol,
        # 'FeO4-': -111685/calmol,
        }
    
    for ion, energy in ions.items():
        comp = Ion.from_formula(ion)
        entry = PourbaixEntry(IonEntry(comp, energy), concentration=1e-6)
        ion_entries.append(entry)

    return ion_entries
    
def plot_pourbaix(entries, png_name):
    pourbaix = PourbaixDiagram(entries, filter_solids=False)
    plotter = PourbaixPlotter(pourbaix)

    ax = plotter.get_pourbaix_plot(limits=[[-2, 16], [-1, 3]])
    
    for line in ax.lines:
        line.set_linewidth(1.0)
    for text in ax.texts:
        text.set_fontsize(14)
    
    ax.set_xlabel("pH", fontsize=14)
    ax.set_ylabel("Potential (V vs SHE)", fontsize=14)
    ax.tick_params(axis='both', labelsize=14)
    
    fig = ax.figure
    fig.set_size_inches((8, 6))
    
    plt.savefig(png_name, dpi=300, bbox_inches='tight') #, transparent=True)
    # plt.close()
    plt.show()

def main():
    print('\n################## Reference Entries ##########################################\n')
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
    plot_pourbaix(all_entries, f'{png_name}_sac.png')
    
    # plot_pourbaix(solid_entries, f'{png_name}_solid.png')
    # plot_pourbaix(ion_entries, f'{png_name}_ion.png')
    # all_entries = solid_entries + ion_entries
    # plot_pourbaix(all_entries, f'{png_name}_exp.png')

    # plot_pourbaix(mpr1_entries, f'{png_name}_mpr1.png')
    # plot_pourbaix(mpr2_entries, f'{png_name}_mpr2.png')
    
    all_entries = ref_entries + sac_entries + solid_entries + ion_entries
    plot_pourbaix(all_entries, f'{png_name}_bulk.png')

if __name__ == "__main__":
    main()

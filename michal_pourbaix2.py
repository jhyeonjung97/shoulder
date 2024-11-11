import os
import re
import glob
import numpy as np
import pandas as pd
from matplotlib import rc
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

# dirs = glob.glob("/pscratch/sd/j/jiuy97/6_MNC/pourbaix/*_*/")
dirs = ["/pscratch/sd/j/jiuy97/6_MNC/pourbaix/1_Fe/",
        "/pscratch/sd/j/jiuy97/6_MNC/pourbaix/2_Co/",
        "/pscratch/sd/j/jiuy97/6_MNC/pourbaix/3_Mo/"]
# main_dirs = ["clean", "h", "oh", "o", 
#              "ohoh", "oh-oh", "ohooh", "oohoh", "oh-ooh", "ooh-oh",
#              "ooh", "oho", "oh-o", "o-oh", "oo", "o-o",
#              "oooh", "ooho", "o-ooh", "ooh-o", "oohooh", "ooh-ooh"]
main_dirs = ["clean", "h", "oh", "o", 
             "ohoh", "oh-oh", "ohooh", "oohoh", "ooh-oh", # "oh-ooh"
             "ooh", "oho", "o-oh", "o-o", # "oh-o", "oo",
             "oooh", "ooho", "ooh-o", "oohooh", "ooh-ooh"] # "o-ooh"
sub_dirs = ["HS1", "HS5", "IS1", "IS5", "LS1", "LS5"]

# Regular expression to match E0 values in scientific notation
e0_pattern = re.compile(r"E0=\s*(-?\.\d+E[+-]?\d+)")
energy_pattern = re.compile(r"energy\s*=\s*(-?\d+\.\d+)")
    
fig_width_pt = 1.8 * 246.0  # LaTeX column width
inches_per_pt = 1.0 / 72.27  # Convert pt to inches
golden_mean = (np.sqrt(5) - 1.0) / 2.0  # Aesthetic ratio
fig_width = fig_width_pt * inches_per_pt  # Width in inches
fig_height = fig_width * golden_mean  # Height in inches
fig_size = [fig_width, fig_height]

font_size = 10
tick_font_size = 10

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['DejaVu Sans']
plt.rcParams['font.size'] = font_size
plt.rcParams['axes.labelsize'] = 2 * font_size
plt.rcParams['legend.fontsize'] = font_size
plt.rcParams['xtick.labelsize'] = tick_font_size
plt.rcParams['ytick.labelsize'] = tick_font_size
plt.rcParams['mathtext.default'] = 'regular'
plt.rcParams['lines.linewidth'] = 1.0

# Constants and variables for plotting
Umin, Umax = -0.5, 2.5
kbt = 0.0256 
const = kbt * np.log(10)
kjmol = 96.485

pH = np.arange(0, 14, 0.10)
U = np.arange(Umin, Umax, 0.05)
Umax2 = Umax + 0.06 * 14
U2 = np.arange(Umin, Umax2, 0.05)

# gas
h2 = -6.77149190
h2o = -14.23091949

zpeh2o = 0.560
cvh2o = 0.103
tsh2o = 0.675

zpeh2 = 0.268
cvh2 = 0.0905
tsh2 = 0.408

gh2o = h2o + cvh2o - tsh2o + zpeh2o
gh2 = h2 + cvh2 - tsh2 + zpeh2

gh = gh2 / 2
go = gh2o - gh2
goh = gh2o - gh2 / 2
gooh = go + goh

dgh2o = zpeh2o + cvh2o - tsh2o
dgh2 = zpeh2 + cvh2 - tsh2

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

dso = dgo - (dgh2o - dgh2)
dsoh = dgoh - (dgh2o - 0.5 * dgh2)
dsooh = dgooh - (2 * dgh2o - 1.5 * dgh2)
dsh = dsoh - dso

color = ['turquoise', 'green', 'red', 'blue', 'gray', 'gold', 'purple', 'pink', 'darkorange',
         'lime', 'olive', 'yellowgreen', 'violet', 'navy', 'brown', 'teal', 'deeppink',
         'cyan', 'dodgerblue', 'steelblue', 'darkslategrey']
pH2 = np.arange(0, 14.01, 0.01)

# Function to find the lowest E0 value in each subdirectory
def find_min_e0(main_dir, sub_dirs):
    min_e0 = None
    for sub_dir in sub_dirs:
        oszicar_path = os.path.join(main_dir, sub_dir, "OSZICAR")
        if os.path.isfile(oszicar_path):
            with open(oszicar_path, 'r') as file:
                for line in file:
                    match = e0_pattern.search(line)
                    if match:
                        e0_value = float(match.group(1))
                        if min_e0 is None or e0_value < min_e0:
                            min_e0 = e0_value
    return min_e0
    
# Function to extract energy from DONE in most_stable or find min_e0 as fallback
def get_energy(main_dir, sub_dirs):
    most_stable_dir = os.path.join(main_dir, "most_stable")
    done_path = os.path.join(most_stable_dir, "DONE")
    if os.path.isfile(done_path):
        with open(done_path, 'r') as file:
            for line in file:
                match = energy_pattern.search(line)
                if match:
                    return float(match.group(1))
    # Fallback to finding minimum E0 value
    return find_min_e0(main_dir, sub_dirs)
    if min_e0 is None:
        print(f"Warning: No valid energy found in {main_dir}.")
    return min_e0
    
def addO(x, y):
    return -(h2o - h2) - 2 * (y + x * const) + dso

def addOH(x, y):
    return -(h2o - 0.5 * h2) - (y + x * const) + dsoh

def addOOH(x, y):
    return -(2 * h2o - 1.5 * h2) - 3 * (y + x * const) + dsooh

def addH(x, y):
    return -0.5 * h2 + 1 * (y + x * const) + dsh

def dg(i, x, y):
    if surfs[i][0] is None:
        return None
    return (surfs[i][0] 
            - surfs[0][0] 
            + surfs[i][1] * addH(x, y) 
            + surfs[i][2] * addO(x, y) 
            + surfs[i][3] * addOH(x, y) 
            + surfs[i][4] * addOOH(x, y))
    
def overpotential_oer(int1, int2, int3, int4, df, overpotentials):
    ints = [int1, int2, int3, int4]
    for i, int in enumerate(ints):
        if isinstance(int, tuple):
            if df.loc[int[0], 'E'] < df.loc[int[1], 'E']:
                ints[i] = int[0]
            else:
                ints[i] = int[1] 
    int1, int2, int3, int4 = ints
    print(int1, int2, int3, int4)
    dG1 = df.loc[int2, 'dG'] - df.loc[int1, 'dG']
    dG2 = df.loc[int3, 'dG'] - df.loc[int2, 'dG']
    dG3 = df.loc[int4, 'dG'] - df.loc[int3, 'dG']
    dG4 = 4.92 - dG1 - dG2 - dG3
    onsetP = max(dG1, dG2, dG3, dG4)
    overP = onsetP - 1.23
    overpotentials.update({'int1': int1, 'int2': int2, 'int3': int3, 'int4': int4, 'overP': overP, 'onsetP': onsetP})
    print(f"Updated overpotentials: {overpotentials}")

for dir in dirs:
    os.chdir(dir)
    print(dir)
    basename = os.path.basename(os.path.normpath(dir))
    A, B = basename.split('_', 1)
    df = pd.DataFrame()
    overpotentials = {'int1': {}, 'int2': {}, 'int3': {}, 'int4': {}, 'overP': {}, 'onsetP': {}}
    
    # Iterate through each main directory to extract E0 values and plot
    for main_dir in main_dirs:
        min_e0 = get_energy(main_dir, sub_dirs)
    
        if min_e0 is None:
            print(f"Missing data in directory '{main_dir}' for plotting.")
            df.loc[main_dir, 'E'] = np.nan
            continue
        else:
            df.loc[main_dir, 'E'] = min_e0

    df.loc['clean', ['#H', '#O', '#OH', '#OOH']] = [0, 0, 0, 0] # [energy, #Hs, #Os, #OHs, #OOHs]
    df.loc['h', ['#H', '#O', '#OH', '#OOH']] = [1, 0, 0, 0]
    df.loc['o', ['#H', '#O', '#OH', '#OOH']] = [0, 0, 1, 0]
    df.loc['oh', ['#H', '#O', '#OH', '#OOH']] = [0, 1, 0, 0]
    df.loc['ohoh', ['#H', '#O', '#OH', '#OOH']] = [0, 0, 2, 0]
    df.loc['oh-oh', ['#H', '#O', '#OH', '#OOH']] = [0, 0, 2, 0]
    df.loc['ohooh', ['#H', '#O', '#OH', '#OOH']] = [0, 0, 1, 1]
    df.loc['oohoh', ['#H', '#O', '#OH', '#OOH']] = [0, 0, 1, 1]
    df.loc['oh-ooh', ['#H', '#O', '#OH', '#OOH']] = [0, 0, 1, 1]
    df.loc['ooh-oh', ['#H', '#O', '#OH', '#OOH']] = [0, 0, 1, 1]
    df.loc['ooh', ['#H', '#O', '#OH', '#OOH']] = [0, 0, 0, 1]
    df.loc['oho', ['#H', '#O', '#OH', '#OOH']] = [0, 1, 1, 0]
    df.loc['oh-o', ['#H', '#O', '#OH', '#OOH']] = [0, 1, 1, 0]
    df.loc['o-oh', ['#H', '#O', '#OH', '#OOH']] = [0, 1, 1, 0]
    df.loc['oo', ['#H', '#O', '#OH', '#OOH']] = [0, 2, 0, 0]
    df.loc['o-o', ['#H', '#O', '#OH', '#OOH']] = [0, 2, 0, 0]
    df.loc['oooh', ['#H', '#O', '#OH', '#OOH']] = [0, 1, 0, 1]
    df.loc['ooho', ['#H', '#O', '#OH', '#OOH']] = [0, 1, 0, 1]
    df.loc['o-ooh', ['#H', '#O', '#OH', '#OOH']] = [0, 1, 0, 1]
    df.loc['ooh-o', ['#H', '#O', '#OH', '#OOH']] = [0, 1, 0, 1]
    df.loc['oohooh', ['#H', '#O', '#OH', '#OOH']] = [0, 0, 0, 2]
    df.loc['ooh-ooh', ['#H', '#O', '#OH', '#OOH']] = [0, 0, 0, 2]

    df['G'] = df['E'] + dgoh * df['#OH'] + dgo * df['#O'] + dgooh * df['#OOH'] # + dgh * df['#H']
    df['dG'] = df['G'] - df.loc['clean', 'E'] - gh * df['#H'] - goh * df['#OH'] - go * df['#O'] - gooh * df['#OOH']

    overpotential_oer('clean', 'oh', 'o', 'ooh', df, overpotentials)
    if A == '1' and B == 'Fe':
        overpotential_oer('oh', 'o', ('o-oh', 'oh-o'), ('ooh-oh', 'oh-ooh'), df, overpotentials)
        overpotential_oer('o', ('o-oh', 'oh-o'), 'o-o', ('ooh-o', 'o-ooh'), df, overpotentials)
        overpotential_oer('o', ('o-oh', 'oh-o'), ('ooh-oh', 'oh-ooh'), ('ooh-o', 'o-ooh'), df, overpotentials)
    elif A == '2' and B == 'Co':
        overpotential_oer('oh', 'oh-oh', ('o-oh', 'oh-o'), ('ooh-oh', 'oh-ooh'), df, overpotentials)
        overpotential_oer('o', ('o-oh', 'oh-o'), 'o-o', ('ooh-o', 'o-ooh'), df, overpotentials)
    elif A == '3' and B == 'Mo':
        overpotential_oer('oh', 'o', 'oho', ('oohoh', 'ohooh'), df, overpotentials)
        overpotential_oer('o', 'oho', 'oo', ('oooh', 'ooho'), df, overpotentials)

    # Define surfaces with extracted E0 values
    surfs = [
        df.loc['clean', ['E', '#H', '#O', '#OH', '#OOH']].values,
        df.loc['h', ['E', '#H', '#O', '#OH', '#OOH']].values,
        df.loc['o', ['E', '#H', '#O', '#OH', '#OOH']].values,
        df.loc['oh', ['E', '#H', '#O', '#OH', '#OOH']].values,
        df.loc['ohoh', ['E', '#H', '#O', '#OH', '#OOH']].values,
        df.loc['oh-oh', ['E', '#H', '#O', '#OH', '#OOH']].values,
        df.loc['ohooh', ['E', '#H', '#O', '#OH', '#OOH']].values,
        df.loc['oohoh', ['E', '#H', '#O', '#OH', '#OOH']].values,
        # df.loc['oh-ooh', ['E', '#H', '#O', '#OH', '#OOH']].values,
        # df.loc['ooh-oh', ['E', '#H', '#O', '#OH', '#OOH']].values,  
        df.loc['ooh', ['E', '#H', '#O', '#OH', '#OOH']].values,
        df.loc['oho', ['E', '#H', '#O', '#OH', '#OOH']].values,
        # df.loc['oh-o', ['E', '#H', '#O', '#OH', '#OOH']].values,  
        df.loc['o-oh', ['E', '#H', '#O', '#OH', '#OOH']].values,
        # df.loc['oo', ['E', '#H', '#O', '#OH', '#OOH']].values,  
        df.loc['o-o', ['E', '#H', '#O', '#OH', '#OOH']].values,
        # df.loc['oooh', ['E', '#H', '#O', '#OH', '#OOH']].values,  
        # df.loc['ooho', ['E', '#H', '#O', '#OH', '#OOH']].values,  
        # df.loc['o-ooh', ['E', '#H', '#O', '#OH', '#OOH']].values,  
        # df.loc['ooh-o', ['E', '#H', '#O', '#OH', '#OOH']].values,  
        # df.loc['oohooh', ['E', '#H', '#O', '#OH', '#OOH']].values,  
        # df.loc['ooh-ooh', ['E', '#H', '#O', '#OH', '#OOH']].values,  
    ]

    if A == '1' and B == 'Fe':
        # surfs.append(df.loc['oh-ooh', ['E', '#H', '#O', '#OH', '#OOH']].values)
        surfs.append(df.loc['ooh-oh', ['E', '#H', '#O', '#OH', '#OOH']].values)
        # surfs.append(df.loc['o-ooh', ['E', '#H', '#O', '#OH', '#OOH']].values)
        surfs.append(df.loc['ooh-o', ['E', '#H', '#O', '#OH', '#OOH']].values)
        surfs.append(df.loc['ooh-ooh', ['E', '#H', '#O', '#OH', '#OOH']].values)
    elif A == '2' and B == 'Co':
        # surfs.append(df.loc['oh-ooh', ['E', '#H', '#O', '#OH', '#OOH']].values)
        surfs.append(df.loc['ooh-oh', ['E', '#H', '#O', '#OH', '#OOH']].values)
        # surfs.append(df.loc['o-ooh', ['E', '#H', '#O', '#OH', '#OOH']].values)
        surfs.append(df.loc['ooh-o', ['E', '#H', '#O', '#OH', '#OOH']].values)
        surfs.append(df.loc['ooh-ooh', ['E', '#H', '#O', '#OH', '#OOH']].values)
    elif A == '3' and B == 'Mo':
        surfs.append(df.loc['ohooh', ['E', '#H', '#O', '#OH', '#OOH']].values)
        surfs.append(df.loc['oohoh', ['E', '#H', '#O', '#OH', '#OOH']].values)
        surfs.append(df.loc['oooh', ['E', '#H', '#O', '#OH', '#OOH']].values)
        surfs.append(df.loc['ooho', ['E', '#H', '#O', '#OH', '#OOH']].values)
        surfs.append(df.loc['oohooh', ['E', '#H', '#O', '#OH', '#OOH']].values)
        
    nsurfs = len(surfs)
    lowest_surfaces = []
    
    for j in U2:
        values = [dg(k, 0, j) for k in range(nsurfs) if dg(k, 0, j) is not None]
        lowest_surfaces.append(np.argmin(values))
    
    crossover = []
    uniquesurf = [lowest_surfaces[0]]
    old_value = lowest_surfaces[0]
    crossover.append(Umin)
    
    for j in range(len(U2)):
        if lowest_surfaces[j] != old_value:
            uniquesurf.append(lowest_surfaces[j])
            crossover.append(U2[j])
            old_value = lowest_surfaces[j]
    
    crossover.append(Umax2)
    
    plt.clf()
    fig = plt.figure(figsize=fig_size, dpi=300)
    ax = fig.add_axes([0.2, 0.2, 0.6, 0.6])
    ax.axis([0, 14, Umin, Umax])
    ax.set_xlabel(r'pH', fontsize='large')
    ax.set_ylabel(r'U/V', fontsize='large')
    current_yticks = list(plt.yticks()[0])  # Get the current y-ticks
    extraticks = [1.23]
    combined_ticks = sorted(set(current_yticks) | set(extraticks))
    plt.yticks(combined_ticks)
    plt.gca().yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

    for i in range(len(uniquesurf)):
        k = uniquesurf[i]
        label = r"S$_{%i}$(H-%i O-%i OH-%i OOH-%i)" % (k, surfs[k][1], surfs[k][2], surfs[k][3], surfs[k][4])
        plt.fill_between(pH2, crossover[i] - pH2 * const, crossover[i + 1] - pH2 * const, 
                         facecolor=color[k], alpha=0.3, lw=0.5, edgecolor='black')
        plt.plot([], [], color=color[k], alpha=0.3, linewidth=5, label=label)

    plt.plot(pH2, 1.23 - pH2 * const, '--', color='blue', lw=1, dashes=(3, 1))
    plt.plot(pH2, onsetpotential_oer['overP'][0] - pH2 * const, '--', color='black', lw=1, dashes=(3, 1))
    if A == '1' and B == 'Fe':
        plt.plot(pH2, overpotentials['overP'][1] - pH2 * const, '--', color='darkorange', lw=1, dashes=(3, 1))
        plt.plot(pH2, overpotentials['overP'][2] - pH2 * const, '--', color='lime', lw=1, dashes=(3, 1))
        plt.plot(pH2, overpotentials['overP'][3] - pH2 * const, '--', color='cyan', lw=1, dashes=(3, 1))
        ax.text(0.2, 0.65, r'2H$_2$O $\leftrightarrow$ 4H$^+$ + O$_2$ + 4e$^-$', color='blue', rotation=-9.5, fontsize=10)
        ax.text(7.7, overpotentials['overP'][1] - 0.96, f"$S_8$ (*OH+*O): {overpotentials['overP'][1]:.2f} eV", color='darkorange', rotation=-9.5, fontsize=10)
        ax.text(7.7, overpotentials['overP'][2] - 0.65, f"$S_9$ (*OOH): {overpotentials['overP'][2]:.2f} eV", color='lime', rotation=-9.5, fontsize=10)
        ax.text(7.7, overpotentials['overP'][3] - 0.65, f"$S_9$ (*OOH): {overpotentials['overP'][3]:.2f} eV", color='cyan', rotation=-9.5, fontsize=10)
    elif A == '2' and B == 'Co':
        plt.plot(pH2, overpotentials['overP'][1] - pH2 * const, '--', color='darkorange', lw=1, dashes=(3, 1))
        plt.plot(pH2, overpotentials['overP'][2] - pH2 * const, '--', color='lime', lw=1, dashes=(3, 1))
        ax.text(0.2, 0.88, r'2H$_2$O $\leftrightarrow$ 4H$^+$ + O$_2$ + 4e$^-$', color='blue', rotation=-9.5, fontsize=10)
        ax.text(7.7, overpotentials['overP'][1] - 0.71, f"$S_8$ (*OH+*O): {overpotentials['overP'][1]:.2f} eV", color='darkorange', rotation=-9.5, fontsize=10)
        ax.text(7.7, overpotentials['overP'][2] - 0.95, f"$S_9$ (*OOH): {overpotentials['overP'][2]:.2f} eV", color='lime', rotation=-9.5, fontsize=10)
    elif A == '3' and B == 'Mo':
        plt.plot(pH2, overpotentials['overP'][1] - pH2 * const, '--', color='darkorange', lw=1, dashes=(3, 1))
        plt.plot(pH2, overpotentials['overP'][2] - pH2 * const, '--', color='lime', lw=1, dashes=(3, 1))
        ax.text(0.2, 0.88, r'2H$_2$O $\leftrightarrow$ 4H$^+$ + O$_2$ + 4e$^-$', color='blue', rotation=-9.5, fontsize=10)
        ax.text(7.7, overpotentials['overP'][1] - 0.71, f"$S_8$ (*OH+*O): {overpotentials['overP'][1]:.2f} eV", color='darkorange', rotation=-9.5, fontsize=10)
        ax.text(7.7, overpotentials['overP'][2] - 0.65, f"$S_9$ (*OOH): {overpotentials['overP'][2]:.2f} eV", color='lime', rotation=-9.5, fontsize=10)
    
    plt.legend(loc='lower left', bbox_to_anchor=(0.0, 1.02), # borderaxespad=17, 
               ncol=1, labelspacing=0.3, handlelength=2, fontsize=10,
               fancybox=True, shadow=True)
    plt.savefig(f'/pscratch/sd/j/jiuy97/6_MNC/figures/{A}{B}_pourbaix_full.png', bbox_inches='tight')
    print(f"Figure saved as {A}{B}_pourbaix_full.png")
    plt.close()
    
    plt.clf()
    fig = plt.figure(figsize=fig_size, dpi=300)
    ax = fig.add_axes([0.2, 0.2, 0.6, 0.6])
    if A=='3' and B=='Mo':
        ax.axis([-1.0, 2.5, -900, 300])
    else:
        ax.axis([-1.0, 2.5, -600, 200])
    ax.set_xlabel(r'RHE (V)', fontsize='large')
    ax.set_ylabel(r'$\Delta$G (kJ/mol)', fontsize='large')
    xx = np.arange(-1.00, 2.55, 0.05)
    for k in range(nsurfs):
        label = r"S$_{%i}$(H: %i O: %i OH: %i OOH: %i)" % (k, surfs[k][1], surfs[k][2], surfs[k][3], surfs[k][4])
        dg_value = dg(k, 0, xx)
        if dg_value is not None:
            ax.plot(xx, dg_value * kjmol, '-', lw=1, c=color[k], label=label)
        else:
            print(f"Skipping plot for surface {k} due to missing data.")
    plt.xlim(-1.0, 2.5)
    plt.legend(loc='upper left', bbox_to_anchor=(1.0, 1.0), # borderaxespad=17, 
               ncol=1, labelspacing=0.3, handlelength=2, fontsize=10,
               fancybox=True, shadow=True)
    # plt.legend(loc='lower left', bbox_to_anchor=(0.0, 1.02), # borderaxespad=17, 
    #            ncol=2, columnspacing=1.0, labelspacing=0.3, handlelength=2, fontsize=10,
    #            fancybox=True, shadow=True)
    plt.savefig(f'/pscratch/sd/j/jiuy97/6_MNC/figures/{A}{B}_pourbaix.png', bbox_inches='tight')
    print(f"Figure saved as {A}{B}_pourbaix.png")
    # plt.show()
    plt.close()
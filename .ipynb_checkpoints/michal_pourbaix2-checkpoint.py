import os
import re
import numpy as np
import matplotlib.pyplot as plt

# List of main directories
main_dirs = ["clean", "h", "o", "oh", "oh-o", "oho", "oh-oh", "ohoh", "o-o", "oo", "o-oh", "ooh"]

# Regular expression to match E0 values in scientific notation
e0_pattern = re.compile(r"E0=\s*(-?\.\d+E[+-]?\d+)")

# Function to find the lowest E0 value in each subdirectory
def find_min_e0(main_dir, sub_dirs):
    min_e0_values = {}
    for sub_dir in sub_dirs:
        oszicar_path = os.path.join(main_dir, sub_dir, "OSZICAR")
        if os.path.isfile(oszicar_path):
            min_e0 = None
            with open(oszicar_path, 'r') as file:
                for line in file:
                    match = e0_pattern.search(line)
                    if match:
                        e0_value = float(match.group(1))
                        if min_e0 is None or e0_value < min_e0:
                            min_e0 = e0_value
            if min_e0 is not None:
                min_e0_values[sub_dir] = min_e0
    return min_e0_values

# Constants and variables for plotting
Umin, Umax = -0.5, 2.5
kbt = 0.0256 
const = kbt * np.log(10)
kjmol = 96.485

pH = np.arange(0, 14, 0.10)
U = np.arange(Umin, Umax, 0.05)
Umax2 = Umax + 0.06 * 14
U2 = np.arange(Umin, Umax2, 0.05)

color = ['turquoise', 'green', 'red', 'blue', 'gray', 'gold', 'purple', 'pink', 'darkorange',
         'lime', 'olive', 'yellowgreen', 'violet', 'navy', 'brown', 'teal', 'deeppink',
         'cyan', 'dodgerblue', 'steelblue', 'darkslategrey']
pH2 = np.arange(0, 14.01, 0.01)

def addO(x, y):
    return -(h2o - h2) - 2 * (y + x * const) + dso

def addOH(x, y):
    return -(h2o - 0.5 * h2) - (y + x * const) + dsoh

def addOOH(x, y):
    return -(2 * h2o - 1.5 * h2) - 3 * (y + x * const) + dsooh

def addH(x, y):
    return -0.5 * h2 + 1 * (y + x * const) + dsh

def dg(i, x, y):
    return surfs[i][0] - surfs[0][0] + surfs[i][1] * addH(x, y) + surfs[i][2] * addO(x, y) + surfs[i][3] * addOH(x, y) + surfs[i][4] * addOOH(x, y)

# Iterate through each main directory to extract E0 values and plot
for main_dir in main_dirs:
    min_e0_values = find_min_e0(main_dir, ["HS1", "HS5", "IS1", "IS5", "LS1", "LS5"])
    print(main_dir, min_e0_values)

    if len(min_e0_values) < 3:
        print(f"Not enough data in directory '{main_dir}' for plotting.")
        continue

    G_clean = min_e0_values.get("clean", None)
    G_O = min_e0_values.get("o", None)
    G_OH = min_e0_values.get("oh", None)

    if G_clean is None or G_O is None or G_OH is None:
        print(f"Missing data in directory '{main_dir}' for plotting.")
        continue

    # Define surfaces with extracted E0 values
    surfs = [
        [G_clean, 0, 0, 0, 0],  # [energy, #Hs, #Os, #OHs, #OOHs]
        [G_O, 0, 1, 0, 0],
        [G_OH, 0, 0, 1, 0]
    ]
    nsurfs = len(surfs)
    lowest_surfaces = []

    for j in U2:
        values = [dg(k, 0, j) for k in range(nsurfs)]
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
    extraticks = [1.23]
    plt.yticks(list(plt.yticks()[0]) + extraticks)

    for i in range(len(uniquesurf)):
        k = uniquesurf[i]
        label = r"S$_{%i}$(H-%i O-%i OH-%i OOH-%i)" % (k, surfs[k][1], surfs[k][2], surfs[k][3], surfs[k][4])
        plt.fill_between(pH2, crossover[i] - pH2 * const, crossover[i + 1] - pH2 * const, 
                         facecolor=color[k], alpha=0.3, lw=0.5, edgecolor='black')
        plt.plot([], [], color=color[k], alpha=0.3, linewidth=5, label=label)

    plt.plot(pH2, 1.23 - pH2 * const, '--', color='blue', lw=1, dashes=(3, 1))
    ax.text(0.2, 1.00, r'2H$_2$O $\leftrightarrow$ 4H$^+$ + O$_2$ + 4e$^-$', color='blue', rotation=-9, fontsize=10)
    plt.legend(bbox_to_anchor=(0.35, 1.4), loc=2, borderaxespad=0.5, ncol=1, fancybox=True, shadow=True, fontsize=10, handlelength=2)
    plt.savefig(f'pourbaix_full_{main_dir}.png', bbox_inches='tight')
    print(f"Figure saved as pourbaix_full_{main_dir}.png")
    plt.close()
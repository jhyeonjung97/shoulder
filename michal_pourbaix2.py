import os
import re
import numpy as np
import pandas as pd
from matplotlib import rc
import matplotlib.pyplot as plt

# List of main directories
main_dirs = ["clean", "h", "o", "oh", "oh-o", "oho", "oh-oh", "ohoh", "o-o", "oo", "o-oh", "ooh"]

# Regular expression to match E0 values in scientific notation
e0_pattern = re.compile(r"E0=\s*(-?\.\d+E[+-]?\d+)")

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

h2 = -6.77149190
h2o = -14.23091949

zpeh2o = 0.560
zpeh2 = 0.268
cvh2o = 0.103
cvh2 = 0.0905
tsh2o = 0.675
tsh2 = 0.408

dgh2o = zpeh2o + cvh2o - tsh2o
dgh2 = zpeh2 + cvh2 - tsh2

dso = 0.064 + 0.034 - 0.060 - (dgh2o - dgh2)
dsoh = 0.376 + 0.042 - 0.066 - (dgh2o - 0.5 * dgh2)
dsooh = 0.471 + 0.077 - 0.134 - (2 * dgh2o - 1.5 * dgh2)
dsh = dsoh - dso

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
    if surfs[i][0] is None:
        return Nonef
    return (surfs[i][0] 
            - surfs[0][0] 
            + surfs[i][1] * addH(x, y) 
            + surfs[i][2] * addO(x, y) 
            + surfs[i][3] * addOH(x, y) 
            + surfs[i][4] * addOOH(x, y))
    
min_e0_values = {}
# Iterate through each main directory to extract E0 values and plot
for main_dir in main_dirs:
    min_e0 = find_min_e0(main_dir, ["HS1", "HS5", "IS1", "IS5", "LS1", "LS5"])
    print(main_dir, min_e0)

    if min_e0 is None:
        print(f"Missing data in directory '{main_dir}' for plotting.")
        continue
    else:
        min_e0_values[main_dir] = min_e0

G_clean = min_e0_values.get("clean", None)
G_H = min_e0_values.get("h", None)
G_O = min_e0_values.get("o", None)
G_OH = min_e0_values.get("oh", None)
G_OH_O = min_e0_values.get("oh-o", None)
G_O_OH = min_e0_values.get("o-oh", None)
G_OH_OH = min_e0_values.get("oh-oh", None)
G_OHOH = min_e0_values.get("ohoh", None)
# G_O_O = min_e0_values.get("o-o", None)
# G_OO = min_e0_values.get("oo", None)
G_OOH = min_e0_values.get("ooh", None)

# Define surfaces with extracted E0 values
surfs = [
    [G_clean, 0, 0, 0, 0],  # [energy, #Hs, #Os, #OHs, #OOHs]
    [G_H, 1, 0, 0, 0],
    [G_O, 0, 1, 0, 0],
    [G_OH, 0, 0, 1, 0],
    [G_OH_O, 0, 1, 1, 0],
    [G_O_OH, 0, 1, 1, 0],
    [G_OH_OH, 0, 0, 2, 0],
    [G_OHOH, 0, 0, 2, 0],
    # [G_O_O, 0, 2, 0, 0],
    # [G_OO, 0, 2, 0, 0],
    [G_O_OH, 0, 1, 1, 0],
    [G_OOH, 0, 0, 0, 1],
]
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
extraticks = [1.23]
plt.yticks(list(plt.yticks()[0]) + extraticks)

for i in range(len(uniquesurf)):
    k = uniquesurf[i]
    label = r"S$_{%i}$(H-%i O-%i OH-%i OOH-%i)" % (k, surfs[k][1], surfs[k][2], surfs[k][3], surfs[k][4])
    plt.fill_between(pH2, crossover[i] - pH2 * const, crossover[i + 1] - pH2 * const, 
                     facecolor=color[k], alpha=0.3, lw=0.5, edgecolor='black')
    plt.plot([], [], color=color[k], alpha=0.3, linewidth=5, label=label)

plt.plot(pH2, 1.23 - pH2 * const, '--', color='blue', lw=1, dashes=(3, 1))
ax.text(0.2, 0.6, r'2H$_2$O $\leftrightarrow$ 4H$^+$ + O$_2$ + 4e$^-$', color='blue', rotation=-9, fontsize=10)
plt.legend(loc='lower left', bbox_to_anchor=(0.0, 1.02), # borderaxespad=17, 
           ncol=1, labelspacing=0.3, handlelength=2, fontsize=10,
           fancybox=True, shadow=True)
plt.savefig(f'pourbaix_full.png', bbox_inches='tight')
print(f"Figure saved as pourbaix_full.png")
plt.close()

plt.clf()
fig = plt.figure(figsize=fig_size, dpi=300)
ax = fig.add_axes([0.2, 0.2, 0.6, 0.6])
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
plt.legend(loc='lower left', bbox_to_anchor=(0.0, 1.02), # borderaxespad=17, 
           ncol=2, columnspacing=1.0, labelspacing=0.3, handlelength=2, fontsize=10,
           fancybox=True, shadow=True)
plt.savefig(f'pourbaix.png', bbox_inches='tight')
print(f"Figure saved as pourbaix.png")
# plt.show()
plt.close()

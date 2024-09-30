#!/usr/bin/env python

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rc

# Set figure size and font properties
fig_width_pt = 1.8 * 246.0  # LaTeX column width
inches_per_pt = 1.0 / 72.27  # Convert pt to inches
golden_mean = (np.sqrt(5) - 1.0) / 2.0  # Aesthetic ratio
fig_width = fig_width_pt * inches_per_pt  # Width in inches
fig_height = fig_width * golden_mean  # Height in inches
fig_size = [fig_width, fig_height]

font_size = 10
tick_font_size = 10

# Matplotlib configuration
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['font.size'] = font_size
plt.rcParams['axes.labelsize'] = 2 * font_size
plt.rcParams['legend.fontsize'] = font_size
plt.rcParams['xtick.labelsize'] = tick_font_size
plt.rcParams['ytick.labelsize'] = tick_font_size
plt.rcParams['mathtext.default'] = 'regular'
plt.rcParams['lines.linewidth'] = 1.0

# Define constants
Umin, Umax = -0.5, 2.5
kbt = 0.0256 
const = kbt * np.log(10)
kjmol = 96.485

# Define pH and U ranges
pH = np.arange(0, 14, 0.10)
U = np.arange(Umin, Umax, 0.05)
Umax2 = Umax + 0.06 * 14
U2 = np.arange(Umin, Umax2, 0.05)

# Read data from scaling_relationship.tsv
data = pd.read_csv('/pscratch/sd/j/jiuy97/6_MNC/figure/scaling_relationship.tsv', sep='\t', header=0, index_col=0)

# Loop through each metal
for metal in data.index:
    # Extract values for calculations
    G_clean = data.loc[metal, 'G_']  # Clean surface
    G_O = data.loc[metal, 'G_O']  # 1*O surface
    G_OH = data.loc[metal, 'G_OH']  # 1*OH surface

    # Define surface data for the current metal
    surfs = [
        [G_clean, 0, 0, 0, 0],  # [energy, #Hs, #Os, #OHs, #OOHs]
        [G_O, 0, 1, 0, 0],
        [G_OH, 0, 0, 1, 0]
    ]
    
    # Find surface energy intersections
    nsurfs = len(surfs)
    lowest_surfaces = []
    for j in U2:
        values = [dg(k, 0, j) for k in range(nsurfs)]
        lowest_surfaces.append(np.argmin(values))
    
    # Identify crossover points
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
    
    # Define colors for different surfaces
    color = ['turquoise', 'green', 'red', 'blue', 'gray', 'gold', 'purple', 'pink', 'darkorange',
             'lime', 'olive', 'yellowgreen', 'violet', 'navy', 'brown', 'teal', 'deeppink',
             'cyan', 'dodgerblue', 'steelblue', 'darkslategrey']
    pH2 = np.arange(0, 14, 0.01)
    
    # Plot surface Gibbs energies
    for i in range(len(uniquesurf)):
        k = uniquesurf[i]
        label = r"S$_{%i}$(H-%i O-%i OH-%i OOH-%i)" % (k, surfs[k][1], surfs[k][2], surfs[k][3], surfs[k][4])
        fill_between(pH2, crossover[i] - pH2 * const, crossover[i + 1] - pH2 * const, facecolor=color[i], alpha=0.3, lw=0.5, edgecolor='black')
        plt.plot([], [], color=color[i], alpha=0.3, linewidth=5, label=label)
    
    # Plot OER line
    Vover = 0.184
    y = 1.23 + Vover - pH2 * const
    plt.plot(pH2, y, '-', color='black', lw=1, dashes=(3, 1), label='$\eta$ OER = ' + repr(Vover) + ' ')
    
    # Plot equilibrium line
    plt.plot(pH2, 1.23 - pH2 * const, '--', color='blue', lw=1, dashes=(3, 1))
    ax.text(0.2, 1.25, r'2H$_2$O $\leftrightarrow$ 4H$^+$ + O$_2$ + 4e$^-$', color='blue', rotation=-7, fontsize='x-small')
    
    # Add legend and save the plot
    plt.legend(bbox_to_anchor=(0.05, 1), loc=2, borderaxespad=0., ncol=2, fancybox=True, shadow=True, fontsize='x-small', handlelength=3)
    plt.savefig(f'pourbaix_{metal}.pdf', bbox_inches='tight')
    
    # Clear the plot and set up for the second figure
    plt.clf()
    fig = plt.figure(figsize=fig_size, dpi=300)
    ax = fig.add_axes([0.2, 0.2, 0.6, 0.6])
    ax.axis([-1.0, 2.5, -800, 200])
    ax.set_xlabel(r'RHE (V)')
    ax.set_ylabel(r'$\Delta$G (kJ/mol)')
    
    # Plot Gibbs energies for surfaces in the second figure
    xx = np.arange(-1.0, 2.5, 0.05)
    for k in range(nsurfs):
        label = r"S$_{%i}$(H: %i O: %i OH: %i OOH: %i)" % (k, surfs[k][1], surfs[k][2], surfs[k][3], surfs[k][4])
        ax.plot(xx, dg(k, 0, xx) * kjmol, '-', lw=1, c=color[k], label=label)
    
    # Add legend and save the second plot
    plt.legend(bbox_to_anchor=(0.05, 1.3), loc=2, borderaxespad=0., ncol=2, fancybox=True, shadow=True, fontsize='x-small', handlelength=3)
    plt.savefig(f'pourbaix2_{metal}.pdf', bbox_inches='tight')

# Show the final plot
plt.show()

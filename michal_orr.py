#! /user/bin/env python
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.markers import MarkerStyle
import csv

# Figure and font settings
fig_width_pt = 1.8 * 246.0
inches_per_pt = 1.0 / 72.27
golden_mean = (np.sqrt(5) - 1.0) / 2.0
fig_width = fig_width_pt * inches_per_pt
fig_height = fig_width * golden_mean
fig_size = [fig_width, fig_height]
fig = plt.figure(figsize=fig_size, dpi=300)

font_size = 9
tick_font_size = 8
plt.rcParams.update({
    'ps.usedistiller': 'xpdf',
    'font.size': font_size,
    'axes.labelsize': font_size,
    'legend.fontsize': font_size,
    'xtick.labelsize': tick_font_size,
    'ytick.labelsize': tick_font_size,
    'lines.linewidth': 1.0
})

def setfont(font='cmss', unicode=True):
    """
    Set Matplotlib rcParams to use LaTeX for font rendering.
    """
    font = font.lower().replace(" ", "")
    font = {'family': 'sans-serif', 'serif': ['cmss']}
    preamble = r"""\usepackage{color} \usepackage[tx]{sfmath}"""
    plt.rc('font', **font)
    plt.rcParams['text.latex.preamble'] = preamble

setfont()

# Plot settings
ax = fig.add_axes([0.2, 0.2, 0.6, 0.6])
zoomx, zoomy = 0.4, 0.6
xcenter, ycenter = 0.9, 3.75
d1, d2, d3 = 4 * zoomx, 3 * zoomy, 5 * zoomx
x1, x2 = xcenter - d1, xcenter + d3
y1, y2 = ycenter - d2, ycenter + d2

ax.axis([x1, x2, y1, y2])
ax.set_xlabel(r'$\Delta$G$_{\sf OH}$ (eV)')
ax.set_ylabel(r'$\Delta$G$_{\sf OOH}$ (eV)')

# Define functions for overpotential calculations
def ooh_oh_scaling(doh):
    return doh + 3.2

def overpotential_orr_for_contour(doh, dooh):
    do = 1.469 * doh + 1.253
    dg14 = [-doh, -do + doh, -dooh + do, -4.92 + dooh]
    return max(dg14) + 1.23

def overpotential_orr(doh, do, dooh):
    dg14 = [-doh, -do + doh, -dooh + do, -4.92 + dooh]
    return max(dg14) + 1.23

def orr_step(i):
    steps = ['O2->OOH*', 'OOH*->O*', 'O*->OH*', 'OH*->H2O']
    return steps[i]

def overpotential_orr_full(doh, do, dooh):
    dg14 = [-4.92 + dooh, -dooh + do, -do + doh, -doh]
    m = max(dg14)
    return [round(m + 1.23, 2), round(-m, 2), orr_step(dg14.index(m))]

# Generate data for contour plot
delta = 0.01
x = np.arange(x1, x2 + delta, delta)
y = np.arange(y1, y2 + delta, delta)
X, Y = np.meshgrid(x, y)

Z = np.array([[overpotential_orr_for_contour(i, j) for i in x] for j in y])

# Plot contour
levels = np.arange(0.1, 1.7, 0.1)
CS = plt.contourf(X, Y, Z, levels, cmap=ListedColormap([
    '#a50026', '#d73027', '#f46d43', '#fdae61', '#fee090', '#ffffbf',
    '#ffffe5', '#ffffff', '#e0f3f8', '#abd9e9', '#74add1', '#4575b4', '#313695'
]), extend='max', origin='lower')

cbar = plt.colorbar(CS, ticks=[np.arange(0.1, 1.6, step=0.1)])
cbar.ax.set_ylabel(r'$\eta_{\sf ORR}$ (V)')
cbar.ax.tick_params(size=3, labelsize=6, labelcolor='black', width=0.5, color='black')

# Calculate and plot overpotentials for each system
calc_systems = [
    [1.456, 2.892, 4.432, 1.752, r'CoCo_M1', '#737373', 0.0, 0.08, 1.5, '8', 'blue', 'o', 'red'],
    # More systems...
]

for i, system in enumerate(calc_systems):
    if i % 2 == 0:
        ax.plot(system[0], system[2], system[9], mec=system[5], mfc=system[11], mew=0.8, zorder=4,
                marker=MarkerStyle('o', fillstyle='left'), markersize=8,
                label=f'{system[4]} : {overpotential_orr(system[0], system[1], system[2]):.2f} V')
    else:
        ax.plot(system[0], system[2], system[9], mec=system[5], mfc=system[11], mew=0.8, zorder=4,
                marker=MarkerStyle('o', fillstyle='right'), markersize=8,
                label=f'{system[4]} : {overpotential_orr(system[0], system[1], system[2]):.2f} V')

ax.plot(x, 0.87 * x + 3.22, '--', lw=1, dashes=(3, 1), c='black')
ax.text(1.2, 2.5, r'$\Delta$G$_{\sf OOH}$=0.87$\Delta$G$_{\sf OH}$', color='black', fontsize=10)
ax.text(1.8, 2.34, '+3.22 eV', color='black', fontsize=10)

ax.legend(bbox_to_anchor=(-0.15, 1.65), loc=2, borderaxespad=0.5, ncol=3, fancybox=True, shadow=False, fontsize='x-small', handlelength=2)
fig.savefig('ORR_contour_plot_v13_full_scaling_trilayer_v3.pdf', bbox_inches='tight')
fig.clf()

# CSV writing for overpotential results
with open('Final_dGs_for_MOOHx_system_ORR.csv', 'w', newline='') as myfile:
    fieldnames = ['Surface name', 'dOH', 'dO', 'dOOH', 'overpotential', 'onset potential', 'PLS']
    writer = csv.DictWriter(myfile, fieldnames=fieldnames)
    writer.writeheader()

    solv_corr_OH, solv_corr_O, solv_corr_OOH = -0.116, -0.083, -0.327
    for system in calc_systems:
        system[0] += solv_corr_OH
        system[1] += solv_corr_O
        system[2] += solv_corr_OOH
        recalculated_over = overpotential_orr_full(system[0], system[1], system[2])
        writer.writerow({
            'Surface name': system[4], 'dOH': system[0], 'dO': system[1],
            'dOOH': system[2], 'overpotential': recalculated_over[0],
            'onset potential': recalculated_over[1], 'PLS': recalculated_over[2]
        })
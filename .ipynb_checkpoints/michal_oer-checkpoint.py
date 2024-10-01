import csv
import matplotlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.markers import MarkerStyle

rows = ['3d', '3d', '3d', '3d', '4d', '5d']
groups = ['5', '6', '7', '8', '4', '4']
metals = ['Mn', 'Fe', 'Co', 'Ni', 'Mo', 'W']

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
xcenter, ycenter = 1.60, 0.81
x1, x2 = xcenter - 1.2, xcenter + 1.2 # 3.2
y1, y2 = ycenter - 2.0, ycenter + 2.0 # 3.6

ax.axis([x1, x2, y1, y2])
ax.set_xlabel(r'$\Delta$G$_{\sf O}$ - $\Delta$G$_{\sf OH}$(eV)', fontsize=10)
ax.set_ylabel(r'$\Delta$G$_{\sf OH}$ (eV)', fontsize=10)

# Define functions for overpotential calculations
def ooh_oh_scaling(doh):
    return doh + 3.2

def oer_step(i):
    steps = ['H2O->OH*', 'OH*->O*', 'O*->OOH*', 'OOH*->O2']
    return steps[i]

def overpotential_oer(doh, do, dooh):
    dg14 = [doh, do - doh, dooh - do, 4.92 - dooh]
    return max(dg14) - 1.23

def overpotential_oer_full(doh, do, dooh):
    dg14 = [doh, do - doh, dooh - do, 4.92 - dooh]
    m = max(dg14)
    return [round(m - 1.23, 2), round(-m, 2), oer_step(dg14.index(m))]
    
# def overpotential_oer_for_contour(do_doh, doh):
#     do = do_doh + doh
#     dooh = ooh_oh_scaling(doh)
#     dg14 = [doh, do - doh, dooh - do, 4.92 - dooh]
#     return max(dg14) - 1.23
    
def overpotential_oer_for_contour(x, y):
    dg14 = [y, x, 3.2-x, 1.72-y]
    return max(dg14) - 1.23
    
# Read data from the TSV file
df = pd.read_csv('/pscratch/sd/j/jiuy97/6_MNC/figure/scaling_relationship.tsv', sep='\t', index_col=0)

# Extract values from the dataframe
doh_values = df['dG_OH']
do_values = df['dG_O']

# Add `dG_OOH` and calculate overpotential for each entry
df['dG_OOH'] = doh_values.apply(ooh_oh_scaling)
df['overpotential'] = df.apply(lambda row: overpotential_oer(row['dG_OH'], row['dG_O'], row['dG_OOH']), axis=1)

# Prepare separate data for each metal
dfs = {}
for m, metal in enumerate(metals):
    row = rows[m]
    group = groups[m]
    dfs[metal] = pd.read_csv(f'/pscratch/sd/j/jiuy97/6_MNC/figure/{row}_{group}{metal}_gibbs.tsv', sep='\t', header=0, index_col=0)
    doh_values = dfs[metal]['dG_OH']
    do_values = dfs[metal]['dG_O']
    dfs[metal]['dG_OOH'] = doh_values.apply(ooh_oh_scaling)
    dfs[metal]['overpotential'] = dfs[metal].apply(lambda row: overpotential_oer(row['dG_OH'], row['dG_O'], row['dG_OOH']), axis=1)

# Generate data for contour plot
delta = 0.01
x = np.arange(x1, x2 + delta, delta)
y = np.arange(y1, y2 + delta, delta)
X, Y = np.meshgrid(x, y)
Z = np.array([[overpotential_oer_for_contour(i, j) for i in x] for j in y])

# Plot contour
levels = np.arange(0.3, 1.6, 0.1)
CS = plt.contourf(X, Y, Z, levels, cmap=ListedColormap([
    '#a50026', '#d73027', '#f46d43', '#fdae61', '#fee090', '#ffffbf',
    '#ffffe5', '#ffffff', '#e0f3f8', '#abd9e9', '#74add1', '#4575b4', '#313695'
]), extend='max', origin='lower')

cbar = plt.colorbar(CS, ticks=np.arange(0.3, 1.6, 0.1))
cbar.ax.set_ylabel(r'$\eta_{\sf OER}$ (V)')
cbar.ax.tick_params(size=3, labelsize=6, labelcolor='black', width=0.5, color='black')

# Plot data points from the TSV file with their calculated overpotentials
markers = ['o', 's', 'd', '^', 'v', '*']  # Different markers for metals
colors = ['blue', 'green', 'orange', 'red', 'purple', 'grey']
color_ranges = [
    plt.cm.Blues(np.linspace(0.3, 0.9, 7)),
    plt.cm.Greens(np.linspace(0.3, 0.9, 7)),
    plt.cm.Oranges(np.linspace(0.3, 0.9, 7)),
    plt.cm.Reds(np.linspace(0.3, 0.9, 7)),
    plt.cm.Purples(np.linspace(0.3, 0.9, 7)),
    plt.cm.Greys(np.linspace(0.3, 0.9, 7)),
    ]

# Plot the general dataset points
for row_num, row in enumerate(df.itertuples(), 1):  # Start row number from 1
    ax.scatter(row.dG_O - row.dG_OH, row.dG_OH, 
               label=f'{row.Index}: {row.overpotential:.2f} V',               
               s = 24, marker='x', 
               # marker=markers[row_num-1],
               linewidths=1.0,
               # facecolors='white',
               # edgecolors=colors[row_num-1],
               color=colors[row_num-1],
               zorder=10)

# Plot the metal-specific data points with colormaps
for m, metal in enumerate(metals):
    for row_num, row in enumerate(dfs[metal].itertuples(), 1):  # Use row number here as well
        ax.scatter(row.dG_O - row.dG_OH, row.dG_OH, 
                   s=24, marker='o', 
                   # marker=markers[m],
                   # linewidths=0.5,
                   facecolors=color_ranges[m][row_num-1],
                   edgecolors='none',
                   zorder=9)

# Add scaling line
# ax.plot(x, x+3.2, '--', lw=1, dashes=(3, 1), c='black')
# ax.text(-0.5, 5.3, r'$\Delta$G$_{\sf OOH}$=', color='black', fontsize=10)
# ax.text(-0.5, 5.1, r'$\Delta$G$_{\sf OH}$+3.2 eV', color='black', fontsize=10)

ax.legend(bbox_to_anchor=(0.5, 1.1), loc='center', borderaxespad=0.0, ncol=3, fancybox=True, shadow=False, fontsize='x-small', handlelength=2)
fig.savefig('contour_OER.png', bbox_inches='tight')
print("Figure saved as contour_OER.png")
fig.clf()

# CSV writing for overpotential results
with open('contour_OER.csv', 'w', newline='') as myfile:
    fieldnames = ['Surface name', 'dOH', 'dO', 'dOOH', 'overpotential', 'onset potential', 'PLS']
    writer = csv.DictWriter(myfile, fieldnames=fieldnames)
    writer.writeheader()
    for idx, row in df.iterrows():
        recalculated_over = overpotential_oer_full(row['dG_OH'], row['dG_O'], row['dG_OOH'])
        writer.writerow({
            'Surface name': row.name, 
            'dOH': row['dG_OH'], 'dO': row['dG_O'], 'dOOH': row['dG_OOH'], 
            'overpotential': recalculated_over[0],
            'onset potential': recalculated_over[1], 
            'PLS': recalculated_over[2]
        })

# TSV writing for overpotential results
with open('contour_OER.tsv', 'w', newline='') as myfile:
    fieldnames = ['Surf.', 'dOH', 'dO', 'dO*', 'diff', 'dOOH', 'overP', 'onsetP', 'PLS']
    writer = csv.DictWriter(myfile, fieldnames=fieldnames, delimiter='\t')  # Change delimiter to '\t'
    writer.writeheader()
    for idx, row in df.iterrows():
        recalculated_over = overpotential_oer_full(row['dG_OH'], row['dG_O'], row['dG_OOH'])
        writer.writerow({
            'Surf.': row.name, 
            'dOH': round(row['dG_OH'], 2),
            'dO': round(row['dG_O'], 2),
            'dO*': round(1.8847*row['dG_OH']+0.7599, 2),
            'diff': round(1.8847*row['dG_OH']+0.7599-row['dG_O'], 2),
            'dOOH': round(row['dG_OOH'], 2),
            'overP': round(recalculated_over[0], 2),
            'onsetP': round(recalculated_over[1], 2),
            'PLS': recalculated_over[2]
        })

# Write results for each metal
for m, metal in enumerate(metals):
    with open(f'contour_OER_{m+1}{metal}.tsv', 'w', newline='') as myfile:
        fieldnames = ['Surf.', 'dOH', 'dO', 'dO*', 'diff', 'dOOH', 'overP', 'onsetP', 'PLS']
        writer = csv.DictWriter(myfile, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        for idx, row in dfs[metal].iterrows():
            recalculated_over = overpotential_oer_full(row['dG_OH'], row['dG_O'], row['dG_OOH'])
            writer.writerow({
                'Surf.': row.name, 
                'dOH': round(row['dG_OH'], 2),
                'dO': round(row['dG_O'], 2),
                'dO*': round(1.8847*row['dG_OH']+0.7599, 2),
                'diff': round(1.8847*row['dG_OH']+0.7599-row['dG_O'], 2),
                'dOOH': round(row['dG_OOH'], 2),
                'overP': round(recalculated_over[0], 2),
                'onsetP': round(recalculated_over[1], 2),
                'PLS': recalculated_over[2]
            })
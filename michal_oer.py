import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.markers import MarkerStyle

# Define functions for calculations
def ooh_oh_scaling(doh):
    return doh + 3.2  # Example scaling relationship for OOH

def overpotential_orr(doh, do, dooh):
    dg14 = [-doh, -do + doh, -dooh + do, -4.92 + dooh]
    return max(dg14) + 1.23

def overpotential_orr_full(doh, do, dooh):
    dg14 = [-4.92 + dooh, -dooh + do, -do + doh, -doh]
    m = max(dg14)
    return [round(m + 1.23, 2), round(-m, 2), orr_step(dg14.index(m))]

# Reading data from 'scaling_relationship.tsv'
df = pd.read_csv('/pscratch/sd/j/jiuy97/6_MNC/figure/scaling_relationship.tsv', sep='\t')

# Extract `doh` and `do` from the file
doh_values = df['dG_OH']
do_values = df['dG_O']

# Calculate `dooh` and `overpotential`
df['dG_OOH'] = doh_values.apply(ooh_oh_scaling)
df['overpotential'] = df.apply(lambda row: overpotential_orr(row['dG_OH'], row['dG_O'], row['dG_OOH']), axis=1)

# Plotting the results
fig, ax = plt.subplots(figsize=(8, 6))

# Contour plot settings
x1, x2 = 0.4, 1.6
y1, y2 = 3.0, 5.0
delta = 0.01
x = np.arange(x1, x2 + delta, delta)
y = np.arange(y1, y2 + delta, delta)
X, Y = np.meshgrid(x, y)

# Generate the Z values based on overpotential
Z = np.array([[overpotential_orr(i, j, ooh_oh_scaling(i)) for i in x] for j in y])

# Create contour plot
levels = np.arange(0.1, 1.7, 0.1)
CS = ax.contourf(X, Y, Z, levels, cmap=ListedColormap([
    '#a50026', '#d73027', '#f46d43', '#fdae61', '#fee090', '#ffffbf',
    '#ffffe5', '#ffffff', '#e0f3f8', '#abd9e9', '#74add1', '#4575b4', '#313695'
]), extend='max', origin='lower')

# Add colorbar
cbar = plt.colorbar(CS, ticks=np.arange(0.1, 1.6, 0.1))
cbar.ax.set_ylabel(r'$\eta_{\sf ORR}$ (V)')
cbar.ax.tick_params(size=3, labelsize=6, labelcolor='black', width=0.5, color='black')

# Plot the calculated points from the dataset
for idx, row in df.iterrows():
    ax.plot(row['dG_OH'], row['dG_OOH'], 'o', label=f'{row["Metal"]}: {row["overpotential"]:.2f} V')

# Add labels and title
ax.set_xlabel(r'$\Delta$G$_{\sf OH}$ (eV)')
ax.set_ylabel(r'$\Delta$G$_{\sf OOH}$ (eV)')
ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)

# Save the plot
plt.tight_layout()
plt.savefig('ORR_contour_with_data_points.pdf')
plt.show()

# Save the results (including calculated overpotentials) back to a CSV file
df.to_csv('scaling_relationship_with_overpotential.csv', index=False)

import numpy as np
from matplotlib import pyplot as plt
import os
import imageio.v2 as imageio

##### INPUTS #####

# Parameter space scanned:
# ref_p = np.arange(230,250,1)
# ref_p = np.arange(237,246,0.1)
ref_p = np.arange(245,246,0.01)

# Number of steps in scan:
iterations = len(ref_p)

##### REFERENCE PARTICLE PLOTTING FUNCTIONS #####

# Define function to plot orbit:
def plot_orbit(x_vals, y_vals, p_val, dir):
    plt.plot(x_vals,y_vals)
    plt.title(f'$p$ = {round(p_val,1)} MeV/c')
    plt.xlabel('$x$ [mm]')
    plt.ylabel('$y$ [mm]')
    plt.savefig(dir+'orbit.png',dpi=300)
    plt.close()

##### MAIN LOOP #####

full_channel_indices = []
for j in range(iterations):

    # Import data:
    dir = f'sim{j+1}/'
    file = f'{dir}AllTracks.txt'
    data = np.loadtxt(file)

    # Values along channel:
    x_vals = []; y_vals = []
    for i in range(data.shape[0]):
        x_vals.append(data[i][0])
        y_vals.append(data[i][1])

    # Plot:
    plt.clf()
    plot_orbit(x_vals, y_vals, ref_p[j], dir)

##### ANIMATIONS #####

# List of sim directories:
out_dirs = [f'sim{i+1}' for i in range(iterations)]

# List of plot files:
paths = []
for i in range(iterations):
    paths.append(out_dirs[i]+'/orbit.png')

# Create animations:
plots = [imageio.imread(img) for img in paths]
frame_duration = 500 # for coarse scans
# frame_duration = 100 # for fine scans
imageio.mimsave('refp_scan.gif', plots, duration=frame_duration, loop=0)
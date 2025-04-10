import numpy as np
from matplotlib import pyplot as plt
import os
import imageio.v2 as imageio

##### INPUTS #####

# Main directory:
main_dir = '/Users/criggall/Documents/muon-cooling/Simplified-HFOFO/sol_current_scan/'
fig_dir = '/Users/criggall/Documents/muon-cooling/Simplified-HFOFO/Figures/'

# Name of parameter scan is over:
param_label = 'solenoid current (amps)'

# Parameter space scanned:
sol_current = np.arange(-80.51,-80.39,0.01)

# Number of steps in scan:
iterations = len(sol_current)

# Specify whether reference particle, beam, or both:
ref_particle = True
# beam = False # <-- Still have to add this option

# RF period (325 MHz frequency):
T = 1/(325*10**6)*10**9 # ns

# Length of period:
len_period = 4.2 # m

##### REFERENCE PARTICLE PLOTTING FUNCTIONS #####

# Define function to plot xy trajectory along z:
def plot_trajectory(x_vals, y_vals, z_vals, param_label, param_val, dir):
    plt.figure(figsize = (12,6))
    plt.clf()
    plt.plot(z_vals,x_vals,color='red',label='x')
    plt.plot(z_vals,y_vals,color='blue',label='y')
    plt.ylim(-4,4)
    plt.title(f'{param_label} = {round(param_val,2)}')
    plt.xlabel('z (m)')
    plt.ylabel('x, y (cm)')
    plt.legend(loc='upper left')
    plt.savefig(dir+'xy_trajectory.png',dpi=300)
    plt.close()

# Define function to plot orbit in xy-plane:
def plot_orbit(x_vals, y_vals, z_vals, param_label, param_val, dir):
    plt.clf()
    plt.plot(x_vals,y_vals)
    plt.title(f'{param_label} = {round(param_val,2)}')
    plt.xlabel('x (m)')
    plt.ylabel('y (cm)')
    plt.savefig(dir+'orbit.png',dpi=300)
    plt.close()

##### BEAM PLOTTING FUNCTIONS #####

# Plot total momentum distributions:
# def plot_ptotal_dist():


##### MAIN LOOP #####

full_channel_indices = []
x_total_residuals = []
y_total_residuals = []
for j in range(iterations):

    if ref_particle == True:

        # Import data:
        dir = f'{main_dir}g4bl-output-sim{j+1}/'
        file = f'{dir}AllTracks.txt'
        data = np.loadtxt(file)

        # Values along channel:
        x_vals = []; y_vals = []; z_vals = []
        px_vals = []; py_vals = []; pz_vals = []; ptotal_vals = []
        t_vals = []; mod_t_vals = []
        Bx_vals = []; By_vals = []; Bz_vals = []
        count = 0; count2 = 0
        next_z = 0; period_start_z_vals = []; period_start_indices = []
        for i in range(data.shape[0]):
            x_vals.append(data[i][0]*0.1) # mm -> cm
            y_vals.append(data[i][1]*0.1)
            z = data[i][2]*0.001 # mm --> m
            z_vals.append(z)
            px = data[i][3]; py = data[i][4]; pz = data[i][5]
            px_vals.append(px) # MeV/c
            py_vals.append(py)
            pz_vals.append(pz)
            ptotal_vals.append(np.sqrt(px**2+py**2+pz**2))
            t = data[i][6]
            t_vals.append(t) # ns
            mod_t = t % T
            mod_t_vals.append(mod_t)
            Bx = data[i][12]*200; By = data[i][13]*200; Bz = data[i][14] # rescale Bx, By for plotting
            Bx_vals.append(Bx)
            By_vals.append(By)
            Bz_vals.append(Bz)
            del px, py, pz, t, Bx, By, Bz

            # Find indices for first value in each period:
            if z > next_z:
                period_start_z_vals.append(z)
                period_start_indices.append(i)
                next_z += len_period

        # Remove last period -- just for reducing end field effects:
        start_last_period = period_start_indices[len(period_start_indices)-1]
        x_vals = x_vals[:start_last_period]
        y_vals = y_vals[:start_last_period]
        z_vals = z_vals[:start_last_period]
        px_vals = px_vals[:start_last_period]
        py_vals = py_vals[:start_last_period]
        pz_vals = pz_vals[:start_last_period]
        Bx_vals = Bx_vals[:start_last_period]
        By_vals = By_vals[:start_last_period]
        Bz_vals = Bz_vals[:start_last_period]

        # Plot:
        # plot_trajectory(x_vals, y_vals, z_vals, param_label, sol_current[j], dir)
        plot_orbit(x_vals, y_vals, z_vals, param_label, sol_current[j], dir)

        del x_vals, y_vals, z_vals

    # if beam == True:


##### ANIMATIONS #####

# # List of sim directories:
# out_dirs = [main_dir+f'g4bl-output-sim{i+1}' for i in range(iterations)]

# # List of plot files:
# xy_plot_paths = []
# for i in range(iterations):
#     xy_plot_paths.append(out_dirs[i]+'/xy_trajectory.png')

# # Create animations:
# xy_images = [imageio.imread(img) for img in xy_plot_paths]
# frame_duration = 50
# imageio.mimsave(fig_dir+'xy_animation_full_channel.gif', xy_images, duration=frame_duration, loop=0)
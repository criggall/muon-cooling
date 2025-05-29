import numpy as np
from matplotlib import pyplot as plt
import os
import imageio.v2 as imageio

##### INPUTS #####

# Main directory:
# main_dir = '/Users/criggall/Documents/muon-cooling/Simplified-HFOFO/delta_p_scan_fine/'
# main_dir = '/Users/criggall/Documents/muon-cooling/Simplified-HFOFO/delta_p_scan_coarse/'
main_dir = '/Users/criggall/Documents/muon-cooling/Simplified-HFOFO/sol_tilt_scan/'

# Figure directory:
fig_dir = '/Users/criggall/Documents/muon-cooling/Simplified-HFOFO/Figures/'

# Name of parameter scan is over:
param_label = 'solenoid tilt'
param_units = '$^{\circ}$'

# Parameter space scanned:
# delta_p = np.arange(-5.0, 5.0, 0.1)
# delta_p = np.arange(-10.0, 10.0, 0.5)
nominal_pitch = -0.0025*180/np.pi
pitch = np.linspace(0,nominal_pitch,10)

# Number of steps in scan:
iterations = len(pitch)

# RF period (325 MHz frequency):
T = 1/(325*10**6)*10**9 # ns

# Length of period:
len_period = 4.2 # m

##### MATCHED REFERENCE PARTICLE DATA #####

# Read in matched reference particle data:
file_ref = '/Users/criggall/Documents/muon-cooling/Simplified-HFOFO/AllTracks_ref.txt'
data_ref = np.loadtxt(file_ref)

# Values along channel:
x_vals_ref = []; y_vals_ref = []; z_vals_ref = []
next_z = 0.5; period_start_z_vals = []; period_start_indices = [] # start solenoid placement at 0.5m
for i in range(data_ref.shape[0]):
    x_vals_ref.append(data_ref[i][0]*0.1) # mm -> cm
    y_vals_ref.append(data_ref[i][1]*0.1)
    z = data_ref[i][2]*0.001 # mm -> m
    z_vals_ref.append(z)

    # Find indices for first value in each period:
    if z > next_z:
        period_start_z_vals.append(z)
        period_start_indices.append(i)
        next_z += len_period

# Remove last period -- just for reducing end field effects:
start_last_period = period_start_indices[len(period_start_indices)-1]
x_vals_ref = x_vals_ref[:start_last_period]
y_vals_ref = y_vals_ref[:start_last_period]
z_vals_ref = z_vals_ref[:start_last_period]

##### PLOTTING FUNCTIONS #####

# Define function to plot xy trajectory along z:
def plot_trajectory(x_vals, y_vals, z_vals, param_label, param_val, dir):
    plt.figure(figsize = (12,6))
    plt.clf()
    plt.plot(z_vals,x_vals,color='red',label='x')
    plt.plot(z_vals,y_vals,color='blue',label='y')
    plt.ylim(-4,4)
    plt.title(f'{param_label} = {round(param_val,2)} {param_units}')
    plt.xlabel('z (m)')
    plt.ylabel('x, y (cm)')
    plt.legend(loc='upper left')
    plt.savefig(dir+'xy_trajectory.png',dpi=300)
    plt.close()

# Define function to plot orbit in xy-plane:
def plot_orbit(x_vals, y_vals, z_vals, param_label, param_val, dir):
    plt.clf()
    plt.plot(x_vals,y_vals)
    plt.title(f'{param_label} = {round(param_val,2)} {param_units}')
    plt.xlabel('x (cm)')
    plt.ylabel('y (cm)')
    plt.xlim(-2.5,2.5)
    plt.ylim(-2.5,2.5)
    plt.savefig(dir+'orbit.png',dpi=300)
    plt.close()

# Define function to plot dispersion along z:
def plot_dispersion(D_vals, z_vals, param_label, param_val, dir):
    plt.clf()
    plt.figure(figsize=(10,4))
    plt.plot(z_vals, D_vals)
    plt.title(f'{param_label} = {round(param_val,2)} {param_units}')
    # plt.ylim(-35,35)
    plt.xlim(0,130)
    plt.xlabel('z (m)')
    plt.ylabel('D (cm)')
    plt.savefig(dir+'dispersion.png',dpi=300)
    plt.close()

# Define function to plot Lz along z:
def plot_angular_momentum(Lz_vals, z_vals, param_label, param_val, dir):
    plt.clf()
    plt.figure(figsize=(10,4))
    plt.plot(z_vals, Lz_vals)
    plt.title(f'{param_label} = {round(param_val,2)} {param_units}')
    plt.ylim(-30,15)
    plt.xlim(0,130)
    plt.xlabel('z (m)')
    plt.ylabel('$L_z$ (cm*MeV/c)')
    plt.savefig(dir+'angular_momentum.png',dpi=300)
    plt.close()

# Define function to plot Bz along z:
# 

# Define function to plot Bx, By along z:
# 

##### MAIN LOOP #####

full_channel_indices = []
x_total_residuals = []
y_total_residuals = []
for j in range(iterations):

    # Import data:
    dir = f'{main_dir}g4bl-output-sim{j+1}/'
    file = f'{dir}AllTracks.txt'
    data = np.loadtxt(file)

    # Values along channel:
    x_vals = []; y_vals = []; z_vals = []
    px_vals = []; py_vals = []
    next_z = 0; period_start_z_vals = []; period_start_indices = []
    Lz_vals = []
    for i in range(data.shape[0]):
        x = data[i][0]*0.1; y = data[i][1]*0.1 # mm -> cm
        x_vals.append(x); y_vals.append(y)
        z = data[i][2]*0.001 # mm --> m
        z_vals.append(z)
        px = data[i][3]; py = data[i][4]
        Lz = x*py - y*px
        Lz_vals.append(Lz)

        # Find indices for first value in each period:
        if z > next_z:
            period_start_z_vals.append(z)
            period_start_indices.append(i)
            next_z += len_period

    # Remove last period -- just for reducing end field effects:
    start_last_period = period_start_indices[len(period_start_indices)-2]
    x_vals = x_vals[:start_last_period]
    y_vals = y_vals[:start_last_period]
    z_vals = z_vals[:start_last_period]
    Lz_vals = Lz_vals[:start_last_period]

    # Compute dispersion:
    D_vals = []
    for i in range(len(z_vals)):
        D_val = (np.sqrt( (x_vals_ref[i] - x_vals[i])**2 + (y_vals_ref[i] - y_vals[i])**2 )) / (pitch[j]/200)
        D_vals.append(D_val)

    # Plot:
    plot_orbit(x_vals, y_vals, z_vals, param_label, pitch[j], dir)
    plot_dispersion(D_vals, z_vals, param_label, pitch[j], dir)
    plot_angular_momentum(Lz_vals, z_vals, param_label, pitch[j], dir)

    del x_vals, y_vals, z_vals, D_vals

##### ANIMATIONS #####

frame_duration = 200

# List of sim directories:
out_dirs = [main_dir+f'g4bl-output-sim{i+1}' for i in range(iterations)]

# List of orbit plot files:
orbit_plot_paths = []
for i in range(iterations):
    orbit_plot_paths.append(out_dirs[i]+'/orbit.png')

# Create animation of orbit over scan:
orbit_plots = [imageio.imread(img) for img in orbit_plot_paths]
# imageio.mimsave(fig_dir+'orbit_scan_dp_fine.gif', orbit_plots, duration=frame_duration, loop=0)
# imageio.mimsave(fig_dir+'orbit_scan_dp_coarse.gif', orbit_plots, duration=frame_duration, loop=0)
imageio.mimsave(fig_dir+'orbit_scan_sol_tilt.gif', orbit_plots, duration=frame_duration, loop=0)

# List of dispersion plot files:
disp_plot_paths = []
for i in range(iterations):
    disp_plot_paths.append(out_dirs[i]+'/dispersion.png')

# Create animation of dispersion over scan:
disp_plots = [imageio.imread(img) for img in disp_plot_paths]
# imageio.mimsave(fig_dir+'dispersion_scan_dp_fine.gif', disp_plots, duration=frame_duration, loop=0)
# imageio.mimsave(fig_dir+'dispersion_scan_dp_coarse.gif', disp_plots, duration=frame_duration, loop=0)
imageio.mimsave(fig_dir+'dispersion_scan_sol_tilt.gif', disp_plots, duration=frame_duration, loop=0)

# List of angular momentum plot files:
Lz_plot_paths = []
for i in range(iterations):
    Lz_plot_paths.append(out_dirs[i]+'/angular_momentum.png')

# Create animation of angular momentum over scan:
Lz_plots = [imageio.imread(img) for img in Lz_plot_paths]
imageio.mimsave(fig_dir+'Lz_scan_sol_tilt.gif', Lz_plots, duration=frame_duration, loop=0)
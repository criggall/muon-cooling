import numpy as np
from matplotlib import pyplot as plt
import os
import imageio.v2 as imageio

##### INPUTS #####

# # Choose to run flipped or not flipped lattice:
# polarity = 'flipped'
# # polarity = 'not_flipped'

# # Main directory:
# main_dir = '/Users/criggall/Documents/muon-cooling/Solenoid-Study/'

# # Define working directory:
# if polarity == 'flipped':
#      main_dir = main_dir+'flipped/'
# elif polarity == 'not_flipped':
#      main_dir = main_dir+'not-flipped/'

# main_dir = '/Users/criggall/Documents/muon-cooling/Solenoid-Study/single-coil/'
main_dir = '/Users/criggall/Documents/muon-cooling/Solenoid-Study/build-channel/'

# Parameter space scanned:
# nominal_tilt = -0.0025*180/np.pi
# tilt = np.linspace(0,nominal_tilt,10)
periods = np.arange(100,1500,100)

# Number of steps in scan:
# iterations = len(tilt)
iterations = len(periods)

# Number of solenoids in lattice:
n = 6

##### PLOTTING FUNCTIONS #####

units = '$^{\circ}$'

# Define function to plot orbit in xy-plane:
def plot_orbit(x_vals, y_vals, z_vals, param_val, dir):
    plt.clf()
    plt.scatter(x_vals,y_vals,s=1)
    # plt.title(f'solenoid tilt = {round(param_val,2)} {units}')
    plt.title(f'period = {param_val} mm')
    plt.xlabel('x (mm)')
    plt.ylabel('y (mm)')
    plt.savefig(dir+'orbit.png',dpi=300)
    plt.close()

# Define function to plot Lz along z:
def plot_angular_momentum(Lz_vals, z_vals, param_val, dir):
    plt.clf()
    plt.figure(figsize=(10,4))
    for i in range(n):
        plt.axvline(x=500+i*param_val,color='black',alpha=0.1)
    plt.scatter(z_vals, Lz_vals,s=0.5,color='black')
    # plt.title(f'solenoid tilt = {round(param_val,2)} {units}')
    plt.title(f'period = {param_val} mm')
    # plt.xlim(-1500,5500) # for 2 solenoids
    plt.xlim(-1500,8500)
    plt.ylim(-0.5,0.5)
    plt.xlabel('z (mm)')
    plt.ylabel('$L_z$ (mm*MeV/c)')
    plt.savefig(dir+'angular_momentum.png',dpi=300)
    plt.close()

# Define function to plot Bz along z:
def plot_B_field(Bz_vals, z_vals, param_val, dir):
    plt.clf()
    plt.figure(figsize=(10,4))
    for i in range(n):
        plt.axvline(x=500+i*param_val,color='black',alpha=0.1)
    plt.scatter(z_vals, Bz_vals,s=0.5,color='red')
    plt.title(f'period = {param_val} mm')
    # plt.xlim(-1500,5500) # for 2 solenoids
    plt.xlim(-1500,8500)
    plt.ylim(-1.5,1.5)
    plt.xlabel('z (mm)')
    plt.ylabel('$B_z$ (T)')
    plt.savefig(dir+'B_field.png',dpi=300)
    plt.close()

##### MAIN LOOP #####

min_Lz_vals = []; max_Lz_vals = []
for j in range(iterations):

    # Import data:
    # dir = f'{main_dir}sol_tilt_scan/g4bl-output-sim{j+1}/'
    dir = f'{main_dir}coil_spacing_scan/g4bl-output-sim{j+1}/'
    file = f'{dir}AllTracks.txt'
    data = np.loadtxt(file)

    # Values along channel:
    x_vals = []; y_vals = []; z_vals = []
    px_vals = []; py_vals = []
    Bz_vals = []
    Lz_vals = []
    for i in range(data.shape[0]):
        id = data[i][8]
        # if id == -2: # reference
        if id == 1: # with initial conditions
            x = data[i][0]; y = data[i][1]; z = data[i][2] # mm
            x_vals.append(x); y_vals.append(y); z_vals.append(z)
            px = data[i][3]; py = data[i][4] # MeV/c
            Bz_vals.append(data[i][14])
            Lz = x*py - y*px
            Lz_vals.append(Lz)
    
    # Find min and max Lz:
    min_Lz_vals.append(abs(np.min(Lz_vals)))
    max_Lz_vals.append(abs(np.max(Lz_vals)))

    # Plot:
    # plot_orbit(x_vals, y_vals, z_vals, tilt[j], dir)
    # plot_angular_momentum(Lz_vals, z_vals, tilt[j], dir)
    plot_angular_momentum(Lz_vals, z_vals, periods[j], dir)
    plot_B_field(Bz_vals, z_vals, periods[j], dir)

    del x_vals, y_vals, z_vals, px_vals, py_vals, Lz_vals

# # Plot magnitude of min and max Lz for each spacing:
# plt.figure()
# plt.plot(periods,min_Lz_vals,label='|min|',marker='.',color='blue')
# plt.plot(periods,max_Lz_vals,label='|max|',marker='.',color='red')
# plt.xlabel('period length (mm)')
# plt.ylabel('$L_z$ (mm*MeV/c)')
# plt.legend()
# plt.savefig(main_dir+'min_max_Lz.png',dpi=300)

##### ANIMATIONS #####

frame_duration = 400

# List of sim directories:
# out_dirs = [main_dir+f'sol_tilt_scan/g4bl-output-sim{i+1}' for i in range(iterations)]
out_dirs = [main_dir+f'coil_spacing_scan/g4bl-output-sim{i+1}' for i in range(iterations)]

# # List of orbit plot files:
# orbit_plot_paths = []
# for i in range(iterations):
#     orbit_plot_paths.append(out_dirs[i]+'/orbit.png')

# # Create animation of orbit over scan:
# orbit_plots = [imageio.imread(img) for img in orbit_plot_paths]
# imageio.mimsave(main_dir+'orbit_scan_sol_tilt.gif', orbit_plots, duration=frame_duration, loop=0)

# List of angular momentum plot files:
Lz_plot_paths = []
for i in range(iterations):
    Lz_plot_paths.append(out_dirs[i]+'/angular_momentum.png')

# Create animation of angular momentum over scan:
Lz_plots = [imageio.imread(img) for img in Lz_plot_paths]
# imageio.mimsave(main_dir+'Lz_scan_sol_tilt.gif', Lz_plots, duration=frame_duration, loop=0)
imageio.mimsave(main_dir+'Lz_scan_coil_spacing.gif', Lz_plots, duration=frame_duration, loop=0)

# List of B field plot files:
Bz_plot_paths = []
for i in range(iterations):
    Bz_plot_paths.append(out_dirs[i]+'/B_field.png')

# Create animation of B field over scan:
Bz_plots = [imageio.imread(img) for img in Bz_plot_paths]
imageio.mimsave(main_dir+'Bz_scan_coil_spacing.gif', Bz_plots, duration=frame_duration, loop=0)
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
# periods = np.arange(100,1500,100)
# tilts = np.arange(0,5,0.1)
# tilts = np.arange(0,0.3,0.01)
tilts = np.arange(0,60,5)

# Number of steps in scan:
# iterations = len(periods)
iterations = len(tilts)

# Number of solenoids in lattice:
# n = 1
n = 20

# Period length (if static):
period = 400

##### FUNCTION DEFINITIONS #####

# Define function to plot lines indicating solenoid centers:
def plot_solenoids(n):
    for i in range(n):
        plt.axvline(x=500+i*period,color='black',alpha=0.1)

# Define function to compute Br:
def calculate_Br(Bx_vals, By_vals, Bz_vals):
    Br_vals = []
    for i in range(len(z_vals)-1):
        x = x_vals[i]; y = y_vals[i]; z = z_vals[i]
        Bx = Bx_vals[i]; By = By_vals[i]; Bz = Bz_vals[i]
        r = np.sqrt(x**2+y**2)
        deltaBz = Bz_vals[i+1] - Bz
        deltaz = z_vals[i+1] - z
        if deltaz != 0:
            dBz_dz = deltaBz / deltaz
            Br = -r/2*dBz_dz
            Br_vals.append(Br)
        else:
            Br_vals.append(np.nan)
    Br_vals.append(np.nan)
    return Br_vals

units = '$^{\circ}$'

# Define function to plot orbit in xy-plane:
def plot_orbit(x_vals, y_vals, z_vals, param_val, dir):
    plt.clf()
    plt.scatter(x_vals,y_vals,s=1)
    # plt.title(f'period = {param_val} mm')
    plt.title(f'solenoid tilt = {round(param_val,3)} {units}')
    plt.xlabel('x (mm)')
    plt.ylabel('y (mm)')
    plt.savefig(dir+'orbit.png',dpi=300)
    plt.close()

# Define function to plot Lz along z:
def plot_angular_momentum(Lz_vals, z_vals, param_val, dir):
    plt.clf()
    plt.figure(figsize=(10,4))
    plot_solenoids(n)
    plt.scatter(z_vals, Lz_vals,s=0.5,color='black')
    # plt.title(f'period = {param_val} mm')
    plt.title(f'solenoid tilt = {round(param_val,3)} {units}')
    plt.xlim(-1500,8500)
    plt.xlabel('$z$ (mm)')
    plt.ylabel('$L_z$ (mm*MeV/c)')
    plt.savefig(dir+'angular_momentum.png',dpi=300)
    plt.close()

# Define function to plot Bx, By along z:
def plot_B_trans(Bx_vals,By_vals, z_vals, param_val, dir):
    plt.clf()
    plt.figure(figsize=(10,4))
    plot_solenoids(n)
    plt.scatter(z_vals, Bx_vals,s=0.5,color='red',label='$B_x$')
    plt.scatter(z_vals, By_vals,s=0.5,color='blue',label='$B_y$')
    plt.legend()
    plt.title(f'solenoid tilt = {round(param_val,3)} {units}')
    plt.xlim(-1500,8500)
    plt.xlabel('$z$ (mm)')
    plt.ylabel('$B$ (T)')
    plt.savefig(dir+'B_trans.png',dpi=300)
    plt.close()

# Define function to plot Bz along z:
def plot_Bz(Bz_vals, z_vals, param_val, dir):
    plt.clf()
    plt.figure(figsize=(10,4))
    plot_solenoids(n)
    plt.scatter(z_vals, Bz_vals,s=0.5,color='green',label='$B_z$')
    # plt.title(f'period = {param_val} mm')
    plt.title(f'solenoid tilt = {round(param_val,3)} {units}')
    plt.xlim(-1500,8500)
    plt.xlabel('$z$ (mm)')
    plt.ylabel('$B_z$ (T)')
    plt.savefig(dir+'B_z.png',dpi=300)
    plt.close()

# Define function to plot Br along z:
def plot_Br(Br_vals, z_vals, param_val, dir):
    plt.clf()
    plt.figure(figsize=(10,4))
    plot_solenoids(n)
    plt.scatter(z_vals, Br_vals, s=0.5,color='darkblue',label='$B_r$')
    plt.title(f'solenoid tilt = {round(param_val,3)} {units}')
    plt.xlim(-1500,8500)
    plt.xlabel('$z$ (mm)')
    plt.ylabel('$B_r$ (T)')
    plt.savefig(dir+'B_r.png',dpi=300)
    plt.close()

##### MAIN LOOP #####

min_Lz_vals = []; max_Lz_vals = []
for j in range(iterations):

    # Import data:
    # dir = f'{main_dir}sol_tilt_scan/g4bl-output-sim{j+1}/'
    # dir = f'{main_dir}coil_spacing_scan/g4bl-output-sim{j+1}/'
    dir = f'{main_dir}coil_tilt_scan/g4bl-output-sim{j+1}/'
    file = f'{dir}AllTracks.txt'
    data = np.loadtxt(file)

    # Values along channel:
    x_vals = []; y_vals = []; z_vals = []
    px_vals = []; py_vals = []
    Bz_vals = []; Bx_vals = []; By_vals = []
    Lz_vals = []
    for i in range(data.shape[0]):
        id = data[i][8]
        if id == -2: # reference
        # if id == 1: # with initial conditions
            x = data[i][0]; y = data[i][1]; z = data[i][2] # mm
            x_vals.append(x); y_vals.append(y); z_vals.append(z)
            px = data[i][3]; py = data[i][4] # MeV/c
            Bz_vals.append(data[i][14])
            Bx_vals.append(data[i][12])
            By_vals.append(data[i][13])
            Lz = x*py - y*px
            Lz_vals.append(Lz)
    
    # Find min and max Lz:
    min_Lz_vals.append(abs(np.min(Lz_vals)))
    max_Lz_vals.append(abs(np.max(Lz_vals)))

    # Compute Br:
    Br_vals = calculate_Br(Bx_vals, By_vals, Bz_vals)

    # Plot:
    # plot_orbit(x_vals, y_vals, z_vals, tilts[j], dir)
    # plot_angular_momentum(Lz_vals, z_vals, periods[j], dir)
    # plot_B_field(Bz_vals, z_vals, periods[j], dir)
    plot_angular_momentum(Lz_vals, z_vals, tilts[j], dir)
    plot_B_trans(Bx_vals, By_vals, z_vals, tilts[j], dir)
    plot_Bz(Bz_vals, z_vals, tilts[j], dir)
    plot_Br(Br_vals, z_vals, tilts[j], dir)

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
out_dirs = [main_dir+f'coil_tilt_scan/g4bl-output-sim{i+1}' for i in range(iterations)]
# out_dirs = [main_dir+f'coil_spacing_scan/g4bl-output-sim{i+1}' for i in range(iterations)]

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
imageio.mimsave(main_dir+'Lz_scan_sol_tilt.gif', Lz_plots, duration=frame_duration, loop=0)
# imageio.mimsave(main_dir+'Lz_scan_coil_spacing.gif', Lz_plots, duration=frame_duration, loop=0)

# List of Bx, By plot files:
Btrans_plot_paths = []
for i in range(iterations):
    Btrans_plot_paths.append(out_dirs[i]+'/B_trans.png')

# Create animation of Bx, By over scan:
Btrans_plots = [imageio.imread(img) for img in Btrans_plot_paths]
imageio.mimsave(main_dir+'Bx_By_scan_sol_tilt.gif', Btrans_plots, duration=frame_duration, loop=0)

# List of Bz plot files:
Bz_plot_paths = []
for i in range(iterations):
    Bz_plot_paths.append(out_dirs[i]+'/B_z.png')

# Create animation of Bz over scan:
Bz_plots = [imageio.imread(img) for img in Bz_plot_paths]
# imageio.mimsave(main_dir+'Bz_scan_coil_spacing.gif', Bz_plots, duration=frame_duration, loop=0)
imageio.mimsave(main_dir+'Bz_scan_sol_tilt.gif', Bz_plots, duration=frame_duration, loop=0)

# List of Br plot files:
Br_plot_paths = []
for i in range(iterations):
    Br_plot_paths.append(out_dirs[i]+'/B_r.png')

# Create animation of Br over scan:
Br_plots = [imageio.imread(img) for img in Br_plot_paths]
imageio.mimsave(main_dir+'Br_scan_sol_tilt.gif', Br_plots, duration=frame_duration, loop=0)
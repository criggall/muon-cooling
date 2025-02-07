import numpy as np
from matplotlib import pyplot as plt
import os
import imageio.v2 as imageio

# Input parameter space:
bls = np.arange(10,31)

# Name of parameter scan is over:
param_label = 'BLS'

# Specify whether reference particle, beam, or both:
ref_particle = True
beam = False

# Number of steps in scan:
iterations = len(bls)

# RF period (325 MHz frequency):
T = 1/(325*10**6)*10**9 # ns

# Length of period:
len_period = 4.2 # m

##### REFERENCE PARTICLE PLOTTING FUNCTIONS #####

# Define function to plot xy trajectory:
def plot_trajectory(x_vals, y_vals, z_vals, param_label, param_val, dir):
    plt.plot(z_vals,x_vals,color='red',label='x')
    plt.plot(z_vals,y_vals,color='blue',label='y')
    plt.title(f'{param_label} = {param_val}')
    plt.xlabel('z (m)')
    plt.ylabel('x, y (cm)')
    plt.legend(loc='upper left')
    plt.savefig(dir+'xy_trajectory.png',dpi=300)
    plt.close()

# Define function to plot B field:
def plot_B_field(x_vals, y_vals, z_vals, Bx_vals, By_vals, Bz_vals, end_period, param_label, param_val, dir):
    # plt.plot(z_vals[0:end_period],Bx_vals[0:end_period],color='green',label='200*Bx') # for first period only
    # plt.plot(z_vals[0:end_period],By_vals[0:end_period],color='blue',label='200*By')
    # plt.plot(z_vals[0:end_period],Bz_vals[0:end_period],color='red',label='Bz')
    plt.plot(z_vals,Bx_vals,color='green',label='200*Bx')
    plt.plot(z_vals,By_vals,color='blue',label='200*By')
    plt.plot(z_vals,Bz_vals,color='red',label='Bz')
    plt.title(f'{param_label} = {param_val}')
    plt.xlabel('z (m)')
    plt.ylabel('B (T)')
    plt.legend(loc='upper left')
    plt.savefig(dir+'B_field.png',dpi=300)
    plt.close()

##### BEAM PLOTTING FUNCTIONS #####

# Plot total momentum distributions:
# def plot_ptotal_dist():


##### MAIN LOOP #####

full_channel_indices = []
for j in range(iterations):

    if ref_particle == True:

        # Import data:
        dir = f'/Users/criggall/Documents/muon-cooling/Automate-G4bl/g4bl-output-sim{j+1}/'
        file = f'{dir}ReferenceParticle.txt'
        data = np.loadtxt(file)

        # Values along channel:
        x_vals = []; y_vals = []; z_vals = []
        px_vals = []; py_vals = []; pz_vals = []; ptotal_vals = []
        t_vals = []; mod_t_vals = []
        Bx_vals = []; By_vals = []; Bz_vals = []
        count = 0
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

            # Find index for last value in first period:
            if z > len_period and count == 0:
                end_period = i
                count += 1
        
        # Check if particle makes it to end of channel:
        if max(z_vals) >= 130:
            full_channel_indices.append(j)

        # Plot:
        plt.clf()
        plot_trajectory(x_vals, y_vals, z_vals, param_label, bls[j], dir)
        plt.clf()
        plot_B_field(x_vals, y_vals, z_vals, Bx_vals, By_vals, Bz_vals, end_period, param_label, bls[j], dir)

    # if beam == True:
        

##### ANIMATIONS #####

# Main directory:
main_dir = '/Users/criggall/Documents/muon-cooling/Automate-G4bl/'

# List of sim directories:
out_dirs = [main_dir+f'g4bl-output-sim{i+1}' for i in range(iterations)]

# List of plot files:
xy_plot_paths = []
B_plot_paths = []
for i in range(iterations):
    xy_plot_paths.append(out_dirs[i]+'/xy_trajectory.png')
    B_plot_paths.append(out_dirs[i]+'/B_field.png')

# Create animations:
xy_images = [imageio.imread(img) for img in xy_plot_paths]
B_images = [imageio.imread(img) for img in B_plot_paths]
imageio.mimsave(main_dir+'xy_animation.gif', xy_images, duration=100, loop=0)
imageio.mimsave(main_dir+'B_animation.gif', B_images, duration=100, loop=0)

# For ONLY configs where particle made it through the full channel:

# List of plot files:
xy_plot_paths = []
B_plot_paths = []
for i in full_channel_indices:
    xy_plot_paths.append(out_dirs[i]+'/xy_trajectory.png')
    B_plot_paths.append(out_dirs[i]+'/B_field.png')

# Create animations:
xy_images = [imageio.imread(img) for img in xy_plot_paths]
B_images = [imageio.imread(img) for img in B_plot_paths]
imageio.mimsave(main_dir+'xy_animation_full_channel.gif', xy_images, duration=1000, loop=0)
imageio.mimsave(main_dir+'B_animation_full_channel.gif', B_images, duration=1000, loop=0)
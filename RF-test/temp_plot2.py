import numpy as np
from matplotlib import pyplot as plt
import os
import imageio.v2 as imageio

##### INPUTS #####

# Main directory:
main_dir = '/Users/criggall/Documents/muon-cooling/RF-test/'

# Name of parameter scan is over:
param_label = 'reference particle momentum (MeV/c)'

# Parameter space scanned:
ref_p = np.arange(220,225.1,0.1)

# Number of steps in scan:
iterations = len(ref_p)

# Specify whether reference particle, beam, or both:
ref_particle = True
# beam = False # <-- Still have to add this option

# Region to plot:
plot_option = 3 # <-- = period number (0 for full channel)

# RF period (325 MHz frequency):
T = 1/(325*10**6)*10**9 # ns

# Length of period:
len_period = 4.2 # m

##### IMPORT DIGITIZED DATA FROM PAPER #####
dir = '/Users/criggall/Documents/muon-cooling/'
xdata = np.genfromtxt(dir+'paper_x_vs_z.csv',delimiter=',')
xdata_x = []; xdata_z = []
for i in range(len(xdata)):
    xdata_z.append(xdata[i][0])
    xdata_x.append(xdata[i][1])
ydata = np.genfromtxt(dir+'paper_y_vs_z.csv',delimiter=',')
ydata_y = []; ydata_z = []
for i in range(len(ydata)):
    ydata_z.append(ydata[i][0])
    ydata_y.append(ydata[i][1])

##### REFERENCE PARTICLE PLOTTING FUNCTIONS #####

# Define function to plot xy trajectory:
def plot_trajectory(x_vals, y_vals, z_vals, param_label, param_val, dir):
    if plot_option == 0:
        plt.figure(figsize = (12,6))
        plt.plot(z_vals,x_vals,color='red',label='x')
        plt.plot(z_vals,y_vals,color='blue',label='y')
        plt.ylim(-4,4)
    else:
        start_index = period_start_indices[plot_option-1]
        end_index = period_start_indices[plot_option]
        plt.figure(figsize = (12,6))
        plt.plot(z_vals[start_index:end_index],x_vals[start_index:end_index],color='red',label='x')
        plt.plot(z_vals[start_index:end_index],y_vals[start_index:end_index],color='blue',label='y')
        plt.ylim(-2,2)
    plt.title(f'{param_label} = {round(param_val,1)}')
    plt.xlabel('z (m)')
    plt.ylabel('x, y (cm)')
    plt.legend(loc='upper left')
    plt.savefig(dir+'xy_trajectory.png',dpi=300)
    plt.close()

# Define function to plot B field:
def plot_B_field(x_vals, y_vals, z_vals, Bx_vals, By_vals, Bz_vals, param_label, param_val, dir):
    if plot_option == 0:
        plt.figure(figsize = (12,6))
        plt.plot(z_vals,Bx_vals,color='green',label='200*Bx')
        plt.plot(z_vals,By_vals,color='blue',label='200*By')
        plt.plot(z_vals,Bz_vals,color='red',label='Bz')
        plt.ylim(-40,40)
    else:
        start_index = period_start_indices[plot_option-1]
        end_index = period_start_indices[plot_option]
        plt.figure(figsize = (12,6))
        plt.plot(z_vals[start_index:end_index],Bx_vals[start_index:end_index],color='green',label='200*Bx')
        plt.plot(z_vals[start_index:end_index],By_vals[start_index:end_index],color='blue',label='200*By')
        plt.plot(z_vals[start_index:end_index],Bz_vals[start_index:end_index],color='red',label='Bz')
        plt.ylim(-25,25)
    plt.title(f'{param_label} = {round(param_val,1)}')
    plt.xlabel('z (m)')
    plt.ylabel('B (T)')
    plt.legend(loc='upper left')
    plt.savefig(dir+'B_field.png',dpi=300)
    plt.close()

def plot_residual(x_vals, y_vals, z_vals, xdata_x, xdata_z, ydata_y, ydata_z, param_label, param_val):
    if plot_option == 0:

        positions = z_vals[start_index:end_index]
        for i in range(len(positions)):
            positions[i] = positions[i]*100 - 2*len_period*100

        x_paper_interp = np.interp(positions,xdata_z,xdata_x)
        y_paper_interp = np.interp(positions,ydata_z,ydata_y)

        x_diff = []
        y_diff = []
        x_total_residual = 0
        y_total_residual = 0
        for i in range(len(positions)):
            x_diff.append(x_vals_periodic[i] - x_paper_interp[i])
            y_diff.append(y_vals_periodic[i] - y_paper_interp[i])
            x_total_residual += x_diff[i]
            y_total_residual += y_diff[i]

    else:
    
        start_index = period_start_indices[plot_option-1]
        end_index = period_start_indices[plot_option]
        x_vals_periodic = x_vals[start_index:end_index]
        y_vals_periodic = y_vals[start_index:end_index]
        positions = z_vals[start_index:end_index]
        for i in range(len(positions)):
            positions[i] = positions[i]*100 - 2*len_period*100

        x_paper_interp = np.interp(positions,xdata_z,xdata_x)
        y_paper_interp = np.interp(positions,ydata_z,ydata_y)

        x_diff = []
        y_diff = []
        x_total_residual = 0
        y_total_residual = 0
        for i in range(len(positions)):
            x_diff.append(x_vals_periodic[i] - x_paper_interp[i])
            y_diff.append(y_vals_periodic[i] - y_paper_interp[i])
            x_total_residual += x_diff[i]
            y_total_residual += y_diff[i]

    plt.plot(positions,x_diff,color='red',label=f'x (total residual = {round(x_total_residual,2)})')
    plt.plot(positions,y_diff,color='blue',label=f'y (total residual = {round(y_total_residual,2)})')
    plt.ylabel('Residual')
    plt.xlabel('z (cm)')
    plt.ylim(-1,1)
    plt.title(f'{param_label} = {round(param_val,1)}')
    plt.legend()
    plt.savefig(dir+'residual.png',dpi=300)
    plt.close()

    return x_total_residual, y_total_residual

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
        file = f'{dir}ReferenceParticle.txt'
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
        
        # Check if particle makes it to end of channel:
        if max(z_vals) >= 130:
            full_channel_indices.append(j)

        # Plot:
        plt.clf()
        plot_trajectory(x_vals, y_vals, z_vals, param_label, ref_p[j], dir)
        plt.clf()
        plot_B_field(x_vals, y_vals, z_vals, Bx_vals, By_vals, Bz_vals, param_label, ref_p[j], dir)
        plt.clf()
        x_total_residual, y_total_residual = plot_residual(x_vals, y_vals, z_vals, xdata_x, xdata_z, ydata_y, ydata_z, param_label, ref_p[j])
        x_total_residuals.append(x_total_residual)
        y_total_residuals.append(y_total_residual)

    # if beam == True:


##### ANIMATIONS #####

# List of sim directories:
out_dirs = [main_dir+f'g4bl-output-sim{i+1}' for i in range(iterations)]

# List of plot files:
xy_plot_paths = []
B_plot_paths = []
residual_plot_paths = []
for i in range(iterations):
    xy_plot_paths.append(out_dirs[i]+'/xy_trajectory.png')
    B_plot_paths.append(out_dirs[i]+'/B_field.png')
    residual_plot_paths.append(out_dirs[i]+'/residual.png')

# Create animations:
xy_images = [imageio.imread(img) for img in xy_plot_paths]
B_images = [imageio.imread(img) for img in B_plot_paths]
residual_images = [imageio.imread(img) for img in residual_plot_paths]
frame_duration = 50
if plot_option == 0:
    imageio.mimsave(main_dir+'xy_animation_full_channel.gif', xy_images, duration=frame_duration, loop=0)
    imageio.mimsave(main_dir+'B_animation_full_channel.gif', B_images, duration=frame_duration, loop=0)
    imageio.mimsave(main_dir+'residual_animation_full_channel.gif', residual_images, duration=frame_duration, loop=0)
else:
    imageio.mimsave(main_dir+f'xy_animation_period_{plot_option}.gif', xy_images, duration=frame_duration, loop=0)
    imageio.mimsave(main_dir+f'B_animation_period_{plot_option}.gif', B_images, duration=frame_duration, loop=0)
    imageio.mimsave(main_dir+f'residual_animation_period_{plot_option}.gif', residual_images, duration=frame_duration, loop=0)

# Plot total residuals:
plt.clf()
plt.plot(ref_p, x_total_residuals, color='red', marker='.', label='x')
plt.plot(ref_p, y_total_residuals, color='blue', marker='.', label='y')
plt.xlabel('reference particle p (MeV/c)')
plt.ylabel('residual')
plt.legend()
plt.savefig(main_dir+'residual_vs_p.png',dpi=300)
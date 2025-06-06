import numpy as np
from matplotlib import pyplot as plt
import os
import math
from scipy.optimize import curve_fit

# Main directory:
main_dir = '/Users/criggall/Documents/muon-cooling/Solenoid-Study/single-coil/'

# Range of coil lengths scanned:
coil_lengths = np.arange(5,400,5) # mm

# Define function to compute focusing length from simulation:
def calc_f(z_vals, x_vals):
    start_f_index = np.nan; end_f_index = np.nan
    count1 = 0; count2 = 0
    for i in range(len(z_vals)):
        if abs(x_vals[i]) < 0.001 and count1 == 0:
            end_f_index = i
            count1 += 1
        if z_vals[i] > 0.5 and count2 == 0:
            start_f_index = i
            count2 += 1
    if math.isnan(start_f_index) == False and math.isnan(end_f_index) == False:
        f = (z_vals[int(end_f_index)] - z_vals[int(start_f_index)])*1000 # --> mm
    else:
        f = np.nan
    return f # mm

# Define function to compute focusing length from analytical approximation:
def pred_f(B,L):
    p = 200 # MeV/c
    p = p*5.344e-22 # MeV/c --> kg*m/s
    e = 1.602e-19 # C
    L = L/1000 # mm --> m
    f = (4*p**2)/(e**2*B**2*L)
    f = f*1000 # m --> mm
    return f # mm

f_vals = []
f_pred_vals = []
B2L_vals = []
for j in range(len(coil_lengths)):

    # Import data:
    dir = f'{main_dir}coil_length_scan/g4bl-output-sim{j+1}/'
    file = f'{dir}AllTracks.txt'
    data = np.loadtxt(file)

    # Values along channel:
    x_vals = []; y_vals = []; z_vals = []
    Bz_vals = []
    for i in range(data.shape[0]):
        id = data[i][8]
        if id == 1:
            x = data[i][0]*0.1; y = data[i][1]*0.1 # mm -> cm
            x_vals.append(x); y_vals.append(y)
            z = data[i][2]*0.001 # mm --> m
            z_vals.append(z)
            Bz_vals.append(data[i][14])
            del x, y, z

    # Compute focusing length:
    peak_Bz = np.max(Bz_vals)
    B2L_vals.append(coil_lengths[j]*peak_Bz**2)
    f_vals.append(calc_f(z_vals, x_vals))
    f_pred_vals.append(pred_f(peak_Bz, coil_lengths[j]))

    del x_vals, y_vals, z_vals

# Remove NaNs:
f_vals_temp = []
B2L_vals_temp = []
f_pred_vals_temp = []
for i in range(len(f_vals)):
    if np.isnan(f_vals[i]) == False:
        f_vals_temp.append(f_vals[i])
        B2L_vals_temp.append(B2L_vals[i])
        f_pred_vals_temp.append(f_pred_vals[i])
f_vals = f_vals_temp; del f_vals_temp
B2L_vals = B2L_vals_temp; del B2L_vals_temp
f_pred_vals = f_pred_vals_temp; del f_pred_vals_temp

# Define function to fit to f vs. B^2*L:
def fit_f(B2L, a):
    return a/B2L

# Perform fit:
popt, pcov = curve_fit(fit_f, B2L_vals, f_vals)
f_fit_vals = fit_f(B2L_vals, *popt)
print(f'Fit result: f = {np.round(*popt,1)}/(B^2*L)')

# # Plot focusing length vs. B^2*L:
# start = 0
# # plt.plot(B2L_vals[start:],f_pred_vals[start:],color='green',label='Theory',zorder=1)
# plt.plot(B2L_vals[start:],f_fit_vals[start:],color='red',label='Fit',zorder=1)
# plt.scatter(B2L_vals[start:],f_vals[start:],s=4,color='black',label='Simulation',zorder=2)
# plt.ylabel('Focusing length (mm)')
# plt.xlabel('$B^2 L$ (T$^2$mm)')
# plt.legend()
# plt.savefig(main_dir+'f_vs_length.png',dpi=300)

# Plot linearly:
B2L_vals_inv = []
for i in range(len(B2L_vals)):
    B2L_vals_inv.append(1/B2L_vals[i])
plt.plot(B2L_vals_inv,f_fit_vals,color='red',label='Fit',zorder=1)
plt.scatter(B2L_vals_inv,f_vals,s=4,color='black',label='Simulation',zorder=2)
plt.ylabel('Focusing length (mm)')
plt.xlabel('$1/(B^2 L)$ (T$^2$mm)')
plt.legend()
plt.savefig(main_dir+'f_vs_length.png',dpi=300)
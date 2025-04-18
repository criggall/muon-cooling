import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt

# Length of period:
len_period = 4.2 # m

# Absolute path to present directory:
main_dir = '/Users/criggall/Documents/muon-cooling/Simplified-HFOFO/'

# Load data for offset particle:
file = main_dir+'offset_particle_output/AllTracks_dxp_dyp_extended.txt'
data = np.loadtxt(file)

# Load data for reference particle:
file_ref = main_dir+'offset_particle_output/AllTracks_nominal_extended.txt'
data_ref = np.loadtxt(file_ref)

# Values along channel:
x_vals = []; y_vals = []; z_vals = []; t_vals = []
x_vals_ref = []; y_vals_ref = []; z_vals_ref = []; t_vals_ref = []
next_z = 0.5; period_start_z_vals = []; period_start_indices = [] # start solenoid placement at 0.5 m
for i in range(data.shape[0]):
    x_vals.append(data[i][0]*0.1) # mm -> cm
    y_vals.append(data[i][1]*0.1)
    z = data[i][2]*0.001 # mm -> m
    z_vals.append(z)
    t_vals.append(data[i][6]) # ns
    x_vals_ref.append(data_ref[i][0]*0.1) # mm -> cm
    y_vals_ref.append(data_ref[i][1]*0.1)
    z_ref = data_ref[i][2]*0.001 # mm -> m
    z_vals_ref.append(z_ref)
    t_vals_ref.append(data_ref[i][6]) # ns

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
t_vals = t_vals[:start_last_period]
x_vals_ref = x_vals_ref[:start_last_period]
y_vals_ref = y_vals_ref[:start_last_period]
z_vals_ref = z_vals_ref[:start_last_period]
t_vals_ref = t_vals_ref[:start_last_period]

# Find displacement:
displacements = []
for i in range(len(z_vals)):
    displacements.append(np.sqrt( (x_vals_ref[i] - x_vals[i])**2 + (y_vals_ref[i] - y_vals[i])**2 ))

# Evenly space z values:
z_uniform = np.linspace(np.min(z_vals), np.max(z_vals), len(z_vals))
displacements_interp_func = interp1d(z_vals, displacements)
displacements_uniform = displacements_interp_func(z_uniform)

# Amplitude and wave number values from FFT:
''' See plot-compare-particles notebook for FFT computation and plots'''
a = [3.4e-2, 6.3e-3]
k = [3.5, 7.0]

# Define functional form of observed dp particle trajectory:
def f(z, phi0, phi1, y):
    return ( a[0]*np.sin(k[0]*z+phi0) +  a[1]*np.sin(k[1]*z+phi1) + y)

# Fit to find phases:
initial = [0, 0, 0]
popt, pcov = curve_fit(f, z_uniform, displacements_uniform, p0=initial)

# Plot:
plt.figure(figsize=(12,5))
plt.plot(z_uniform, displacements_uniform, label='data',color='black')
plt.plot(z_uniform, f(z_uniform, *popt), label='fit', linestyle='--',color='green')
plt.xlabel('z (m)')
plt.ylabel('displacement (cm)')
plt.legend()
plt.savefig(main_dir+'Figures/fit_dxp_dyp.png')

# Print fitted functional form:
print(f'f(z) = {a[0]}sin({k[0]}z+{round(popt[0],2)}) +  {a[1]}sin({k[1]}z+{round(popt[1],2)}) + {round(popt[2],2)}')
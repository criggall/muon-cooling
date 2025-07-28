import numpy as np
from scipy.optimize import leastsq
import matplotlib.pyplot as plt

import sys
import os
parent_dir = os.path.abspath(os.path.join(os.getcwd(), '..'))
sys.path.append(parent_dir)
from field_map_g4bl import readFieldMapData

##### INPUTS #####

# z value at which to evaluate:
z_val = 0.0

# Path to field map file:
file = 'fieldmap.txt'

##### DEFINE FUNCTIONS #####

def fit_multipoles(x, y, Bx, By, max_order=10):

    ''' Fit B = Bx + i By to multipole expansion: B = sum_n c_n (x + i y)^n
    Inputs:
        x, y --> 2D arrays of positions
        Bx, By --> 2D arrays of transverse magnetic field components
        max_order --> maximum multipole order (positive integer)
    Returns:
        coefficients --> dictionary {n: c_n} for -max_order <= n <= max_order '''
    
    z = x + 1j * y
    z_flat = z.flatten()

    B = Bx + 1j * By
    B_flat = B.flatten()

    mask = np.abs(z_flat) > 1e-10
    z_valid = z_flat[mask]
    B_valid = B_flat[mask]

    start_n = 0
    n_values = list(range(start_n, max_order + 1))

    A = np.array([z_valid**n for n in n_values]).T
    coeffs, _, _, _ = np.linalg.lstsq(A, B_valid, rcond=None)

    return {n: coeffs[i] for i, n in enumerate(n_values)}

##### APPLY TO SIMULATION DATA #####

data = readFieldMapData(file)

x = np.unique(data['x'].values)
y = np.unique(data['y'].values)
X, Y = np.meshgrid(x, y)

data_slice = data[data['z'] == z_val]
Bx = data_slice.pivot_table(index='y', columns='x', values='Bx')
Bx = Bx.values
By = data_slice.pivot_table(index='y', columns='x', values='By')
By = By.values

coeffs = fit_multipoles(X, Y, Bx, By, max_order=5)

print("Multipole coefficients:")
for n in sorted(coeffs.keys()):
    print(f"n={n:2d} : {coeffs[n]:.4e}")

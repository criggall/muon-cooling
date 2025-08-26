import numpy as np
from scipy.optimize import leastsq
import cmath

def fitMultipoles(x, y, Bx, By, max_order=5):

    ''' Fit B = Bx + i By to multipole expansion:
        B(z) = sum_{n=1}^{max_order} c_n * (x + i y)^(n-1)
    Inputs:
        x, y --> 2D arrays of positions
        Bx, By --> 2D arrays of transverse magnetic field components
        max_order --> maximum multipole order (integer >= 1)
    Returns:
        coefficients --> dict {n: c_n} for 1 <= n <= max_order '''
    
    z = x + 1j * y
    z_flat = z.ravel()

    B = Bx + 1j * By
    B_flat = B.ravel()

    mask = np.abs(z_flat) > 1e-10
    z_valid = z_flat[mask]
    B_valid = B_flat[mask]

    n_values = list(range(1, max_order+1))
    A = np.array([z_valid**(n-1) for n in n_values]).T

    coeffs, _, _, _ = np.linalg.lstsq(A, B_valid, rcond=None)

    return {n: coeffs[i] for i, n in enumerate(n_values)}


def computeField(x, y, coeffs):

    ''' Compute transverse B field using coefficients from multipole expansion, as per:
        B(z) = sum_{n=1}^{max_order} c_n * (x + i y)^(n-1)
    Inputs:
        x, y --> 2D arrays of position values
        coeffs --> list of coefficients in multipole expansion, ordered by n
    Returns:
        Bx, By --> 2D array of transverse field values '''
    
    B = np.zeros((len(x), len(y)))
    phi = np.zeros((len(x), len(y)))
    for i in range(len(x)):
        for j in range(len(y)):

            z = x[i] + 1j * y[j]
    
            B_val = 0
            for n in range(len(coeffs)):
                B_val += coeffs[n] * z**n
            
            B[i][j] = abs(B_val)
            phi[i][j] = cmath.phase(B_val)

    return B, phi
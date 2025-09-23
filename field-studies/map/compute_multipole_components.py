import numpy as np
from scipy.optimize import leastsq


R_ref = 420 # HFOFO solenoid inner radius (mm)


def fitMultipoles(x, y, Bx, By, max_order=6, R=R_ref):

    ''' Fit B = Bx + i By to multipole expansion:
        B(z) = sum_{n=0}^{max_order} c_n * (z/R)^(n-1) where z = x + i y
    Inputs:
        x, y --> 2D arrays of positions
        Bx, By --> 2D arrays of transverse magnetic field components
        max_order --> maximum multipole order (integer >= 1)
        R --> reference radius (default is HFOFO solenoid inner radius)
    Returns:
        coefficients --> dict {n: c_n} for 1 <= n <= max_order '''
    
    z = x + 1j * y
    z_flat = z.ravel()

    B = By + 1j * Bx
    B_flat = B.ravel()

    mask = np.abs(z_flat) > 1e-10
    z_valid = z_flat[mask]
    B_valid = B_flat[mask]

    n_values = list(range(1, max_order+1))
    A = np.array([(z_valid/R)**(n-1) for n in n_values]).T

    coeffs, _, _, _ = np.linalg.lstsq(A, B_valid, rcond=None)

    return {n: coeffs[i] for i, n in enumerate(n_values)}


def computeField(x, y, coeffs_mag, coeffs_phase, R=R_ref):

    ''' Compute transverse B field using coefficients from multipole expansion, as per:
        B(z) = sum_{n=1}^{max_order} c_n * (z/R)^(n-1) for n >= 1
    Inputs:
        x, y --> lists of position values
        coeffs_mag --> list of magnitudes of coefficients in multipole expansion, ordered by n
        coeffs_phase --> list of phases of coefficients in multipole expansion, ordered by n
        R --> reference radius (default is HFOFO solenoid inner radius)
    Returns:
        Bx, By --> 2D array of transverse field values, as defined by B = B_y + i B_x convention '''
    
    Bx = np.zeros((len(x), len(y)))
    By = np.zeros((len(x), len(y)))
    for i in range(len(x)):
        for j in range(len(y)):

            z = x[i] + 1j * y[j]

            if np.isnan(z) == False:
                B_val = 0
                for n in range(1, len(coeffs_mag)+1):
                    coeff = coeffs_mag[n-1] * np.exp(1j * coeffs_phase[n-1])
                    B_val += coeff * (z/R)**(n-1)
            else:
                B_val = np.nan

            Bx[i][j] = B_val.imag
            By[i][j] = B_val.real

    return Bx, By
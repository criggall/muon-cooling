import numpy as np

def approximateBr(Bz_vals, x_vals, y_vals, z_vals):
    
    ''' Returns array of B_r values obtained via approximation for B_phi=0 case '''

    Br_vals = []
    for i in range(len(x_vals)-1):
        x = x_vals[i]; y = y_vals[i]
        r = np.sqrt(x**2+y**2)
        deltaBz = Bz_vals[i+1] - Bz_vals[i]
        deltaz = z_vals[i+1] - z_vals[i]
        if deltaz != 0:
            dBz_dz = deltaBz / deltaz
            Br = -r/2*dBz_dz
            Br_vals.append(Br)
        else:
            Br_vals.append(np.nan)
    Br_vals.append(np.nan)
    return Br_vals

def calculateBrBphi(Bx_vals, By_vals, x_vals, y_vals):

    ''' Returns arrays of B_r and B_phi values obtained via rotation of input B_x and B_y values '''

    Br_vals = []; Bphi_vals = []
    for i in range(len(x_vals)-1):
        x = x_vals[i]; y = y_vals[i]
        r = np.sqrt(x**2+y**2)
        theta = np.cos(x/r)
        Bx = Bx_vals[i]; By = By_vals[i]
        Br = Bx*np.cos(theta) + By*np.sin(theta)
        Bphi = Bx*np.sin(theta) + By*np.cos(theta)
        Br_vals.append(Br)
        Bphi_vals.append(Bphi)
    Br_vals.append(np.nan)
    Bphi_vals.append(np.nan)
    return Br_vals, Bphi_vals
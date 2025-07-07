import pandas as pd
import numpy as np
from scipy.optimize import minimize_scalar
import numpy.linalg as la

import matplotlib.pyplot as plt

##### INPUTS #####

# Parameters:
f_rf = 325e6 # RF frequency (Hz)
p0 = 250  # reference momentum (MeV/c)
Bz0 = 2  # peak on-axis B field (T)

# # Path to input file:
# file_path = "det.txt"

# # Plotting options:
# p_histo = False
# scatter = False
# phi_histo = False

##### CONSTANTS #####

m = 105.65837  # MeV/c^2
c = 299.792458  # mm/ns (speed of light)

omega_rf = 2 * np.pi * 325e6 * 1e-9 # RF frequency (Hz --> 1/ns)
period_rf = 2 * np.pi / omega_rf # RF period (ns)

kappa = 3 / p0
gamma0 = np.sqrt(1 + (p0 / m)**2)
beta0 = p0 / m / gamma0
v0 = c * beta0  # mm/ns

# Define 6x6 symplectic matrix:
S6 = np.array([
    [0, 1, 0, 0, 0, 0],
    [-1, 0, 0, 0, 0, 0],
    [0, 0, 0, 1, 0, 0],
    [0, 0, -1, 0, 0, 0],
    [0, 0, 0, 0, 0, 1],
    [0, 0, 0, 0, -1, 0]
])

def calculateEmittance(file_path, p_histo=False, scatter=False, phi_histo=False):

    ##### LOAD DATA #####

    # Read data from file:
    df = pd.read_csv(file_path, sep=r'\s+', comment="#",
                    names=["x", "y", "z", "Px", "Py", "Pz", "t", "PDGid", "EventID", "TrackId", "ParentID", "Weight"])
    df[["x", "y", "z"]] *= 10 # cm --> mm

    ''' The standardized units for this script are mm for position, MeV/c for momentum, and ns for time '''

    # Save only mu plus:
    muons = df[df["PDGid"] == -13].copy()
    muons["ptotal"] = np.sqrt(muons["Px"]**2 + muons["Py"]**2 + muons["Pz"]**2)

    # Plot total momentum distribution:
    if p_histo:
        plt.hist(muons["ptotal"], bins=80)
        plt.xlabel("$p_{total}$ (MeV/c)")
        plt.show()

    # Apply momentum cut:
    low_cut, high_cut = 150, 400
    good_muons = muons[(muons["ptotal"] > low_cut) & (muons["ptotal"] < high_cut)].copy()

    # Scatter plot of p_total vs. p_z:
    if scatter:
        plt.scatter(good_muons["Pz"], good_muons["ptotal"], s=1)
        plt.xlabel("$p_z$ (MeV/c)")
        plt.ylabel("$p_{total}$ (MeV/c)")
        plt.show()

    ##### RF TIMING #####

    t_good = good_muons["t"].to_numpy()

    # Compute optimal RF phase:
    def cosine_sum(phi_rf):
        return -np.sum(np.cos(omega_rf * t_good - phi_rf))
    result = minimize_scalar(cosine_sum, bounds=(0, 2*np.pi), method='bounded')
    optimal_phi_rf = result.x
    phi0 = np.pi - optimal_phi_rf

    # Find RF phases for all muons:
    phases_rad = (omega_rf * t_good + phi0) % (2 * np.pi)
    phases_deg = np.degrees(phases_rad)

    # Plot distribution of RF phases:
    if phi_histo:
        plt.hist(phases_deg, bins=30)
        plt.xlabel("$\phi_{RF}$ (degrees)")
        plt.tight_layout()
        plt.show()

    ##### COMPUTE VECTORS AND COVARIANCE #####

    psvec = []
    for i, row in enumerate(good_muons.itertuples()):
        x = row.x
        y = row.y
        z_mm = 0
        px = row.Px
        py = row.Py
        pz = row.Pz
        t = phases_rad[i] / omega_rf
        psvec.append([x, y, z_mm, px, py, pz, t])

    u_vectors = []
    outer_sum = np.zeros((6, 6))
    u_mean = np.zeros(6)
    n_good = len(good_muons)

    for ps in psvec:
        x, y, _, px, py, pz, t = ps

        # Compute components of vector potential:
        Ax = -0.5 * Bz0 * y
        Ay =  0.5 * Bz0 * x
        Az = 0.0

        # Construct u vector:
        u = np.array([
            x,
            px / p0 + kappa * Ax,
            y,
            py / p0 + kappa * Ay,
            -v0 * t,
            (np.sqrt(1 + (px**2 + py**2 + pz**2) / m**2) / gamma0 - 1) / (beta0**2)
        ])
        u_vectors.append(u)
        u_mean += u
        outer_sum += np.outer(u, u)

    # Compute covariance matrix:
    u_mean /= n_good
    cov_matrix = outer_sum / n_good - np.outer(u_mean, u_mean)

    ##### RMS EMITTANCES #####

    rms_emittances = []
    for i in range(3):
        i1 = 2 * i # 0, 2, 4
        i2 = 2 * i + 1 # 1, 3, 5
        sigma11 = cov_matrix[i1, i1]
        sigma22 = cov_matrix[i2, i2]
        sigma12 = cov_matrix[i1, i2]
        emittance = np.sqrt(sigma11 * sigma22 - sigma12**2)
        rms_emittances.append(emittance)

    # print(f"Emittances from RMS sigma matrix: {np.round(rms_emittances[0],1)}, {np.round(rms_emittances[1],1)}, {np.round(rms_emittances[2],1)} mm")





    ##### EMITTANCES FROM GAUSSIAN FIT #####

    eta0 = 1.0
    dampf = 0.85

    apr = u_mean.copy()      # initial mean (uav)
    sgpr = cov_matrix.copy() # initial covariance (sgmav)
    n_good = len(u_vectors)

    max_iter = 400

    for itr in range(max_iter):
        sgin = np.linalg.inv(sgpr)
        trpr = np.trace(sgpr)
        
        s1 = np.zeros(6)
        s2 = np.zeros((6, 6))
        s3 = 0.0
        
        for u in u_vectors:
            dz = u - apr
            exf = np.exp(-0.5 * dz @ (sgin @ dz))
            s1 += u * exf
            s2 += np.outer(dz, dz) * exf
            s3 += exf
        
        avr = s1 / s3
        sgnw = s2 / (s3 - eta0 * n_good / 16)
        
        # Check convergence based on trace of covariance matrix
        if abs(np.trace(sgnw) / trpr - 1) < 1e-7:
            break
        
        apr = (1 - dampf) * apr + dampf * avr
        sgpr = (1 - dampf) * sgpr + dampf * sgnw

    emittances = []
    for im in range(3):
        i = 2*im
        a = sgpr[i, i]
        d = sgpr[i+1, i+1]
        b = sgpr[i, i+1]
        emit = np.sqrt(a*d - b**2)
        emittances.append(emit)

    # print(f"Emittances from fitted sigma matrix: {np.round(emittances[0],1)}, {np.round(emittances[1],1)}, {np.round(emittances[2],1)} mm")

    ##### NORMALIZED EMITTANCES FROM RMS SIGMA MATRIX #####

    # Calculate eigenvalues and eigenvectors
    vals, vecs = la.eig(-cov_matrix @ S6)

    # Extract indices of eigenvalues with positive imaginary parts
    mdns = [i for i, val in enumerate(vals) if np.imag(val) > 0]

    # Compute normalized emittances (β0 * γ0) * Im(eigenvalues)
    normalized_rms_emittances = [beta0 * gamma0 * np.imag(vals[i]) for i in mdns]

    # print(f"Normalized emittances from RMS sigma matrix: {np.round(normalized_rms_emittances[0],1)}, {np.round(normalized_rms_emittances[1],1)}, {np.round(normalized_rms_emittances[2],1)} mm")

    ##### NORMALIZED EMITTANCES FROM FITTED SIGMA MATRIX #####

    # Compute eigenvalues and eigenvectors of -sgnw @ S6
    vals, vecs = np.linalg.eig(-sgnw @ S6)

    # Find indices where imaginary part of eigenvalue is positive
    mdns = [i for i, val in enumerate(vals) if np.imag(val) > 0]

    # Calculate normalized emittances from imaginary parts of eigenvalues
    normalized_emittances = beta0 * gamma0 * np.array([np.imag(vals[i]) for i in mdns])

    # print(f"Normalized emittances from fitted sigma matrix: {np.round(normalized_emittances[0],1)}, {np.round(normalized_emittances[1],1)}, {np.round(normalized_emittances[2],1)} mm")

    ##### BEAM SIZE #####

    # Beam sizes from sgmav (RMS covariance matrix)
    beam_sizes_sgmav = np.sqrt([cov_matrix[2*im - 2, 2*im - 2] for im in range(1, 4)])

    # Beam sizes from sgnw (Gaussian fit covariance matrix)
    beam_sizes_sgnw = np.sqrt([sgnw[2*im - 2, 2*im - 2] for im in range(1, 4)])

    # print(f"Beam sizes from RMS sigma matrix: {np.round(beam_sizes_sgmav[0],1)}, {np.round(beam_sizes_sgmav[1],1)}, {np.round(beam_sizes_sgmav[2],1)} mm")
    # print(f"Beam sizes from fitted sigma matrix: {np.round(beam_sizes_sgnw[0],1)}, {np.round(beam_sizes_sgnw[1],1)}, {np.round(beam_sizes_sgnw[2],1)} mm")

    ##### MISCELLANEOUS OTHER PARAMETERS #####

    # Effective betas:
    effective_betas = np.sqrt([sgnw[2*im - 2, 2*im - 2] / sgnw[2*im - 1, 2*im - 1] for im in range(1, 4)])
    # print(f"Effective betas: {np.round(effective_betas[0],1)}, {np.round(effective_betas[1],1)}, {np.round(effective_betas[2],1)} mm")

    # Beam central position:
    # print("Beam central position:", avr)

    # Mean momentum:
    mean_momentum = p0 * np.sqrt(1 + 2*avr[5] + (beta0**2) * (avr[5]**2))
    # print(f"Mean momentum: {np.round(mean_momentum,1)} MeV/c")




    return rms_emittances, emittances, normalized_rms_emittances, normalized_emittances
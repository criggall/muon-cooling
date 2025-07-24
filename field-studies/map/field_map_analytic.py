from math import pi
import numpy as np
from scipy import special
import mpmath
import math

### INPUT SOLENOID ROTATIONS ###

# Rotation about x-axis:
pitch = 0.0 # degrees

# Rotation about z-axis:
roll = 0.0 # degrees

### DEFINE CONSTANTS ###

mu = 1.256637e-6 # N/A^2

# Length of solenoid:
L = 300e-3 # m

# Solenoid inner radius:
a = 360e-3 # m

# Current density (A/m^2):
J = 80.46e-6

# Current/length (A/m):
I = 2*pi*a*J

### DEFINE USEFUL VARIABLES ###

def compute_xi(r, z):
    xi = np.sqrt((a+r)**2 + z**2)
    return xi

def compute_k(r, z):
    k = np.sqrt(4*a*r /((a+r)**2 + z**2))
    return k

def compute_c(r):
    c = -4*a*r/(a+r)**2
    return c

### DEFINE SPECIAL FUNCTIONS ###

def K(k):
    # K = special.ellipk(k)
    K = mpmath.ellipk(k**2)
    return K

def E(k):
    # E = special.ellipe(k)
    E = mpmath.ellipe(k**2)
    return E

def Pi(k, c):
    # Pi = mpmath.ellippi(k, c**2)
    # return Pi
    try:
        return mpmath.ellippi(k, c**2)
    except ValueError:
        print(f"Skipping invalid input: k={k}, c={c}")
        return mpmath.mpf(0)

### DEFINE FIELD FUNCTIONS ###

def bz(r, z):
    xi = compute_xi(r, z); k = compute_k(r, z); c = compute_c(r)
    bz = mu*I/pi * z*a/(xi * (a+r)) * (K(k) + (a-r)/(2*a) * (Pi(k,c) - K(k)))
    return bz

def br(r, z):
    xi = compute_xi(r, z); k = compute_k(r, z)
    br = mu*I/pi * xi/(4*r) * (2*(K(k) - E(k)) - k**2*K(k))
    return br

def Bz(r, z):
    Bz = -bz(r, z-L) + bz(r, z+L)
    return Bz

def Br(r, z):
    Br = br(r, z-L) + br(r, z+L)
    return Br

### RE-DEFINE FIELD FUNCTIONS WITH ROTATION ###

def Bz_rot(r, z, theta, pitch=0.0):
    Bz = Bz(r, z)
    Br = Br(r, z)
    Bzp = Br*np.sin(theta)*np.sin(pitch) + Bz*np.cos(pitch)
    return Bzp

def Br_rot(r, z, theta, pitch=0.0):
    Bz = Bz(r, z)
    Br = Br(r, z)
    Brp = Br*(np.cos(theta))**2 - Br*(np.sin(theta))**2*np.cos(pitch) - Bz*np.sin(theta)*np.sin(pitch)
    return Brp
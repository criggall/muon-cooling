import numpy as np
import pandas as pd


def readTraceData(file):

    ''' Returns a DataFrame with positions, momenta, time, ID, and B-field components from G4beamline trace output. 
    
    Inputs:
    file = G4beamline track file (ascii) '''
    
    data = np.loadtxt(file)
    
    columns1 = ['x', 'y', 'z', 'px', 'py', 'pz', 't', 'PDGid', 'EventID']
    columns2 = ['Bx', 'By', 'Bz']
    df = pd.DataFrame(np.column_stack([data[:, 0:9], data[:, 12:15]]), columns=columns1+columns2)

    df['ptotal'] = np.sqrt(df['px']**2 + df['py']**2 + df['pz']**2)
    df['xp'] = df['px'] / df['pz']
    df['yp'] = df['py'] / df['pz']
    df['r'] = np.sqrt(df['x']**2 + df['y']**2)
    df['theta'] = np.arctan2(df['y'], df['x'])
    df['rp'] = (df['x']*df['xp'] + df['y']*df['yp']) / df['r']
    df['thetap'] = df['x']*df['yp'] - df['y']*df['xp']
    df['Lz'] = df['x']*df['py'] - df['y']*df['px']

    return df


def readDetData(file, cuts=False, low_p_cut=0, high_p_cut=400):

    ''' Returns a DataFrame with positions, momenta, time, and ID from G4beamline detector output.

    Inputs:
    file = G4beamline detector output file (ascii) 
    cuts = Boolean, whether to apply cuts on total momentum
    low_p_cut = lower limit on total momentum (MeV/c)
    high_p_cut = upper limit on total momentum (MeV/c) '''
    
    data = np.loadtxt(file)

    columns = ['x','y','z','px','py','pz','t','PDGid','EventID']
    df = pd.DataFrame(data[:, :9], columns=columns)
    
    df['ptotal'] = np.sqrt(df['px']**2 + df['py']**2 + df['pz']**2)
    df['xp'] = df['px'] / df['pz']
    df['yp'] = df['py'] / df['pz']

    if cuts:
        df = df[(df['ptotal'] > low_p_cut) & (df['ptotal'] < high_p_cut)]

    return df
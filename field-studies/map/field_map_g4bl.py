import numpy as np
import pandas as pd


def readFieldMapData(file):

    ''' Returns a DataFrame with positions, time, and magnetic field components. From the fieldntuple G4bl command.

    The default units are mm for position and Tesla for B field.
    
    Inputs:
    file = G4beamline output file from fieldntuple (ascii) '''
    
    data = np.loadtxt(file)
    
    columns = ['x', 'y', 'z', 't', 'Bx', 'By', 'Bz']
    df = pd.DataFrame(data[:, 0:7], columns=columns)

    df['Br'] = np.sqrt(df['Bx']**2 + df['By']**2)

    return df
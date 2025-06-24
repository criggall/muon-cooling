import pandas as pd


def readECALC9(file):

    ''' Returns dictionary with longitudinal position and B field, emittances, and Twiss parameters from ECALC9 output 
        Note that the default unit for position is meters '''

    keys = ['z', 'Bz', 'eperp', 'elong', 'e6D', 'beta', 'alpha']
    values = { k : [] for k in keys }

    column_names = ['regn #', 'Z', 'Bz', 'eperp', 'elong', 'e6D', 'Ldim', 'Pzavg', 'beta', 'alpha', 'betaL', 'alphaL', 'n0', 'n1', 'n2', 'Lcan(eVs)', 'sigmaE', 'sigmaT', 'corrE', 'corrT', 'sigmaE_c', 'avg', 'yavg', 'Dx', 'Dy', 'Dr', 'Dr2']
    df = pd.read_csv(file, delim_whitespace=True, skiprows=12, names=column_names)

    values['z'].append(df['Z'].values)
    values['Bz'].append(df['Bz'].values)
    values['eperp'].append(df['eperp'].values)
    values['elong'].append(df['elong'].values)
    values['e6D'].append(df['e6D'].values)
    values['beta'].append(df['beta'].values)
    values['alpha'].append(df['alpha'].values)

    return values
import pandas as pd


def readECALC9(file):

    ''' Returns dictionary with longitudinal position and B field, emittances, and Twiss parameters from ECALC9 output.

        Note that the default unit for position is meters.
        
        Inputs:
        file = ecalc9.out '''

    keys = ['z', 'Bz', 'eperp', 'elong', 'e6D', 'beta', 'betaL']
    values = { k : [] for k in keys }

    column_names = ['regn #', 'Z', 'Bz', 'eperp', 'elong', 'e6D', 'Ldim', 'Pzavg', 'beta', 'alpha', 'betaL', 'alphaL', 'n0', 'n1', 'n2', 'Lcan(eVs)', 'sigmaE', 'sigmaT', 'corrE', 'corrT', 'sigmaE_c', 'avg', 'yavg', 'Dx', 'Dy', 'Dr', 'Dr2']
    df = pd.read_csv(file, sep='\s+', skiprows=12, names=column_names)

    values['z'].append(df['Z'].values)
    values['Bz'].append(df['Bz'].values)
    values['eperp'].append(df['eperp'].values)
    values['elong'].append(df['elong'].values)
    values['e6D'].append(df['e6D'].values)
    values['beta'].append(df['beta'].values)
    values['betaL'].append(df['betaL'].values)

    return values


def readECALC9F(file):

    ''' Returns dictionary with longitudinal position and B field, emittances, and Twiss parameters from ECALC9F output.

        Note that the default unit for position is meters.
        
        Inputs:
        file = ecalc9f.dat '''

    keys = ['z', 'Bz', 'eperp', 'elong', 'e6D', 'beta', 'betaL']
    values = { k : [] for k in keys }

    column_names = ['regn #', 'Z', 'Bz', 'eperp', 'elong', 'e6D', 'Ldim', 'Pzavg', 'beta', 'alpha', 'betaL', 'alphaL', 'n0', 'n1', 'n2', 'Lcan(eVs)', 'sigmaE', 'sigmaT', 'corrE', 'corrT', 'sigmaE_c', 'avg', 'yavg', 'Dx', 'Dy', 'Dr', 'Dr2']
    df = pd.read_csv(file, sep='\s+', skiprows=13, names=column_names)

    values['z'].append(df['Z'].values)
    values['Bz'].append(df['Bz'].values)
    values['eperp'].append(df['eperp'].values)
    values['elong'].append(df['elong'].values)
    values['e6D'].append(df['e6D'].values)
    values['beta'].append(df['beta'].values)
    values['betaL'].append(df['betaL'].values)

    return values


def readEMITCALC(file):

    ''' Returns dictionary with longitudinal position and B field, emittances, and Twiss parameters from EMITCALC output.

        Note that the default unit for position is meters.
        
        Inputs:
        file = emitcalc.out '''
    
    keys = ['z', 'Bz0', 'emit1', 'emit2', 'emit3', 'eperp', 'elong', 'beta', 'betaL', 'Dx', 'Dy']
    values = { k : [] for k in keys }

    column_names = ['regn #', 'Z', 'Bz0', '<Pz>', '<Lcanon>', '<x>', '<y>', '<t>', 'total #', '# in cuts', 'emit_1', 'emit_2', 'emit_3', 'e_perp', 'e_long', 'beta_perp', 'alpha_perp', 'Ldim', 'beta_long', 'alpha_long', 'D_x', 'D_y', 'Corr_E', 'Corr_t', 'sigma_E', 'asymmetry']
    df = pd.read_csv(file, sep='\s+', skiprows=15, names=column_names)

    values['z'].append(df['Z'].values)
    values['Bz0'].append(df['Bz0'].values)
    values['emit1'].append(df['emit_1'].values)
    values['emit2'].append(df['emit_2'].values)
    values['emit3'].append(df['emit_3'].values)
    values['eperp'].append(df['e_perp'].values)
    values['elong'].append(df['e_long'].values)
    values['beta'].append(df['beta_perp'].values)
    values['betaL'].append(df['beta_long'].values)
    values['Dx'].append(df['D_x'].values)
    values['Dy'].append(df['D_y'].values)

    return values
import numpy as np

def calculatePhaseAdvance(p, Bz0):

    ''' Inputs:
        p = reference momentum, in MeV/c (G4bl native units)
        Bz0 = peak longitudinal B field '''
    
    p = 200*10**-3 # --> GeV
    m = 105.7e-3 # GeV
    e = 0.303
    gamma = np.sqrt(1+(p/m)**2)
    beta = p/m/gamma
    phi = e*Bz0/(2*m*gamma*beta) # rad
    phi = phi*180/np.pi # --> deg
    return phi
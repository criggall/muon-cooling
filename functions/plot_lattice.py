def plot_solenoids(zoffset, n, d, L, ax):

    ''' Inputs:
        n = number of solenoids
        zoffset = z placement of first solenoid (mm)
        d = distance between solenoids (mm)
        L = length of solenoids (mm) '''

    for i in range(n):
        center = zoffset + i*d
        ax.axvspan(xmin=center-L/2, xmax=center+L/2, color='gray', alpha=0.2)
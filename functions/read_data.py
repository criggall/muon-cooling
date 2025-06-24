import numpy as np


def readTraceData(file):

    ''' Returns dictionary with positions, momenta, time, ID, and B field components from trace output '''

    data = np.loadtxt(file)
    keys = ['x','y','z','px','py','pz','ptotal','t','PDGid','EventID','Bx','By','Bz']
    values = { k : [] for k in keys }

    for i in range(data.shape[0]):
        x, y, z = data[i][0], data[i][1], data[i][2] # mm
        values['x'].append(x)
        values['y'].append(y)
        values['z'].append(z)
        px, py, pz = data[i][3], data[i][4], data[i][5] # MeV/c
        values['px'].append(px)
        values['py'].append(py) 
        values['pz'].append(pz)
        values['ptotal'].append(np.sqrt(px**2+py**2+pz**2))
        values['t'].append(data[i][6]) # ns
        Bx, By, Bz = data[i][12], data[i][13], data[i][14] # T
        values['Bx'].append(Bx)
        values['By'].append(By)
        values['Bz'].append(Bz)
        values['PDGid'].append(data[i][7])
        values['EventID'].append(data[i][8])
        del x, y, z, px, py, pz, Bx, By, Bz

    return values


def readDetData(file):

    ''' Returns dictionary with positions, momenta, time, and ID from detector output '''

    data = np.loadtxt(file)
    keys = ['x','y','z','px','py','pz','ptotal','t','PDGid','EventID']
    values = { k : [] for k in keys }

    for i in range(data.shape[0]):
        x, y, z = data[i][0], data[i][1], data[i][2] # mm
        values['x'].append(x)
        values['y'].append(y)
        values['z'].append(z)
        px, py, pz = data[i][3], data[i][4], data[i][5] # MeV/c
        values['px'].append(px)
        values['py'].append(py) 
        values['pz'].append(pz)
        values['ptotal'].append(np.sqrt(px**2+py**2+pz**2))
        values['t'].append(data[i][6]) # ns
        values['PDGid'].append(data[i][7])
        values['EventID'].append(data[i][8])
        del x, y, z, px, py, pz

    return values
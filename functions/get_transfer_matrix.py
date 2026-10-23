import numpy as np


def getInitVals(data_dvar, data_ref):
    out = np.zeros((6,1))
    out[0] = data_dvar['x'].values[0]-data_ref['x'].values[0]
    out[1] = data_dvar['xp'].values[0]-data_ref['xp'].values[0]
    out[2] = data_dvar['y'].values[0]-data_ref['y'].values[0]
    out[3] = data_dvar['yp'].values[0]-data_ref['yp'].values[0]
    out[4] = data_dvar['t'].values[0]-data_ref['t'].values[0]
    out[5] = (data_dvar['ptotal'].values[0]-data_ref['ptotal'].values[0])/data_ref['ptotal'].values[0]
    return out


def getFinVals(data_dvar, data_ref):
    out = np.zeros((6,1))
    out[0] = data_dvar['x'].values[420]-data_ref['x'].values[420]
    out[1] = data_dvar['xp'].values[420]-data_ref['xp'].values[420]
    out[2] = data_dvar['y'].values[420]-data_ref['y'].values[420]
    out[3] = data_dvar['yp'].values[420]-data_ref['yp'].values[420]
    out[4] = data_dvar['t'].values[420]-data_ref['t'].values[420]
    out[5] = (data_dvar['ptotal'].values[420]-data_ref['ptotal'].values[420])/data_ref['ptotal'].values[420]
    return out


def getTransferMatrix(data_ref, data_dx_pos, data_dx_neg, data_dxp_pos, data_dxp_neg, data_dy_pos, data_dy_neg, data_dyp_pos, data_dyp_neg, data_dt_pos, data_dt_neg, data_dp_pos, data_dp_neg):

    ''' Inputs: 
            - reference trajectory 
            - twelve tracks, w/ positive and negative offsets in each of the six phase space coordinates
    '''

    xstep = data_dx_pos['x'].values[0] - data_dx_neg['x'].values[0]
    xpstep = data_dxp_pos['xp'].values[0] - data_dxp_neg['xp'].values[0]
    ystep = data_dy_pos['y'].values[0] - data_dy_neg['y'].values[0]
    ypstep = data_dyp_pos['yp'].values[0] - data_dyp_neg['yp'].values[0]
    tstep = data_dt_pos['t'].values[0] - data_dt_neg['t'].values[0]
    deltastep = ((data_dp_pos['ptotal'].values[0]-data_ref['ptotal'].values[0]) - (data_dp_neg['ptotal'].values[0]-data_ref['ptotal'].values[0]))/data_ref['ptotal'].values[0]

    M = np.zeros((6,6))
    # column 1: x
    for i in range(6):
        M[i][0] = float((getFinVals(data_dx_pos,data_ref)[i] - getFinVals(data_dx_neg,data_ref)[i])/xstep)
    # column 2: xp
    for i in range(6):
        M[i][1] = float((getFinVals(data_dxp_pos,data_ref)[i] - getFinVals(data_dxp_neg,data_ref)[i])/xpstep)
    # column 3: y
    for i in range(6):
        M[i][2] = float((getFinVals(data_dy_pos,data_ref)[i] - getFinVals(data_dy_neg,data_ref)[i])/ystep)
    # column 4: yp
    for i in range(6):
        M[i][3] = float((getFinVals(data_dyp_pos,data_ref)[i] - getFinVals(data_dyp_neg,data_ref)[i])/ypstep)
    # column 5: t
    for i in range(6):
        M[i][4] = float((getFinVals(data_dt_pos,data_ref)[i] - getFinVals(data_dt_neg,data_ref)[i])/tstep)
    # column 6: delta
    for i in range(6):
        M[i][5] = float((getFinVals(data_dp_pos,data_ref)[i] - getFinVals(data_dp_neg,data_ref)[i])/deltastep)
    
    return M
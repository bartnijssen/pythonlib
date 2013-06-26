import datetime as dt
import numpy as np
import pandas as pd
import tempfile


def getsmlimits(soilfile, cells, celltemplate, nlayers=3):
    '''Determine minimum and maximum soil moisture storage'''
    # read lat, lon, soil depth, soil bulk density, soil mineral density,
    # and residual moisture for the deepest soil layer
    soil = np.loadtxt(soilfile, usecols=[1, 2, 3, (5*nlayers+10)-1, (8*nlayers+12)-1,
                                         (9*nlayers+12)-1, (12*nlayers+16)-1])
    maxsoilmoist = soil[:, 3] * (1 - soil[:, 4]/soil[:, 5]) * 1000.
    minsoilmoist = soil[:, 3] * soil[:, -1] * 1000.
    maxsoilmoist = maxsoilmoist.reshape(soil.shape[0], 1)
    minsoilmoist = minsoilmoist.reshape(soil.shape[0], 1)
    soil = np.hstack((soil[:, 0:3], minsoilmoist, maxsoilmoist))
    sm = {}
    for row in soil:
        cell = celltemplate.format(row[1], row[2])
        if cell in cells:
            sm[cell] = (row[0], row[3], row[4])
    return sm


def normalize_soilmoisture(df):
    '''Normalize soil moisture by substracting the minimum value'''
    for var in df:
        df[var] = df[var]-df[var].min()
    return df


def parsestatefile(statefile):
    '''Parse statefile for VIC 4.0.5 and return as a dictionary'''
    state = {}
    with open(statefile, 'r') as f:
        state['date'] = dt.datetime(*[int(i) for i in f.readline().split()])
        (state['nlayers'], state['nnodes']) = [int(i) for i in f.readline().split()]
        state['cells'] = {}
        state['order'] = []
        while 1:
            line = f.readline()
            if not line:
                break
            fields = line.split()
            (cellid, nveg, nbands) = [int(i) for i in fields[:3]]
            state['order'].append(cellid)
            state['cells'][cellid] = {}
            state['cells'][cellid]['nveg'] = nveg
            state['cells'][cellid]['nbands'] = nbands
            state['cells'][cellid]['nodedepths'] = [float(i) for i in fields[3:]]
            (state['cells'][cellid]['stillstorm'], state['cells'][cellid]['drytime']) = f.readline().split()
            state['cells'][cellid]['mu'] = []
            state['cells'][cellid]['vegband_data'] = []
            for i in range(nveg+1):
                state['cells'][cellid]['mu'].append(f.readline().strip())
                state['cells'][cellid]['vegband_data'].append([])
                for j in range(nbands):
                    state['cells'][cellid]['vegband_data'][i].append(f.readline().split())

    return state


def readfluxfile(filename):
    '''Read soil moisture from VIC flux file and return data frame'''
    colnames = ['year', 'month', 'day', 'c4', 'c5', 'c6', 'c7', 'c8', 'sm1', 'sm2', 'sm3', 'c12']
    df = pd.read_table(filename, usecols=[0, 1, 2, 8, 9, 10], parse_dates={'date': [0, 1, 2]},
                       header=None, names=colnames, index_col='date', delim_whitespace=True)
    return df


def readflowfile(filename, skiprows=0, conversion=1):
    '''Read routed or observed flows and return data frame'''
    # Pandas cannot handle the leading space in these files. Read the data,
    # strip the leading whitespace, write to temporary file and read using pandas
    colnames = ['year', 'month', 'day', 'flow']
    with open(filename, 'r') as f:
        read_data = f.readlines()
    read_data = [line.strip() for line in read_data]
    tmp = tempfile.TemporaryFile()
    try:
        for line in read_data:
            tmp.write("{}\n".format(line))
        tmp.seek(0)
        df = pd.read_table(tmp, parse_dates={'date': [0, 1, 2]}, header=None, names=colnames,
                           index_col='date', skiprows=skiprows, delim_whitespace=True)
    finally:
        tmp.close()
    df[df['flow'] < 0] = np.nan
    df['flow'] *= conversion
    return df


def writestatefile(statefile, state):
    '''Write statefile for VIC 4.0.5 based on a dictionary of state'''
    with open(statefile, 'w') as f:
        f.write(state['date'].strftime('%Y %m %d')+'\n')
        f.write('{} {}\n'.format(state['nlayers'], state['nnodes']))
        for cellid in state['order']:
            f.write('{} {} {} '.format(cellid, state['cells'][cellid]['nveg'], state['cells'][cellid]['nbands']))
            f.write(' '.join(['{}'.format(i) for i in state['cells'][cellid]['nodedepths']])+'\n')
            f.write('{} {}\n'.format(state['cells'][cellid]['stillstorm'], state['cells'][cellid]['drytime']))
            for i in range(state['cells'][cellid]['nveg']+1):
                f.write('{}\n'.format(state['cells'][cellid]['mu'][i]))
                for j in range(state['cells'][cellid]['nbands']):
                    f.write(' '.join(['{}'.format(k) for k in state['cells'][cellid]['vegband_data'][i][j]])+'\n')

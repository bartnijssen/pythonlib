#!/usr/bin/env python
"""Convert VIC ASCII output data into a netcdf file
"""

import argparse
import sys
from netCDF4 import Dataset
from netCDF4 import default_fillvals
import numpy as np
import glob
import os
import datetime
import time

def main():
    (inpath, ncfile, domain, precision, start, end, ncformat) = get_args()
    filelist = glob.glob(inpath + "*")
    ndays = end - start + datetime.timedelta(1)
    ndays = ndays.days
    # set up matrices for writing to netcdf file
    print 'ny: {} nx: {}'.format(domain['ny'], domain['nx'])
    data = {}
    data['evap'] = np.ones((ndays, domain['ny'], domain['nx']),
                           dtype='float32') * default_fillvals['f4']
    data['runoff'] = np.ones((ndays, domain['ny'], domain['nx']),
                             dtype='float32') * default_fillvals['f4']
    data['baseflow'] = np.ones((ndays, domain['ny'], domain['nx']),
                               dtype='float32') * default_fillvals['f4']
    data['soilmoist1'] = np.ones((ndays, domain['ny'], domain['nx']),
                                 dtype='float32') * default_fillvals['f4']
    data['soilmoist2'] = np.ones((ndays, domain['ny'], domain['nx']),
                                 dtype='float32') * default_fillvals['f4']
    data['soilmoist3'] = np.ones((ndays, domain['ny'], domain['nx']),
                                 dtype='float32') * default_fillvals['f4']
    data['swe'] = np.ones((ndays, domain['ny'], domain['nx']),
                          dtype='float32') * default_fillvals['f4']
    data['wdew'] = np.ones((ndays, domain['ny'], domain['nx']),
                           dtype='float32') * default_fillvals['f4']
    data['sensible'] = np.ones((ndays, domain['ny'], domain['nx']),
                               dtype='float32') * default_fillvals['f4']
    data['latent'] = np.ones((ndays, domain['ny'], domain['nx']),
                             dtype='float32') * default_fillvals['f4']
    data['grndflux'] = np.ones((ndays, domain['ny'], domain['nx']),
                               dtype='float32') * default_fillvals['f4']
    data['rnet'] = np.ones((ndays, domain['ny'], domain['nx']),
                           dtype='float32') * default_fillvals['f4']
    data['Trad'] = np.ones((ndays, domain['ny'], domain['nx']),
                           dtype='float32') * default_fillvals['f4']
    data['prcp'] = np.ones((ndays, domain['ny'], domain['nx']),
                           dtype='float32') * default_fillvals['f4']

    for file in filelist:
        print "Ingesting: {}".format(file)
        indata = np.loadtxt(file, dtype=None)
        fields = file.rsplit('_')
        lon = float(fields[-1])
        lat = float(fields[-2])
        x = int((lon-domain['west'])/domain['lonres'])
        y = int((lat-domain['south'])/domain['latres'])
        print 'lon: {} lat: {} x: {} y: {}'.format(lon, lat, x, y)
        data['evap'][:,y,x] = indata[:,3]
        data['runoff'][:,y,x] = indata[:,4]
        data['baseflow'][:,y,x] = indata[:,5]
        data['soilmoist1'][:,y,x] = indata[:,6]
        data['soilmoist2'][:,y,x] = indata[:,7]
        data['soilmoist3'][:,y,x] = indata[:,8]
        data['swe'][:,y,x] = indata[:,9]
        data['wdew'][:,y,x] = indata[:,10]
        data['sensible'][:,y,x] = indata[:,11]
        data['latent'][:,y,x] = indata[:,12]
        data['grndflux'][:,y,x] = indata[:,13]
        data['rnet'][:,y,x] = indata[:,14]
        data['Trad'][:,y,x] = indata[:,15]
        data['prcp'][:,y,x] = indata[:,16]

    write_netcdf(ncfile, ncformat, start, ndays, domain, data)

    return

def get_args():
    parser = argparse.ArgumentParser(description='Convert VIC ASCII met files '
                                     'to NetCDF')
    parser.add_argument('--netcdf', required=True, metavar='<netcdf file>',
                        help='netcdf file (output)')
    parser.add_argument('--input', required=True, metavar='<input path>',
                        help='input path with VIC met files, including leading '
                        'part of filename (before the latitude)')
    parser.add_argument('--domain', required=True, metavar='<domain>',
                        help='domain in following format: '
                        'west/east/resolution/south/north/resolution '
                        'all in degrees')
    parser.add_argument('--precision', required=True, metavar='<precision>',
                        help='number of digits after the decimal in the '
                        'VIC file names')
    parser.add_argument('--start', required=True, metavar='<date>',
                        help='start date in following format: YYYY/MM/DD')
    parser.add_argument('--end', required=True, metavar='<date>',
                        help='end date in following format: YYYY/MM/DD')
    parser.add_argument('--format', metavar='<NETCDF4|NETCDF4_CLASSIC|'
                        'NETCDF3_CLASSIC|NETCDF3_64BIT>',
                        help='NetCDF format. NETCDF4 is default')
    args = parser.parse_args()
    domain = parse_domain(args.domain)
    ncformat = 'NETCDF4'
    if (args.format):
        ncformat = args.format
    start = parse_date(args.start, '/')
    end = parse_date(args.end, '/')
    try:
        precision = int(args.precision)
        if (precision < 0):
            raise ValueError
    except:
        print "Error: precision must be an integer >= 0: {}".format(args.precision)

    return (args.input, args.netcdf, domain, precision, start, end, ncformat)

def parse_date(datestr, sep='-'):
    try:
        date = datetime.date(*[int(x) for x in datestr.split(sep)])
    except (ValueError, TypeError):
        sys.exit('Not a valid date: {}'.format(datestr))
    except:
        sys.exit('Unexpected error: {}'.format(sys.exc_info()[0]))
    return date

def parse_domain(str):
    domain = {}

    fields = str.split('/')
    if len(fields) != 6:
        sys.exit('Error: Domain {} should have 6 fields'.format(str))

    fields[:] = [float(field) for field in fields]
    domain['west'] = fields[0]
    domain['east'] = fields[1]
    domain['lonres'] = fields[2]
    domain['south'] = fields[3]
    domain['north'] = fields[4]
    domain['latres'] = fields[5]
    domain['nx'] = int((domain['east']-domain['west'])/domain['lonres'])
    domain['ny'] = int((domain['north']-domain['south'])/domain['latres'])

    return domain

def write_netcdf(ncfile, format, start, ndays, domain, data):
    print "Writing {}".format(ncfile)
    nc = Dataset(ncfile, 'w', format=format)
    nc.createDimension('time', None)
    nc.createDimension('lon', domain['nx'])
    nc.createDimension('lat', domain['ny'])

    # coordinate variables
    time = nc.createVariable('time', 'f4', 'time')
    lon = nc.createVariable('lon', 'f4', 'lon')
    lat = nc.createVariable('lat', 'f4', 'lat')

    time[:] = np.arange(0, ndays, 1)
    lon[:] = domain['west'] + np.arange(0, domain['nx']) * domain['lonres']
    lat[:] = domain['south'] + np.arange(0, domain['ny']) * domain['latres']

    for var in data.iterkeys():
        ncvar = nc.createVariable(var, 'f4', ('time','lat','lon',),
                                  fill_value=default_fillvals['f4'])
        ncvar[:] = data[var]

    attrtable = {
        'time': {
            'long_name': 'time',
            'units': 'days since {}'.format(start),
            'calendar': 'standard',
            },
        'lon': {
            'long_name': 'longitude',
            'units': 'degrees_east',
            },
        'lat': {
            'long_name': 'latitude',
            'units': 'degrees_north',
            },
        'evap': {
            'long_name': 'evapotransiration',
            'units': 'kg/m2/s'
            },
        'runoff': {
            'long_name': 'surface runoff',
            'units': 'kg/m2/s'
            },
        'baseflow': {
            'long_name': 'baseflow',
            'units': 'kg/m2/s'
            },
        'soilmoist1': {
            'long_name': 'soil moisture layer 1',
            'units': 'kg/m2'
            },
        'soilmoist2': {
            'long_name': 'soil moisture layer 2',
            'units': 'kg/m2'
            },
        'soilmoist3': {
            'long_name': 'soil moisture layer 3',
            'units': 'kg/m2'
            },
        'swe': {
            'long_name': 'snow water equivalent',
            'units': 'kg/m2'
            },
        'wdew': {
            'long_name': 'No idea',
            'units': '-'
            },
        'sensible': {
            'long_name': 'sensible heat flux',
            'units': 'W/m2'
            },
        'latent': {
            'long_name': 'latent heat flux',
            'units': 'W/m2'
            },
        'grndflux': {
            'long_name': 'ground heat flux',
            'units': 'W/m2'
            },
        'rnet': {
            'long_name': 'net radiation',
            'units': 'W/m2'
            },
        'Trad': {
            'long_name': 'radiative surface temperature',
            'units': 'K'
            },
        'prcp': {
            'long_name': 'precipitation',
            'units': 'mm/day',
            },
            }

    for var in attrtable.iterkeys():
        for attr in attrtable[var].iterkeys():
            nc.variables[var].setncattr(attr, attrtable[var][attr])

    nc.close()

if __name__ == "__main__":
    main()

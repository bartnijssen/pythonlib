#!/usr/bin/env python
"""Convert VIC meteorological ASCII data into a netcdf file
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
    data = {}
    data['prcp'] = np.ones((ndays, domain['ny'], domain['nx']),
                           dtype='float32') * default_fillvals['f4']
    data['tmax'] = np.ones((ndays, domain['ny'], domain['nx']),
                           dtype='float32') * default_fillvals['f4']
    data['tmin'] = np.ones((ndays, domain['ny'], domain['nx']),
                           dtype='float32') * default_fillvals['f4']
    data['wind'] = np.ones((ndays, domain['ny'], domain['nx']),
                           dtype='float32') * default_fillvals['f4']

    for file in filelist:
        print "Ingesting: {}".format(file)
        indata = np.loadtxt(file, dtype=None)
        fields = file.rsplit('_')
        lon = float(fields[-1])
        lat = float(fields[-2])
        x = int((lon-domain['west'])/domain['lonres'])
        y = int((lat-domain['south'])/domain['latres'])
        data['prcp'][:,y,x] = indata[:,0]
        data['tmax'][:,y,x] = indata[:,1]
        data['tmin'][:,y,x] = indata[:,2]
        data['wind'][:,y,x] = indata[:,3]

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
        'prcp': {
            'long_name': 'precipitation',
            'units': 'mm/day',
            'valid_min': 0.,
            'valid_max': 500.,
            },
        'tmin': {
            'long_name': 'daily minimum temperature',
            'units': 'degrees C',
            'valid_min': -60.,
            'valid_max': 40.,
            },
        'tmax': {
            'long_name': 'daily maximum temperature',
            'units': 'degrees C',
            'valid_min': -50.,
            'valid_max': 50.,
            },
        'wind': {
            'long_name': 'daily mean wind speed',
            'units': 'm/s',
            'valid_min': 0.,
            'valid_max': 50.,
            },
            }

    for var in attrtable.iterkeys():
        for attr in attrtable[var].iterkeys():
            nc.variables[var].setncattr(attr, attrtable[var][attr])

    nc.close()

if __name__ == "__main__":
    main()

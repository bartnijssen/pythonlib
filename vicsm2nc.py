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
import time as tm

def main():
    (inpath, ncfile, domain, start, end, model, ncformat) = get_args()
    filelist = glob.glob(inpath + "*")
    ndays = end - start + datetime.timedelta(1)
    ndays = ndays.days
    # set up matrices for writing to netcdf file
    data = {}
    data['sm'] = np.ones((ndays, domain['ny'], domain['nx']),
                           dtype='float32') * default_fillvals['f4']

    for file in filelist:
        fields = file.rsplit('_')
        lon = float(fields[-1])
        lat = float(fields[-2])
        x = int((lon-domain['west'])/domain['lonres'])
        y = int((lat-domain['south'])/domain['latres'])
        if (x >= 0) & (x < domain['nx']) & (y >= 0) & (y < domain['ny']):
            print "Ingesting: {}".format(file)
            indata = np.loadtxt(file, dtype=None)
            data['sm'][:,y,x] = indata[:,3]

    write_netcdf(ncfile, ncformat, start, ndays, domain, model, data)

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
    parser.add_argument('--start', required=True, metavar='<date>',
                        help='start date in following format: YYYY/MM/DD')
    parser.add_argument('--end', required=True, metavar='<date>',
                        help='end date in following format: YYYY/MM/DD')
    parser.add_argument('--model', required=True, metavar='<model>',
                        help='model - this is added as a global attribute')
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

    return (args.input, args.netcdf, domain, start, end, args.model, ncformat)

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
    domain['nx'] = int((domain['east']-domain['west'])/domain['lonres'])+1
    domain['ny'] = int((domain['north']-domain['south'])/domain['latres'])+1

    return domain

def write_netcdf(ncfile, format, start, ndays, domain, model, data):
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
        'sm': {
            'long_name': 'Soil moisture percentile',
            'units': '-',
            'valid_min': 0.,
            'valid_max': 1.,
            },
            }

    for var in attrtable.iterkeys():
        for attr in attrtable[var].iterkeys():
            nc.variables[var].setncattr(attr, attrtable[var][attr])

    nc.model = model
    nc.history = 'Created: {}\n'.format(tm.ctime(tm.time()))
    nc.history += ' '.join(sys.argv) + '\n'

    nc.close()

if __name__ == "__main__":
    main()

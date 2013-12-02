#!/usr/bin/env python
"""Convert VIC ASCII output data into a netcdf file that can be used by the
   python routing routines
"""

import argparse
import collections
import ConfigParser
import sys
from netCDF4 import Dataset
from netCDF4 import default_fillvals
import numpy as np
import glob
import datetime as dt
import time as tm
import re
from utils.parse import parse_date

verbose = False

def main():
    configfile = get_args()
    config = parse_config(configfile)
    if verbose:
        for section in config.iterkeys():
            print '{}'.format(section)
            for key, value in config[section].iteritems():
                print '{:18s}: {}'.format(key, value)
    inpath = config['OPTIONS']['inpath']
    if inpath[-1] != '/':
        inpath += '/'
    filelist = glob.glob(inpath + config['OPTIONS']['infiletemplate'] + "*")
    domainfile = config['OPTIONS']['domainfile']
    period = {}
    period['start'] = parse_date(config['OPTIONS']['start'])
    period['end'] = parse_date(config['OPTIONS']['end'])
    period['interval'] = dt.timedelta(int(config['OPTIONS']['timeinterval']))
    ncfile = config['OPTIONS']['ncfile']
    ncformat = config['OPTIONS']['ncformat']
    ncvarinfo = {}
    varlist = filterPick(config['OPTIONS'].keys(), re.compile('^ncvar.*_name$'))
    varlist = [x.split('_name')[0] for x in varlist]
    for var in varlist:
        ncvarinfo[var] = {}
        ncvarinfo[var]['name'] = config['OPTIONS'][var + '_name']
        ncvarinfo[var]['longname'] = config['OPTIONS'][var + '_longname']
        ncvarinfo[var]['column'] = config['OPTIONS'][var + '_column']
        ncvarinfo[var]['units'] = config['OPTIONS'][var + '_units']
        ncvarinfo[var]['divideby'] = float(config['OPTIONS'][var + '_divideby'])

    vicasc2nc(filelist, domainfile, period, ncfile,
              ncvarinfo, ncformat)

    return

def filterPick(mylist, regex, pos=0):
    matches = map(re.compile(regex).match, mylist)
    return [x.group(pos) for x in matches if x]

def get_args():
    parser = argparse.ArgumentParser(description='Convert VIC ASCII flux files '
                                     'to NetCDF for routing')
    parser.add_argument('--cfg', required=True, type=str,
                        metavar='<configuration file>',
                        help='configuration file')
    args = parser.parse_args()

    return args.cfg

def parse_config(configfile):
    configparser = ConfigParser.SafeConfigParser(allow_no_value=True)
    configparser.optionxform = str # preserve case of configuration keys
    if verbose:
        print 'Reading {}'.format(configfile)
    configparser.read(configfile)
    config = collections.OrderedDict() # preserve order of entries
    for section in configparser.sections():
        config[section] = collections.OrderedDict()
        for key, value in configparser.items(section):
            config[section][key] = value
    return config

def parse_domainfile(domainfile):
    domain = {}
    nc = Dataset(domainfile, 'r')
    lon = nc.variables['lon'][:]
    lat = nc.variables['lat'][:]
    domain['lonres'] = (lon.max()-lon.min())/(len(lon)-1)
    domain['west'] = lon.min() - 0.5 * domain['lonres']
    domain['east'] = lon.max() + 0.5 * domain['lonres']
    domain['latres'] = (lat.max()-lat.min())/(len(lat)-1)
    domain['south'] = lat.min() - 0.5 * domain['latres']
    domain['north'] = lat.max() + 0.5 * domain['latres']
    domain['nx'] = len(lon)
    domain['ny'] = len(lat)
    if verbose:
        for key, value in domain.iteritems():
            print '{:18s}: {}'.format(key, value)

    return domain

def vicasc2nc(filelist, domainfile, period, ncfile,
              ncvarinfo, ncformat):
    """Convert VIC ASCII to netCDF for use by the python routing program"""
    domain = parse_domainfile(domainfile)
    period['nsteps'] = ((period['end']-period['start']).days /
                        period['interval'].days + 1)
    varlist = ncvarinfo.keys()
    # set up matrices for writing to netcdf file
    data = {}
    for var in varlist:
        data[ncvarinfo[var]['name']] = (np.ones((period['nsteps'], domain['ny'],
                                        domain['nx']), dtype='float32') *
                                        default_fillvals['f4'])
    for file in filelist:
        if verbose:
            print "Ingesting: {}".format(file)
        indata = np.loadtxt(file, dtype=None)
        fields = file.rsplit('_')
        lon = float(fields[-1])
        lat = float(fields[-2])
        x = int((lon-domain['west'])/domain['lonres'])
        y = int((lat-domain['south'])/domain['latres'])
        if verbose:
            print 'lon: {} lat: {} x: {} y: {}'.format(lon, lat, x, y)
            print indata.shape
        for var in varlist:
            data[ncvarinfo[var]['name']][:,y,x] = (indata[:,ncvarinfo[var]['column']] /
                                                   ncvarinfo[var]['divideby'])

    # Loop over the days and write for each day.
    write_netcdf(ncfile, domain, period, ncvarinfo, ncformat, data)

    return

def write_netcdf(ncfile, domain, period, ncvarinfo, format, data):
    if verbose:
        print "Writing {}".format(ncfile)
    nc = Dataset(ncfile, 'w', format=format)
    nc.createDimension('time', None)
    nc.createDimension('lon', domain['nx'])
    nc.createDimension('lat', domain['ny'])

    # coordinate variables
    time = nc.createVariable('time', 'f4', 'time')
    lon = nc.createVariable('lon', 'f4', 'lon')
    lat = nc.createVariable('lat', 'f4', 'lat')

    time[:] = np.arange(0, period['nsteps'], period['interval'].days)
    lon[:] = domain['west'] + np.arange(0.5, domain['nx']) * domain['lonres']
    lat[:] = domain['south'] + np.arange(0.5, domain['ny']) * domain['latres']

    for var in data.iterkeys():
        if verbose:
            print 'var: {}'.format(var)
        ncvar = nc.createVariable(var, 'f4', ('time','lat','lon',),
                                  fill_value=default_fillvals['f4'])
        ncvar[:] = data[var]

    attrtable = {
        'time': {
            'long_name': 'time',
            'units': 'days since {}'.format(period['start']),
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
        }

    for var in ncvarinfo.keys():
        attrtable[ncvarinfo[var]['name']] = {
            'long_name': ncvarinfo[var]['longname'],
            'units': ncvarinfo[var]['units'],
            '_Fillvalue': default_fillvals['f4']
            }

    for var in attrtable.iterkeys():
        for attr in attrtable[var].iterkeys():
            nc.variables[var].setncattr(attr, attrtable[var][attr])

    nc.history = 'Created: {}\n'.format(tm.ctime(tm.time()))
    nc.history += ' '.join(sys.argv) + '\n'

    nc.close()

if __name__ == "__main__":
    main()

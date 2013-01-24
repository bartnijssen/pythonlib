#!/usr/bin/env python
"""Convert ASCII xyzz files from the surface water monitor to NetCDF
"""

import argparse
import sys
from netCDF4 import Dataset
from netCDF4 import default_fillvals
import numpy as np
import time

def main():
    parser = argparse.ArgumentParser(description='Convert ASCII xyzz files '
                                     'from the surface water monitor to '
                                     'NetCDF')
    parser.add_argument('--xyzz', required=True, metavar='<xyzz file>',
                        help='xyzz file from surface water monitor')
    parser.add_argument('--netcdf', required=True, metavar='<netcdf file>',
                        help='netcdf file (output)')
    parser.add_argument('--domain', required=True, metavar='<domain>',
                        help='domain in following format: '
                        'west/east/resolution/south/north/resolution '
                        'all in degrees')
    parser.add_argument('--date', required=True, metavar='<date>',
                        help='date in following format: YYYYMMDD')
    parser.add_argument('--format', metavar='<NETCDF4|NETCDF4_CLASSIC|'
                        'NETCDF3_CLASSIC|NETCDF3_64BIT>',
                        help='NetCDF format. NETCDF4 is default')


    args = parser.parse_args()

    domain = parse_domain(args.domain)
    domain['time'] = args.date
    format = 'NETCDF4'
    if (args.format):
        format = args.format
    xyzz = parse_xyzzfile(args.xyzz, domain)
    write_netcdf(args.netcdf, format, domain, xyzz)


def parse_domain(str):
    domain = {}

    fields = str.split('/')
    if len(fields) != 6:
        sys.exit('Error: Domain {} does not have enough fields'.format(str))

    fields[:] = [float(field) for field in fields]
    domain['west'] = fields[0]
    domain['east'] = fields[1]
    domain['lonres'] = fields[2]
    domain['south'] = fields[3]
    domain['north'] = fields[4]
    domain['latres'] = fields[5]
    domain['nx'] = (domain['east']-domain['west'])/domain['lonres']
    domain['ny'] = (domain['north']-domain['south'])/domain['latres']

    return domain

def parse_xyzzfile(infile, domain):
    xyzz = {}
    with open(infile, 'r') as f:
        contents = f.readlines()

    # Strip white space and only keep non-empty lines
    contents[:] = [line.strip() for line in contents]
    contents[:] = [line for line in contents if line]

    xyzz['nav_lon'] = np.ones((domain['ny'], domain['nx']),
                              dtype=np.float32) * default_fillvals['f4']
    xyzz['nav_lat'] = np.ones((domain['ny'], domain['nx']),
                              dtype='float32') * default_fillvals['f4']
    xyzz['soilsat_inst'] = np.ones((domain['ny'], domain['nx']),
                                 dtype='float32') * default_fillvals['f4']
    xyzz['soilsat_mean'] = np.ones((domain['ny'], domain['nx']),
                                 dtype='float32') * default_fillvals['f4']
    xyzz['soilsat_anom'] = np.ones((domain['ny'], domain['nx']),
                                 dtype='float32') * default_fillvals['f4']
    xyzz['soilsat_perc'] = np.ones((domain['ny'], domain['nx']),
                                 dtype='float32') * default_fillvals['f4']

    for line in contents:
        fields = [float(field) for field in line.split()]
        (lon, lat, soilsat_inst, soilsat_mean, soilsat_anom, zzz, soilsat_perc) = fields
        if (lon < domain['west'] or lon > domain['east'] or
            lat < domain['south'] or lat > domain['north']):
            sys.exit('Error: ({}E, {}N) outside domain {}/{}/{}/{}'.
                     format(lon, lat, domain['west'], domain['east'],
                            domain['south'], domain['north']))
        x = ((lon - 0.5 * domain['lonres']) - domain['west']) / domain['lonres']
        y = ((lat - 0.5 * domain['latres']) - domain['south']) / domain['latres']
        xyzz['nav_lon'][y, x] = lon
        xyzz['nav_lat'][y, x] = lat
        xyzz['soilsat_inst'][y, x] = soilsat_inst
        xyzz['soilsat_mean'][y, x] = soilsat_mean
        xyzz['soilsat_anom'][y, x] = soilsat_anom
        xyzz['soilsat_perc'][y, x] = soilsat_perc

    return xyzz

def write_netcdf(ncfile, format, domain, data):
    nc = Dataset(ncfile, 'w', format=format)
    nc.createDimension('x', domain['nx'])
    nc.createDimension('y', domain['ny'])

    nav_lon = nc.createVariable('nav_lon', 'f4', ('y','x',))
    nav_lat = nc.createVariable('nav_lat', 'f4', ('y','x',))
    soilsat_inst = nc.createVariable('soilsat_inst', 'f4', ('y','x',))
    soilsat_mean = nc.createVariable('soilsat_mean', 'f4', ('y','x',))
    soilsat_anom = nc.createVariable('soilsat_anom', 'f4', ('y','x',))
    soilsat_perc = nc.createVariable('soilsat_perc', 'f4', ('y','x',))

    nav_lon[:] = data['nav_lon']
    nav_lat[:] = data['nav_lat']
    soilsat_inst[:] = data['soilsat_inst']
    soilsat_mean[:] = data['soilsat_mean']
    soilsat_anom[:] = data['soilsat_anom']
    soilsat_perc[:] = data['soilsat_perc']

    nav_lon.units = 'degrees_east'
    nav_lat.units = 'degrees_north'
    soilsat_inst.units = '-'
    soilsat_mean.units = '-'
    soilsat_anom.units = '-'
    soilsat_perc.units = '-'

    nav_lon.long_name = 'Longitude'
    nav_lat.long_name = 'Latitude'
    soilsat_inst.long_name = 'Soil saturation'
    soilsat_mean.long_name = 'Soil saturation mean'
    soilsat_anom.long_name = 'Soil saturation anomaly'
    soilsat_perc.long_name = 'Soil saturation percentile'

    nav_lon.FillValue = default_fillvals['f4']
    nav_lat.FillValue = default_fillvals['f4']
    soilsat_inst.FillValue = default_fillvals['f4']
    soilsat_mean.FillValue = default_fillvals['f4']
    soilsat_anom.FillValue = default_fillvals['f4']
    soilsat_perc.FillValue = default_fillvals['f4']

    nav_lon.missing_value = 1.e+20
    nav_lat.missing_value = 1.e+20
    soilsat_inst.missing_value = 1.e+20
    soilsat_mean.missing_value = 1.e+20
    soilsat_anom.missing_value = 1.e+20
    soilsat_perc.missing_value = 1.e+20

    nav_lon.axis = 'YX'
    nav_lat.axis = 'YX'
    soilsat_inst.axis = 'YX'
    soilsat_mean.axis = 'YX'
    soilsat_anom.axis = 'YX'
    soilsat_perc.axis = 'YX'

    nav_lon.description = 'Longitude of grid cell center'
    nav_lat.description = 'Latitude of grid cell center'
    soilsat_inst.description = ('Simulated total column soil saturation for '
                              'a specific date')
    soilsat_mean.description = ('Long-term mean simulated total colunm soil '
                              'saturation for a specific day of the year. '
                              'The long-term mean is calculated as the 5 day '
                              'moving average centered on the current day. '
                              'The averaging pariod is 1916-2004')
    soilsat_anom.description = ('Total column soil saturation anomaly. Calculated '
                              'as soilsat_inst - soilsat_anom')
    soilsat_perc.description = ('Total column soil saturation percentile. This '
                              'value shows how often during the 1916-2004 '
                              'reference period the soil saturation on this '
                              'day of the year (using a 5 day centered '
                              'window) was less than soilsat_inst')

    nav_lon.valid_min = -180.
    nav_lat.valid_min = -90.
    soilsat_inst.valid_min = 0.
    soilsat_mean.valid_min = 0.
    soilsat_anom.valid_min = -1.
    soilsat_perc.valid_min = 0.

    nav_lon.valid_max = 180.
    nav_lat.valid_max = 90.
    soilsat_inst.valid_max = 1.
    soilsat_mean.valid_max = 1.
    soilsat_anom.valid_max = 1.
    soilsat_perc.valid_max = 1.

    nav_lon.modulo = 360.

    nc.history = 'Created ' + time.ctime(time.time())
    nc.description = ('Soil saturation data from the Variable Infiltration Model '
                      'as part of the operational surface water monitor. Note '
                      'that since this data is from the operational surface '
                      'water monitor, there may be occassional data problems. '
                      'Please check the data carefully')
    nc.source = ('Surface Water Monitor, Surface Water Hydrology Group, '
                 'University of Washington, Seattle, Washington, USA')
    nc.website = 'http://www.hydro.washington.edu/forecast/monitor'
    nc.contact = 'Bart Nijssen, email: nijssen@uw.edu'
    nc.history = ' '.join(sys.argv)
    nc.projection = 'Geographic'
    nc.resolution = ('Spatial resolution: Longitude ({} degrees), latitude ({} '
                     'degrees)'.format(domain['lonres'], domain['latres']))
    nc.date = 'YYYYMMDD: {}'.format(domain['time'])

    nc.close()

if __name__ == "__main__":
    main()

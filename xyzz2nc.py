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
    xyzz['soilm_inst'] = np.ones((domain['ny'], domain['nx']),
                                 dtype='float32') * default_fillvals['f4']
    xyzz['soilm_mean'] = np.ones((domain['ny'], domain['nx']),
                                 dtype='float32') * default_fillvals['f4']
    xyzz['soilm_anom'] = np.ones((domain['ny'], domain['nx']),
                                 dtype='float32') * default_fillvals['f4']
    xyzz['soilm_perc'] = np.ones((domain['ny'], domain['nx']),
                                 dtype='float32') * default_fillvals['f4']
    
    for line in contents:
        fields = [float(field) for field in line.split()]
        (lon, lat, soilm_inst, soilm_mean, soilm_anom, zzz, soilm_perc) = fields
        if (lon < domain['west'] or lon > domain['east'] or
            lat < domain['south'] or lat > domain['north']):
            sys.exit('Error: ({}E, {}N) outside domain {}/{}/{}/{}'.
                     format(lon, lat, domain['west'], domain['east'],
                            domain['south'], domain['north']))
        x = ((lon - 0.5 * domain['lonres']) - domain['west']) / domain['lonres']
        y = ((lat - 0.5 * domain['latres']) - domain['south']) / domain['latres']
        xyzz['nav_lon'][y, x] = lon
        xyzz['nav_lat'][y, x] = lat
        xyzz['soilm_inst'][y, x] = soilm_inst
        xyzz['soilm_mean'][y, x] = soilm_mean
        xyzz['soilm_anom'][y, x] = soilm_anom
        xyzz['soilm_perc'][y, x] = soilm_perc
        
    return xyzz

def write_netcdf(ncfile, format, domain, data):
    nc = Dataset(ncfile, 'w', format=format)
    nc.createDimension('x', domain['nx'])
    nc.createDimension('y', domain['ny'])

    nav_lon = nc.createVariable('nav_lon', 'f4', ('y','x',))
    nav_lat = nc.createVariable('nav_lat', 'f4', ('y','x',))
    soilm_inst = nc.createVariable('soilm_inst', 'f4', ('y','x',))
    soilm_mean = nc.createVariable('soilm_mean', 'f4', ('y','x',))
    soilm_anom = nc.createVariable('soilm_anom', 'f4', ('y','x',))
    soilm_perc = nc.createVariable('soilm_perc', 'f4', ('y','x',))

    nav_lon[:] = data['nav_lon']
    nav_lat[:] = data['nav_lat']
    soilm_inst[:] = data['soilm_inst']
    soilm_mean[:] = data['soilm_mean']
    soilm_anom[:] = data['soilm_anom']
    soilm_perc[:] = data['soilm_perc']

    nav_lon.units = 'degrees_east'
    nav_lat.units = 'degrees_north'
    soilm_inst.units = 'mm'
    soilm_mean.units = 'mm'
    soilm_anom.units = 'mm'
    soilm_perc.units = '-'

    nav_lon.long_name = 'Longitude'
    nav_lat.long_name = 'Latitude'
    soilm_inst.long_name = 'Soil moisture'
    soilm_mean.long_name = 'Soil moisture mean'
    soilm_anom.long_name = 'Soil moisture anomaly'
    soilm_perc.long_name = 'Soil moisture percentile'

    nav_lon.FillValue = default_fillvals['f4']
    nav_lat.FillValue = default_fillvals['f4']
    soilm_inst.FillValue = default_fillvals['f4']
    soilm_mean.FillValue = default_fillvals['f4']
    soilm_anom.FillValue = default_fillvals['f4']
    soilm_perc.FillValue = default_fillvals['f4']

    nav_lon.missing_value = 1.e+20
    nav_lat.missing_value = 1.e+20
    soilm_inst.missing_value = 1.e+20
    soilm_mean.missing_value = 1.e+20
    soilm_anom.missing_value = 1.e+20
    soilm_perc.missing_value = 1.e+20

    nav_lon.axis = 'YX'
    nav_lat.axis = 'YX'
    soilm_inst.axis = 'YX'
    soilm_mean.axis = 'YX'
    soilm_anom.axis = 'YX'
    soilm_perc.axis = 'YX'

    nav_lon.description = 'Longitude of grid cell center'
    nav_lat.description = 'Latitude of grid cell center'
    soilm_inst.description = ('Simulated total column soil moisture for '
                              'a specific date')
    soilm_mean.description = ('Long-term mean simulated total colunm soil '
                              'moisture for a specific day of the year. '
                              'The long-term mean is calculated as the 5 day '
                              'moving average centered on the current day. '
                              'The averaging pariod is 1916-2004')
    soilm_anom.description = ('Total column soil moisture anomaly. Calculated '
                              'as soilm_inst - soilm_anom')
    soilm_perc.description = ('Total column soil moisture percentile. This '
                              'value shows how often during the 1916-2004 '
                              'reference period the soil moisture on this '
                              'day of the year (using a 5 day centered '
                              'window) was less than soilm_inst')

    nav_lon.valid_min = -180.
    nav_lat.valid_min = -90.
    soilm_inst.valid_min = 0.
    soilm_mean.valid_min = 0.
    soilm_perc.valid_min = 0.

    nav_lon.valid_max = 180.
    nav_lat.valid_max = 90.
    soilm_perc.valid_max = 1.

    nav_lon.modulo = 360.

    nc.history = 'Created ' + time.ctime(time.time())
    nc.description = ('Soil moisture data from the Variable Infiltration Model '
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

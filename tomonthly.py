#!/usr/bin/env python
"""Convert submonthly netcdf files to monthly by averaging
"""

import argparse
import calendar
import datetime
from netCDF4 import Dataset
from netCDF4 import default_fillvals
import netcdftime
import numpy as np
import os
import sys
import subprocess
import tempfile

ncra = '/opt/local/bin/ncra'
ncrcat = '/opt/local/bin/ncrcat'
ncks = '/opt/local/bin/ncks'

def main():
    ncin, ncout = get_args()
    start, end = get_daterange(ncin)
    now = datetime.datetime(start.year, start.month, 1)
    tobedeleted = []
    while now < end:
        eom = now + datetime.timedelta(days =
                                       calendar.monthrange(now.year,
                                                           now.month)[1]-1)
        fname = get_tempfilename('.nc')
        fargs = [ncks, '-O', '-dtime,{:%Y-%m-%d},{:%Y-%m-%d}'.format(now, eom),
                 ncin, fname]
        returnval = subprocess.call(fargs)
        if returnval != 0:
            sys.exit('Error executing: {}'.format(' '.join(fargs)))

        tmpname = 'tmp-{:%Y-%m}.nc'.format(now)
        tobedeleted.append(tmpname)
        fargs = [ncra, '-O', fname, tmpname]
        returnval = subprocess.call(fargs)
        if returnval != 0:
            sys.exit('Error executing: {}'.format(' '.join(fargs)))

        remove_file(fname)
        now += datetime.timedelta(days=calendar.monthrange(now.year,
                                                           now.month)[1])

    fargs = list(flatten([ncrcat, '-O', tobedeleted, ncout]))
    returnval = subprocess.call(fargs)
    if returnval != 0:
        sys.exit('Error executing: {}'.format(' '.join(fargs)))
    for file in tobedeleted:
        try:
            os.remove(file)
        except OSError as e:
            print 'Error removing file: {}'.format(e)
            pass
            
    return

def flatten(*args):
    for x in args:
        if hasattr(x, '__iter__'):
            for y in flatten(*x):
                yield y
        else:
            yield x

def get_args():
    parser = argparse.ArgumentParser(description='Aggregate to monthly. '
                       'The netCDF file must have time units that are '
                       'understood by nco')
    parser.add_argument('--input', required=True, metavar='<input netcdf file>',
                        help='netcdf file (input)')
    parser.add_argument('--output', required=True, metavar='<output netcdf file>',
                        help='netcdf file (output)')
    args = parser.parse_args()
            
    return (args.input, args.output)


def get_daterange(ncin):
    nc = Dataset(ncin, 'r')

    # get the time variable (note that it must be named 'time'
    time = nc.variables['time']

    # initialize utime, which allows conversion of the time axis
    nctime = netcdftime.utime(time.units)

    start = nctime.num2date(time[0])
    end = nctime.num2date(time[-1])
    nc.close()
    return (start, end)

def get_tempfilename(suffix='.tmp'):
        fh, fname = tempfile.mkstemp(suffix=suffix)
        outsock = os.fdopen(fh,'w')
        outsock.close()
        return fname
        
def remove_file(file):
    try:
        os.remove(file)
    except OSError as e:
        print 'Error removing file: {}'.format(e)
        pass

    
if __name__ == "__main__":
    main()

#!/usr/bin/env python

import os
import subprocess
import tempfile
import sys

startyear = 1950
endyear = 1999

extract = '/home/nijssen/apps/bin/extractmonthly.py'
exepath = '/usr/bin'
ncks = exepath + '/ncks'
ncra = exepath + '/ncra'
ncrcat = exepath + '/ncrcat'
ncrename = exepath + '/ncrename'
ncatted = exepath + '/ncatted'

ncvars = ['Prcp', 'Tmin', 'Tmax', 'Wind']
mfiletemplate = 'gridded_obs.daily.{}.{:4d}.nc'
tmpfiletemplate = 'tmp.{}.{:4d}'
lfiletemplate = 'maurer.daily_data_{:4d}{:02d}.nc'
lmfiletemplate = 'maurer.monthly_data_{:4d}{:02d}.nc'
verbose = 1

def main():
    for year in range(startyear, endyear+1):
        if verbose:
            print 'Processing year: {}'.format(year)
        for var in ncvars:
            if verbose:
                print '\tExtracting months for {}'.format(var)
            fargs = [extract, '--input', mfiletemplate.format(var, year),
                     '--output', tmpfiletemplate.format(var, year)]
            command(fargs)
        for month in range(1,13):
            if verbose:
                print '\tProcessing month {:02d}'.format(month)
            if verbose:
                print '\t\tCombining variables'
            ncfile = lfiletemplate.format(year, month)
            for var in ncvars:
                mfile = tmpfiletemplate.format(var, year) + '_{:02d}.nc'.format(month)
                fargs = [ncks, '-A', mfile, ncfile]
                command(fargs)
                remove_file(mfile)
            # now modify the files
            if verbose:
                print '\t\tEditing dimensions, variables and attributes'
            tmpfname = get_tempfilename(suffix='.nc')
            fargs = [ncks, '-O', '-x', '-vbounds_latitude,bounds_longitude', ncfile, tmpfname]
            command(fargs)
            fargs = ['ncrename', '-dlongitude,lon', '-dlatitude,lat',
                     '-vlatitude,lat', '-vlongitude,lon', '-vPrcp,Prec', tmpfname]
            command(fargs)
            fargs = ['ncatted', '-O', '-a_FillValue,,o,f,1.e20', '-aaxis,,d,,',
                     '-aassociate,,d,,', '-abounds,,d,,', tmpfname, ncfile]
            command(fargs)
            # calculating monthly average
            if verbose:
                print '\t\tCalculating monthly average file'
            fargs = ['ncra', ncfile, lmfiletemplate.format(year, month)]
            command(fargs)
            remove_file(tmpfname)

def command(fargs):
    returnval = subprocess.call(fargs)
    if returnval != 0:
        sys.exit('Error executing: {}'.format(' '.join(fargs)))

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

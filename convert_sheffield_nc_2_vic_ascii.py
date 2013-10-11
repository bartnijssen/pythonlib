#!/usr/bin/env python

from netCDF4 import Dataset
import numpy as np
import os.path

startyear = 1948
endyear = 2008
basedir = '/raid3/nijssen/sheffield'
localdir = '/state/partition1/nijssen/sheffield'
ncdir = basedir
ascdir = localdir + '/' + 'ascii_masked'
infiletemplate = '{}/{}_masked/{}_masked_daily_{:04d}-{:04d}.nc'
outfiletemplate = '{}/data_{:0.4f}_{:0.4f}'

fvars = ['prcp', 'tmax', 'tmin', 'wind']
latvar = 'latitude'
lonvar = 'longitude'
timevar = 'time'

verbose = 0

def main():
    if (os.path.isdir(ascdir)):
        pass
    else:
        os.makedirs(ascdir)


    filename = infiletemplate.format(ncdir, fvars[0], fvars[0], startyear, startyear)
    if verbose:
        print filename
    nc = Dataset(filename, 'r')
    lons = nc.variables[lonvar][:]
    lats = nc.variables[latvar][:]
    nc.close()

    for year in range(startyear, endyear+1):
        print 'Processing year: {}'.format(year)
        count = -1
        filename = infiletemplate.format(ncdir, fvars[0], fvars[0], year, year)
        dim = read_nc_dimension(filename)
        forc_data = np.ma.empty((len(fvars), dim[0], dim[1], dim[2]))
        for varname in fvars:
            count+=1
            filename = infiletemplate.format(ncdir, varname, varname, year, year)
            if verbose:
                print "Reading %d %s" %(count, filename)
            dim = read_nc_dimension(filename)
            forc_data[count, :, :, :] = read_nc_data(filename, varname)
            if verbose:
                print "Finished Reading %s" %filename

        if verbose:
            print "Ascii data now"
        for i in range(len(lats)):
            for j in range(len(lons)):
                if not np.isfinite(forc_data[0, 0, i, j]):
                    if verbose:
                        print '\tSkipping {:.4f} {:.4f}'.format(lats[i], lons[j])
                    continue
                outfile = outfiletemplate.format(ascdir, lats[i], lons[j])
                if verbose:
                    print "Going to write: {}".format(outfile)
                if os.path.isfile(outfile):
                    ## Overwrite if this is the first year of forcings
                    if (year == startyear):
                        file_out = file(outfile, 'w')
                    else:
                        file_out = file(outfile, 'a')
                else:
                    file_out = file(outfile, 'w')

                #print "Writing %s" %outfile
                temp_data = np.empty((len(fvars), dim[0]))
                temp_data[0:len(fvars),:] = forc_data[0:len(fvars), :, i, j]
                np.savetxt(file_out, np.transpose(temp_data), fmt="%.6e", newline='\n')
                file_out.close()

def read_nc_dimension(nc_filename):
    nc = Dataset(nc_filename, 'r')
    lats = nc.variables[latvar][:]
    lons = nc.variables[lonvar][:]
    time = nc.variables[timevar][:]
    dimensions=[len(time), len(lats), len(lons)]
    nc.close()
    return dimensions

def read_nc_data(nc_filename, varname):
    nc = Dataset(nc_filename, 'r')
    data = nc.variables[varname][:]
    nc.close()
    return data


if __name__ == "__main__":
    main()
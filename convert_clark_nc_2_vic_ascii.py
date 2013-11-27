#!/usr/bin/env python

from netCDF4 import Dataset
import numpy as np
import pandas as pd
import sys

startyear = 1980
endyear = 2010
nsamples = 100
basedir = '/Users/nijssen/data/vicpf/data/gunnison/vic/randomfields'
ncdir = basedir + '/nc'
ascdir = basedir + '/asc/id_{:03d}'
infiletemplate = 'gunnison_prcp_temp_{}_{}.nc_{}'.format(startyear, endyear, '{:03d}')
outfiletemplate = 'data_{:.4f}_{:.4f}'
windfile = '/Users/nijssen/Dropbox/vicpf/data/gunnison/vic/retro/input/gunnison_retro_mean_wind.asc'

latvar = 'latitude'
lonvar = 'longitude'
timevar = 'time'

verbose = 0

def main():

    df = pd.read_table(windfile, sep='\s+', header=None, names=('lat', 'lon', 'wind'))

    ncfilename = ncdir + '/' + infiletemplate.format(1)
    nc = Dataset(ncfilename, 'r')
    lons = nc.variables[lonvar][:, 0]
    lats = nc.variables[latvar][0, :]
    nc.close()
    lonres = (lons.max() - lons.min()) / (lons.shape[0]-1)
    latres = (lats.max() - lats.min()) / (lats.shape[0]-1)
    if lonres != latres:
        print 'Internal error: Dimension mismatch'
        raise RuntimeError

    for sample in range(1, nsamples+1, 1):
#    for sample in range(59, 60, 1):
        outdir = ascdir.format(sample)
        ncfilename = ncdir + '/' + infiletemplate.format(sample)
        print ncfilename
        #print outdir
        nc = Dataset(ncfilename, 'r')
        prcp = nc.variables['pcp'][:]
        tmean = nc.variables['t_mean'][:]
        trange = nc.variables['t_range'][:]
        for i in range(len(df.index)):
            x = int((df['lon'][i] - lons.min()) / lonres)
            #print '\nx: ', x, df['lon'][i]
            if lons[x] != df['lon'][i]:
                print 'Internal error: Dimension mismatch'
                raise RuntimeError
            y = int((df['lat'][i] - lats.min()) / latres)
            #print 'y: ', y, df['lat'][i]
            if lats[y] != df['lat'][i]:
                print 'Internal error: Dimension mismatch'
                raise RuntimeError
            cellprcp = prcp[:,x,y]
            celltmax = tmean[:,x,y] + 0.5 * trange[:,x,y]
            celltmin = tmean[:,x,y] - 0.5 * trange[:,x,y]
            cellwind = np.ones(len(celltmax)) * df['wind'][i]
            #print cellprcp.mean(), celltmax.mean(), celltmin.mean(), cellwind.mean()
            outfile = outdir + '/' + outfiletemplate.format(df['lat'][i], df['lon'][i])
            #print outfile
            np.savetxt(outfile,
                       np.transpose([cellprcp, celltmax, celltmin, cellwind]),
                       fmt='%.2f', delimiter=' ', newline='\n')
        nc.close()

if __name__ == "__main__":
    main()
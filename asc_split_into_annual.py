#!/usr/bin/env python
import numpy as np
import sys
import os
import glob

if len(sys.argv) < 4:
    print 'Must provide inpath, outpath, and filemask'
    sys.exit()
    
inpath=sys.argv[1]
outpath=sys.argv[2]
filemask=sys.argv[3]
filelist = glob.glob('{}/{}*'.format(inpath, filemask))

# all the files are the same, so we need to determine the indices only once
infile = filelist[0]
data = np.loadtxt(infile,
                  dtype={'names': ('year', 'month', 'day', 'perc'),
                         'formats': ('i4', 'i4', 'i4', 'f4')})
startyear = data['year'][0]
endyear = data['year'][-1]
year = startyear

# determine the number of months
nyears = endyear - startyear + 1

# initialize empty arrays
idx = np.empty((data['year'].size, nyears), dtype='bool')*False
years = np.ones(nyears, dtype='i4')

# get the indices, year and months
year = startyear
for i in range(nyears):
    idx[:,i] = data['year'] == year
    years[i] = year
    year += 1

# loop over all the files
for infile in filelist:
    print infile
    data = np.loadtxt(infile, dtype='S')
    filename = os.path.basename(infile)
    for i in range(nyears):
        x = data[idx[:,i]]
        np.savetxt('{}/{:04d}/{}'.format(outpath, years[i], filename), x, '%s')


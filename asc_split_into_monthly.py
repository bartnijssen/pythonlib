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
startmonth = data['month'][0]
endyear = data['year'][-1]
endmonth = data['month'][-1]
year = startyear
month = startmonth

# determine the number of months
nmonths = 0
while year*100+month <= endyear*100+endmonth:
    month += 1
    nmonths += 1
    if month > 12:
        year += 1
        month = 1

# initialize empty arrays
idx = np.empty((data['year'].size, nmonths), dtype='bool')*False
months = np.ones(nmonths, dtype='i4')
years = np.ones(nmonths, dtype='i4')

# get the indices, year and months
year = startyear
month = startmonth
for i in range(nmonths):
    idx[:,i] = (data['year'] == year) & (data['month'] == month)
    months[i] = month
    years[i] = year
    month += 1
    if month > 12:
        year += 1
        month = 1

# loop over all the files
for infile in filelist:
    print infile
    data = np.loadtxt(infile, dtype='S')
    filename = os.path.basename(infile)
    for i in range(nmonths):
        x = data[idx[:,i]]
        np.savetxt('{}/{:04d}_{:02d}/{}'.format(outpath, years[i], months[i],
                                                filename),
                    x, '%s')


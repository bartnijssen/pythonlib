#!/usr/bin/env python
import numpy as np
import sys
import os
import glob
import pandas as pd
from utils.parse import parse_date

if len(sys.argv) < 6:
    print 'Must provide inpath, outpath, filemask, startdate and enddate'
    sys.exit()

inpath = sys.argv[1]
outpath = sys.argv[2]
filemask = sys.argv[3]
startdate = sys.argv[4]
enddate = sys.argv[5]

start = parse_date(startdate, '/')
end = parse_date(enddate, '/')
dates = pd.date_range(start, end, freq='D')

filelist = glob.glob('{}/{}*'.format(inpath, filemask))

# all the files are the same, so we need to determine the indices only once
infile = filelist[0]
data = np.loadtxt(infile,
                  dtype={'names': ('prcp', 'tmax', 'tmin', 'wind'),
                         'formats': ('f4', 'f4', 'f4', 'f4')})
# determine the number of months
nyears = end.year - start.year + 1

# initialize empty arrays
idx = np.empty((data['prcp'].size, nyears), dtype='bool')*False
years = np.ones(nyears, dtype='i4')

# get the indices, year and months
year = start.year
for i in range(nyears):
    idx[:, i] = dates.year == year
    years[i] = year
    year += 1

# loop over all the files
for infile in filelist:
    # print infile
    data = np.loadtxt(infile, dtype='S')
    filename = os.path.basename(infile)
    for i in range(nyears):
        x = data[idx[:, i]]
        np.savetxt('{}/{:04d}/{}'.format(outpath, years[i], filename), x, '%s')

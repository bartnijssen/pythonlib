import argparse
import ConfigParser
import datetime as dt
import numpy as np
import os
import pygrib
import sys

# modify the fields between MODIFY START and MODIFY END

# MODIFY START

# indices for summing the precipitation totals for each lead time
# note that in python the index of the first element is 0
# from the existing grib2 files (index - forecast interval):
# 00 -- 0
# 01 -- 0-3
# 02 -- 0-6
# 03 -- 6-9
# 04 -- 6-12
# 05 -- 12-15
# 06 -- 12-18
# 07 -- 18-21
# 08 -- 18-24
# 09 -- 24-27
# 10 -- 24-30
# 11 -- 30-33
# 12 -- 30-36
# 13 -- 36-39
# 14 -- 36-42
# 15 -- 42-45
# 16 -- 42-48
# 17 -- 48-51
# 18 -- 48-54
# 19 -- 54-57
# 20 -- 54-60
# 21 -- 60-63
# 22 -- 60-66
# 23 -- 66-69
# 24 -- 66-72
# 25 -- 72-78
# 26 -- 78-84
# 27 -- 84-90
# 28 -- 90-96
# 29 -- 96-102
# 30 -- 102-108
# 31 -- 108-114
# 32 -- 114-120
# 33 -- 120-126
# 34 -- 126-132
# 35 -- 132-138
# 36 -- 138-144
# 37 -- 144-150
# 38 -- 150-156
# 39 -- 156-162
# 40 -- 162-168
# 41 -- 168-174
# 42 -- 174-180
# 43 -- 180-186
# 44 -- 186-192

day_idx = [
 [2, 4, 6, 8],     # lead time 1
 [10, 12, 14, 16], # lead time 2
 [18, 20, 22, 24], # lead time 3
 [25, 26, 27, 28], # lead time 4
 [29, 30, 31, 32], # lead time 5
 [33, 34, 35, 36], # lead time 6
 [37, 38, 39, 40], # lead time 7
 ]

# MODIFY END

verbose = False

def parse_commandline():
    parser = argparse.ArgumentParser(description='Convert NCEP GFS grib2 '
                                     'forecast files for PNNL forecast'
                                     'project')
    parser.add_argument('--config', required=True,
                        metavar='<configuration file>',
                        help='Configuration file')
    parser.add_argument('--verbose', '-v', action='store_true',
                        default=False, help='verbose')

    args = parser.parse_args()

    global verbose
    if args.verbose:
        verbose = True

    configparser = ConfigParser.ConfigParser()
    configparser.read(os.path.expanduser(args.config))
    config = {}
    for section in configparser.sections():
        for key, value in configparser.items(section):
            config[key] = value

    config['multiplier'] = int(config['multiplier'])
    config['start_date'] = parse_date(config['start_date'], '-')
    config['end_date'] = parse_date(config['end_date'], '-')
    config['ensemble'] = [x.strip() for x in config['ensemble'].split(',')]

    return config

def parse_date(datestr, sep='-'):
    try:
        date = dt.date(*[int(x) for x in datestr.split(sep)])
    except (ValueError, TypeError):
        sys.exit('Not a valid date: {}'.format(datestr))
    except:
        sys.exit('Unexpected error: {}'.format(sys.exc_info()[0]))
    return date

def readmaskfile(filename):
    cells = []
    with open(filename, 'r') as f:
        for line in f:
            cells.append([float(i) for i in line.split()])
    return cells

def read_grib2(filename, parameter, mask, accumulator=None):
    var = []
    grbs = pygrib.open(filename)
    try:
        if grbs[1].parameterName != parameter:
            raise ValueError
    except:
        sys.exit('Grib parameter \'{}\' does not match \'{}\''.
                 format(grbs[1].parameterName, parameter))
    for grb in grbs:
        var.append(grb.values[mask[:,0],mask[:,1]])
    var = np.array(var)
    lats, lons = grb.latlons()
    if accumulator:
        tmpvar = []
        for subset in accumulator:
            tmpvar.append(var[subset, :].sum(axis=0))
        var = np.array(tmpvar)
    grbs.close()
    return var

def read_grib2_latlons(filename):
    grbs = pygrib.open(filename)
    lats, lons = grbs[1].latlons()
    grbs.close()
    return np.array(lats), np.array(lons)

def daterange(start_date, end_date):
    for n in range(int ((end_date - start_date).days)):
        yield start_date + dt.timedelta(n)

def main():

    config = parse_commandline()
    mask_file = config['mask_file']
    start_date = config['start_date']
    end_date = config['end_date']
    out_path = config['out_path']
    grib2_path = config['grib2_path']
    grib2_filename_template = config['grib2_filename_template']
    grib2_parameter_name = config['grib2_parameter_name']
    ensemble = config['ensemble']
    multiplier = config['multiplier']

    # read the mask file which has a list of latitude, longitude pairs that
    # need to be extracted from the grib file. Note that the program only
    # looks for exact matches between the coordinates in the mask file and
    # in the grid file. This means that the longitudes have to be the same.
    # So if the longitudes in the grib file go from 0 to 360, then the mask
    # file needs to show the same.
    mask = readmaskfile(mask_file)

    # read the latitudes and longitudes from the first grib file. The
    # program assumes that ALL grib files have the same grid. This is
    # not checked for each file independently.
    now = start_date
    grib2_filename = os.path.join(grib2_path,
                                  grib2_filename_template.format(now.year,
                                                                 now.month,
                                                                 now.day,
                                                                 ensemble[0]))
    lats, lons = read_grib2_latlons(grib2_filename)

    # for each location in the mask file, determine the corresponding
    # indices in the grib file (mask_idx) and set up the arrays to contain
    # the reforecast data (fcst)
    mask_idx = []
    fcst = []
    # total number of days in period
    ndays = (end_date - start_date).days+1
    for location in mask:
        mask_idx.append(np.where((lats==location[0]) &
                                 (lons==location[1])))
        fcst.append(np.zeros([len(day_idx), ndays, len(ensemble)],
                             dtype=np.uint16))
    mask_idx = np.squeeze(mask_idx)

    # loop over the days and the ensemble members and read the grib files.
    # For each location store the results in fcst. Note that the grib file
    # values are summed to days in the read_grib2() function by using the
    # day_idx values
    for now in daterange(start_date, end_date + dt.timedelta(1)):
        # day of year (0-36[45])
        yday = now.timetuple().tm_yday-1
        # day since start of record (0-ndays-1)
        iday = (now - start_date).days
        # year since start of record
        year_idx = now.year - start_date.year
        for ens_idx, member in enumerate(ensemble):
            grib2_filename = os.path.join(grib2_path,
                                          grib2_filename_template.format(now.year,
                                                                         now.month,
                                                                         now.day,
                                                                         member))
            if verbose:
                print 'Processing: {}'.format(grib2_filename)

            var = read_grib2(grib2_filename, grib2_parameter_name,
                             mask_idx, day_idx)
            var = np.uint16(var * multiplier)
            for loc_idx in range(len(fcst)):
                fcst[loc_idx][:,iday,ens_idx] = var[:,loc_idx]

    # write the resulting arrays to a binary file, with one file per location
    for i, location in enumerate(mask):
        outfilename = os.path.join(out_path,
                                   'reforecast_{:.1f}_{:.1f}'.
                                   format(location[0], location[1]))
        fcst[i].tofile(outfilename)

if __name__ == "__main__":
    main()



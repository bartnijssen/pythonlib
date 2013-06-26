#!/usr/bin/env python
"""Postprocess the daily VIC files from RASM simulations

The script does the following to a series of daily RASM VIC files:
 * Fixes the time stamp of the daily VIC files
 * Calculates monthly and annual averages
 * Concatenates the daily files into a annual file that can be used for further
   processing
 * Fixes the time axis and the time attributes of the annual files so that
   ncview displays the correct date
 * Deletes the original daily files to preserve space. The original files are
   only deleted if all other steps complete successfully.

Note that any existing file will be overwritten.

The script relies on the nco utilities for some of its steps: http://nco.sourceforge.net
"""

import argparse
import calendar
import datetime
import glob
import os
import sys
import subprocess

ncra = '/opt/local/bin/ncra'
ncrcat = '/opt/local/bin/ncrcat'
ncap = '/opt/local/bin/ncap'
ncatted = '/opt/local/bin/ncatted'


def main():
    parser = argparse.ArgumentParser(description='Postprocess the daily VIC '
                                     'files from RASM simulations')
    parser.add_argument('--startdate', required=True, metavar='<start date in YYYY-MM-DD>',
                        help='Date stamp of first VIC file from simulation')
    parser.add_argument('--enddate', required=True, metavar='<end date in YYYY-MM-DD>',
                        help='Date stamp of last VIC file from simulation')
    parser.add_argument('--srcdir', required=True, metavar='<source directory>',
                        help='Directory with RASM VIC files')
    parser.add_argument('--destdir', required=True, metavar='<target directory>',
                        help='Directory with processed RASM VIC files')
    parser.add_argument('--casename', required=True, metavar='<casename>',
                        help='RASM case name')
    parser.add_argument('--cal', required=True, metavar='<leap|noleap>',
                        help='calendar includes leap days or not')
    parser.add_argument('--delete', action='store_true',
                        help='delete files that have been updated or merged')

    args = parser.parse_args()

    startdate = parse_date(args.startdate)
    enddate = parse_date(args.enddate)
    leapcal = args.cal.lower()
    if leapcal not in ['leap', 'noleap']:
        sys.exit('Not a valid calendar: {}'.format(args.calendar))

    if args.srcdir == args.destdir:
        sys.exit('Source directory and target directory must be different')
    if not os.path.exists(args.srcdir):
        sys.exit('Source path does not exist: {}'.format(args.srcdir))
    if not os.path.exists(args.destdir):
        os.makedirs(args.destdir)

    tobedeleted = []
    # Change the time stamp of the daily files and fix time axis
    filedate = startdate
    newfiledate = filedate + datetime.timedelta(days=-1)
    while filedate <= enddate:
        if leapcal == 'noleap' and filedate.month == 2 and filedate.day == 29:
            filedate += datetime.timedelta(days=1)
        filename = ('{0}.vic.ha.{1:04d}-{2:02d}-{3:02d}.nc'.
                    format(args.casename, filedate.year, filedate.month,
                           filedate.day))
        newfilename = ('{0}.vic.ha.{1:04d}-{2:02d}-{3:02d}.nc'.
                       format(args.casename, newfiledate.year,
                              newfiledate.month, newfiledate.day))
        print '{:40s} -> {:40s}'.format(filename, newfilename)
        fargs = [ncap, '-O', '-s', 'time=time-366']
        fargs.append(os.path.join(args.srcdir, filename))
        fargs.append(os.path.join(args.destdir, newfilename))
        returnval = subprocess.call(fargs)
        if returnval != 0:
            sys.exit('Error executing: {}'.format(' '.join(fargs)))
        fargs = [ncatted]
        fargs.append(['-a', 'units,time,a,c,days since 1-1-1'])
        fargs.append(['-a', 'long_name,time,a,c,time'])
        fargs.append(['-a', 'dimensions,time,a,c,1'])
        fargs.append(['-a', 'calendar,time,a,c,noleap'])
        fargs.append(['-a', 'type_preferred,time,a,c,int'])
        fargs.append(['-a', 'title,global,o,c,{}'.format(newfilename)])
        fargs.append(os.path.join(args.destdir, newfilename))
        fargs = list(flatten(fargs))
        returnval = subprocess.call(fargs)
        if returnval != 0:
            sys.exit('Error executing: {}'.format(' '.join(fargs)))
        tobedeleted.append(os.path.join(args.srcdir, filename))
        newfiledate = filedate
        filedate += datetime.timedelta(days=1)

    # Create monthly and annual average using ncra (note that only whole months
    # and years will be averaged)
    firstdate = startdate + datetime.timedelta(days=-1)
    if (firstdate.day != 1):
        day = 1
        month = firstdate.month + 1
        year = firstdate.year
        if month > 12:
            month = 1
            year += 1
        firstdate = datetime.date(year, month, day)
    lastdate = enddate + datetime.timedelta(days=-1)
    day = lastdate.day
    daysinmonth = calendar.monthrange(lastdate.year, lastdate.month)[1]
    if leapcal == 'noleap' and lastdate.month == 2:
        daysinmonth = 28
        day = 28
    if day != daysinmonth:
        month = lastdate.month - 1
        year = lastdate.year
        if month < 1:
            year -= 1
            month = 12
        day = calendar.monthrange(year, month)[1]
        if leapcal == 'noleap' and month == 2:
            day = 28
        lastdate = datetime.date(year, month, day)

    # monthly means
    currentdate = firstdate
    firstmonth = currentdate
    while currentdate <= lastdate:
        fileglob = ('{0}.vic.ha.{1:04d}-{2:02d}-??.nc'.
                    format(args.casename, currentdate.year, currentdate.month))
        filelist = glob.glob(os.path.join(args.destdir, fileglob))
        outfile = ('{0}.vic.ha.{1:04d}-{2:02d}.mean.nc'.
                   format(args.casename, currentdate.year, currentdate.month))
        fargs = list(flatten([ncra, '-O', filelist,
                              os.path.join(args.destdir, outfile)]))
        print '{} [...] -> {}'.format(fargs[0], fargs[-1])
        returnval = subprocess.call(fargs)
        if returnval != 0:
            sys.exit('Error executing: {}'.format(' '.join(fargs)))
        ndays = calendar.monthrange(currentdate.year, currentdate.month)[1]
        lastmonth = currentdate
        currentdate += datetime.timedelta(days=ndays)

    # time series of monthly means
    fileglob = '{0}.vic.ha.????-??.mean.nc'.format(args.casename)
    filelist = glob.glob(os.path.join(args.destdir, fileglob))
    outfile = ('{0}.vic.ha.{1:04d}{2:02d}-{3:04d}{4:02d}.monthly.mean.nc'.
               format(args.casename, firstmonth.year, firstmonth.month,
                      lastmonth.year, lastmonth.month))
    fargs = list(flatten([ncrcat, '-O', filelist,
                          os.path.join(args.destdir, outfile)]))
    print '{} [...] -> {}'.format(fargs[0], fargs[-1])
    returnval = subprocess.call(fargs)
    if returnval != 0:
        sys.exit('Error executing: {}'.format(' '.join(fargs)))
    tobedeleted.append(filelist)

    # annual means: only complete years
    firstyear = firstdate.year
    lastyear = lastdate.year
    if firstdate.month != 1:
        firstyear += 1
    if lastdate.month != 12:
        lastyear -= 1
    currentyear = firstyear
    while currentyear <= lastyear:
        fileglob = ('{0}.vic.ha.{1:04d}-??-??.nc'.
                    format(args.casename, currentyear))
        filelist = glob.glob(os.path.join(args.destdir, fileglob))
        outfile = ('{0}.vic.ha.{1:04d}.mean.nc'.
                   format(args.casename, currentyear))
        fargs = list(flatten([ncra, '-O', filelist,
                              os.path.join(args.destdir, outfile)]))
        print '{} [...] -> {}'.format(fargs[0], fargs[-1])
        returnval = subprocess.call(fargs)
        if returnval != 0:
            sys.exit('Error executing: {}'.format(' '.join(fargs)))
        currentyear += 1

    # annual time series: incomplete years as well
    firstyear = firstdate.year
    lastyear = lastdate.year
    currentyear = firstyear
    while currentyear <= lastyear:
        fileglob = ('{0}.vic.ha.{1:04d}-??-??.nc'.
                    format(args.casename, currentyear))
        filelist = glob.glob(os.path.join(args.destdir, fileglob))
        outfile = ('{0}.vic.ha.{1:04d}.daily.nc'.
                   format(args.casename, currentyear))
        fargs = list(flatten([ncrcat, '-O', filelist,
                              os.path.join(args.destdir, outfile)]))
        print '{} [...] -> {}'.format(fargs[0], fargs[-1])
        returnval = subprocess.call(fargs)
        if returnval != 0:
            sys.exit('Error executing: {}'.format(' '.join(fargs)))
        tobedeleted.append(filelist)
        currentyear += 1

    # cleanup
    if args.delete is True:
        print 'Cleanup: Deleting original and merged files'
        tobedeleted = list(flatten(tobedeleted))
        for file in tobedeleted:
            try:
                os.remove(file)
            except OSError as e:
                print 'Error removing file: {}'.format(e)
                pass


def flatten(*args):
    for x in args:
        if hasattr(x, '__iter__'):
            for y in flatten(*x):
                yield y
        else:
            yield x


def parse_date(datestr, sep='-'):
    try:
        date = datetime.date(*[int(x) for x in datestr.split(sep)])
    except (ValueError, TypeError):
        sys.exit('Not a valid date: {}'.format(datestr))
    except:
        sys.exit('Unexpected error: {}'.format(sys.exc_info()[0]))
    return date

if __name__ == "__main__":
    main()

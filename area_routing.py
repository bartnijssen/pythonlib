#!/usr/bin/env python
"""Determine the contributing area for each station in a station file

Calculate the upstream contributing area for each station in a station file given a flow direction, flow fraction and station file.

The final output is simply printed to the screen. The area is reported in acres, mi2 and km2.
"""

import argparse
import math
import numpy as np

archeaderlines = 6
arcdirections = {1: 'E', 2: 'SE', 4: 'S', 8: 'SW', 16: 'W', 32: 'NW', 64: 'N', 128: 'NE'}
vicdirections = {1: 'N', 2: 'NE', 3: 'E', 4: 'SE', 5: 'S', 6: 'SW', 7: 'W', 8: 'NW'}
diroffset = {'E': (0,1), 'SE': (1,1), 'S': (1,0), 'SW': (1,-1),
             'W': (0,-1), 'NW': (-1,-1), 'N': (-1,0), 'NE': (-1,1)}
diropposite = {'E': 'W', 'SE': 'NW', 'S': 'N', 'SW': 'NE',
               'W': 'E', 'NW': 'SE', 'N': 'S', 'NE': 'SW'}

class FileError(Exception):
    """Exception raised for errors in the input.

    Attributes:
        filename -- file in which the error occured
        msg      -- explanation of the error
    """

    def __init__(self, filename, msg):
        self.filename = filename
        self.msg = msg

    def __str__(self):
        return repr("Error in \"" + self.filename + "\": " + self.msg)
        
def main():
    parser = argparse.ArgumentParser(description='Calculate contributing area '
                                     'for each station in a station file. All '
                                     'input files have to conform to the '
                                     'specifications of the input files for the '
                                     'VIC routing program. See '
                                     'http://www.hydro.washington.edu/Lettenmaier/'
                                     'Models/VIC/Documentation/Routing/'
                                     'RoutingInput.shtml')
    parser.add_argument('stationfile', metavar='<station file>',
                        help='station file')
    parser.add_argument('directionfile', metavar='<direction file>',
                        help='flow direction file')
    parser.add_argument('fractionfile', metavar='<fractionfile>',
                        help='fraction file')
    
    args = parser.parse_args()

    # Read the files
    stations = parse_stationfile(args.stationfile)
    header = parse_gridasciiheader(args.directionfile)
    directions = parse_gridasciifile(args.directionfile)
    fractions = parse_gridasciifile(args.fractionfile)

    # Calculate the area of each grid cell and multiply it by the fraction
    areas = [(x + 0.5) * float(header['cellsize']) + float(header['yllcorner'])
             for x in range(int(header['nrows']))]
    areas = [calculate_area_vicrouting(x, float(header['cellsize'])) for x in areas]
    # Invert the area array, since the top row of the header and fraction tables
    # are north
    areas = areas[::-1]
    # Create a grid that has the area for each cell. Use the broadcast
    # functionality of numpy for the multiplication. This does require
    # transposing and re-transposing
    areas = (fractions.transpose()*areas).transpose()/1000000.
    areas.tofile('areas.save.txt', sep = ' ')

    # For each station, determine the contributing grid cells and calculate the
    # upstream area
    for id, record in stations.iteritems():
        mask = mask_contributing(record, directions, vicdirections)
        totalarea = (mask * areas).sum()
        print id, ":", totalarea/4046.85642, "acres,", totalarea/2589988.11, "miles^2,", totalarea/1000000., "km^2"
        mask.tofile(id+'.save.txt', sep=' ')
        
def calculate_area_vicrouting(centerlat, resolution):
    """Calculate the cell area according to the VIC routing model"""
    earth_radius = 6371229.0
    resolution = math.radians(resolution)
    centerlat = math.radians(centerlat)
    area = (math.pow(earth_radius, 2) * resolution * 
            abs(math.sin(centerlat - 0.5*resolution) -
                math.sin(centerlat + 0.5*resolution)))
    return area
    
def mask_contributing(record, directions, flowdirections=arcdirections):
    directions = directions.astype(np.int32)
    alloweddirections = flowdirections.keys()
    mask = np.zeros(directions.shape)
    nrows, ncols = mask.shape
    #    print nrows, ncols
    #print record['row'], record['col']
    outrow = nrows - record['row'] # VIC station file counts up from bottom
                                        # with first record as 1
    outcol = record['col'] - 1
    #    print outrow, outcol
    contrib = [];
    upstream = [(outrow, outcol)]
    while upstream:
        (row, col) = upstream.pop(0)
        if ((row,col) in contrib):
            #            print "Loop", (row, col), contrib
            raise FileError("Loop", "Arghhhh")
        #        print "Current cell: ", row, col
        mask[row][col] = 1;
        contrib.append((row, col))
        for dir, offset in diroffset.iteritems():
            orow = row + offset[0]
            ocol = col + offset[1]
            #            print "\tNeighbor cell", dir, offset, orow, ocol
            if (orow >= 0 and orow < nrows and ocol >= 0 and ocol < ncols and
                directions[orow][ocol] in alloweddirections):
                #                print "\t\tDirection: ", directions[orow][ocol], flowdirections[directions[orow][ocol]], diropposite[dir]
                if (flowdirections[directions[orow][ocol]] ==
                    diropposite[dir]):
                    upstream.append((orow, ocol))
                    #                    print "\t\tAppended", orow, ocol
        
    return mask
    
def parse_stationfile(infile):
    """Parse a VIC station file and return a dict of stations"""
    stations = {};
    
    with open(infile, 'r') as f:
        contents = f.readlines()

    # strip whitespace and beginning and end of line. Note that by assigning to
    # a slice of the original list, the reference to that list does not change
    # (this does not really matter here, but is good practice if you are
    # essentially replacing "in-place"
    contents[:] = [line.strip() for line in contents]
    # only keep the non-empty lines
    contents[:] = [line for line in contents if line]

    if len(contents)%2 != 0:
        raise FileError(infile, "Expectation is two lines per station")
    
    # Need to process two lines at a time
    for i in range(0,len(contents),2):
        fields = contents[i].split()
        if (len(fields) != 5):
            raise FileError(infile, "Station line must have five fields")
        station = fields[1]
        stations[station] = {}
        stations[station]['col'] = int(fields[2])
        stations[station]['row'] = int(fields[3])
        stations[station]['uhs'] = contents[i+1]        
    return stations

def parse_gridasciifile(infile):
    header = parse_gridasciiheader(infile)
    # Note that the missing_values part in genfromtxt does not appear to work at
    # all, so I am not relying on it here
    data = np.genfromtxt(infile, skip_header=archeaderlines)
    # Because the missing_values part np.genfromtxt is not working as expected,
    # I do it myself here
    missing = np.fromstring(header['nodata_value'], data.dtype, sep=' ')
    data[data==missing]=np.nan
    if (data.shape != (int(header['nrows']), int(header['ncols']))):
        raise FileError(infile, "Number of rows and columns different from header")
    return data

def parse_gridasciiheader(infile, headerlines=archeaderlines):
    """Parse the header of an ArcInfo grid ascii file and return as dict"""
    contents = []
    with open(infile, 'r') as f:
        for i in range(archeaderlines):
            contents.append(f.readline())
    contents[:] = [line.strip() for line in contents]
    return dict([line.lower().split() for line in contents])

    
    
    
if __name__ == "__main__":
    main()

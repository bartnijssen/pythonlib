#!/usr/bin/env python
"""Given a HU-12 code find all the HU-12 elements that flow into it

Find all the HU-12 elements that contribute to a given HU-12 unit. This is
determined by following the HU_12_DS field in the shapefiles downstream. An
extra field CONTRIBUTES is added that is either 1 or 0 (true or false). The
output is saved as a new file with the name of the most downstream element.

The final output is saved as both a shapefile and a KML file. Note that all
shapefiles are combined into one new one. The original files are not
overwritten.
"""

import argparse
import glob
import shapefile
import subprocess
import sys

parser = argparse.ArgumentParser(description='Find all units that contribute to '
                                 'a given HU-12 unit in WBC HUC12 shapefiles')
parser.add_argument('--indir', required=True,
                    help='input directory with WBD HUC12 shapefiles')
parser.add_argument('--outdir', required=True,
                    help='output directory for new WBD HUC12 shapefile')
parser.add_argument('--id', required=True,
                    help='HUC12 ID of the most downstream element. Files will'
                    'be selected from <indir> based on this ID')
args = parser.parse_args()

if len(args.id) != 12:
    parser.error("Invalid ID: \"{}\" -- Needs to have 12 digits".
                 format(args.id))

# I am simply going to assume that the first two digits determine which files I
# need to read in. In almost all cases this will be an area that is much too
# large, but so be it. This can be improved as need be
infilenames = glob.glob(args.indir + '/wbdhu12_a_' + args.id[:2] + '*' + '.shp')

# hus is a dict with the hu_12 id as the key. 
hus = {}
shapes = {}
for infile in infilenames:
    sf = shapefile.Reader(infile)
    for i in range(sf.numRecords):
        record = sf.record(i)
        shape = sf.shape(i)
        record.append(False)
        if (record[2] == args.id):
            record[-1] = True
        hus[record[2]] = record
        shapes[record[2]] = shape

if args.id not in hus:
    print "ID: {} not found in any of the files:".format(args.id)
    for filename in infilenames:
        print "\t{}".format(filename)
    sys.exit(1)
    
fields = sf.fields[1:]
fields.append(['CONTRIBUTES', 'N', 12, 0])

for hu, record in hus.iteritems():
    # follow the hus downstream
    while record[11] in hus and record[11] != record[2] and record[11] != args.id:
        record = hus[record[11]]
    if record[11] == args.id:
        hus[hu][-1] = True

# retain only the upstream elements
hus = {hu : record for hu, record in hus.iteritems() if record[-1] == True}
shapes = {hu : shape for hu, shape in shapes.iteritems() if hu in hus}
for hu, record in hus.iteritems():
    record[-1] = args.id

# Write a new shapefile with all the elements
w = shapefile.Writer(shapeType=sf.shapeType)
w.autoBalance = 1
for field in fields:
    w.field(*field)

for hu in sorted(hus.keys()):
    w.poly(shapeType=shapes[hu].shapeType, parts=[shapes[hu].points])
    w.record(*(hus[hu]))

w.save(args.outdir + '/' + args.id)

# Translate to kml using ogr2ogr
subprocess.call(['ogr2ogr', '-f', 'KML',
                 args.outdir + '/' + args.id + '.kml',
                 args.outdir + '/' + args.id + '.shp'])

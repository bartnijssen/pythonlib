#!/usr/bin/env python
"""Calculate upstream area for HUC12 in WBD

Add an additional field to the HUC12 shapefiles to account for the total
upstream area. This is calculated by following the HU_12_DS field in the
shapefiles downstream.

The final output is saved as both a shapefile and a KML file. Note that all
shapefiles are combined into one new one. The original files are not
overwritten.
"""

import argparse
import glob
import shapefile
import subprocess
import sys

parser = argparse.ArgumentParser(description='Accumulate upstream area in WBC '
                                 'HUC12 shapefiles')
parser.add_argument('--indir', required=True,
                    help='input directory with WBD HUC12 shapefiles')
parser.add_argument('--outdir', required=True,
                    help='output directory for new WBD HUC12 shapefiles')
parser.add_argument('--id', required=True,
                    help='HUC ID of the region to be processed. Files will '
                    'be selected from <indir> based on this ID')
args = parser.parse_args()

infilenames = glob.glob(args.indir + '/wbdhu12_a_' + args.id + '*' + '.shp')

# hus is a dict with the hu_12 id as the key. The upstream area is appended as a
# new element
hus = {}
shapes = {}
for infile in infilenames:
    sf = shapefile.Reader(infile)
    for i in range(sf.numRecords):
        record = sf.record(i)
        shape = sf.shape(i)
        # area of the hu minus the non-contributing area
        area = record[3] - record[4]
        record.append(area)
        hus[record[2]] = record
        shapes[record[2]] = shape
        
fields = sf.fields[1:]
fields.append(['ACCAREA', 'N', 20, 1])

# outlets is dict of all terminal hu_12 ids. Mainly used for checking for
# errors
outlets = {}
for hu, record in hus.iteritems():
    area = record[3] - record[4]
    # follow the hus downstream
    while record[11] in hus and record[11] != record[2]:
        record = hus[record[11]]
        # add the area of the hu to all downstream hus
        record[-1] += area
    if record[11] not in outlets:
        outlets[record[11]] = 0
    outlets[record[11]] += 1

print 'Unique outlets:'
for outlet, n in outlets.iteritems():
    print outlet, ': ', n
    for hu, record in hus.iteritems():
        if record[11] == outlet:
            print '\t', record[2]

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

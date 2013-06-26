#!/usr/bin/env python
import argparse
import ConfigParser
import os
import vic.assimilate as va
from utils.parse import parse_date


def parse_commandline():
    parser = argparse.ArgumentParser(description='Update VIC state file '
                                     'based on assimilation of observed '
                                     'streamflow')
    parser.add_argument('--date', required=True,
                        metavar='<date>',
                        help='Forecast date: YYYY-MM-DD')
    parser.add_argument('--stations', required=True,
                        metavar='<stations>',
                        help='Comma-delimited list of stations to use '
                        'in the update')
    parser.add_argument('--config', required=True,
                        metavar='<configuration file>',
                        help='Configuration file')

    args = parser.parse_args()

    forecastdate = parse_date(args.date, '-')

    configparser = ConfigParser.ConfigParser()
    configparser.read(os.path.expanduser(args.config))
    config = {}
    for section in configparser.sections():
        for key, value in configparser.items(section):
            config[key] = value

    config['timeofconcentration'] = int(config['timeofconcentration'])
    config['latlondigits'] = int(config['latlondigits'])

    stations = dict((station, config['timeofconcentration'])
                    for station in args.stations.split(','))

    return (forecastdate, stations, config)


def main():

    (forecastdate, stations, config) = parse_commandline()

    stationinfo = va.assimilate_streamflow(forecastdate=forecastdate,
                                           stationdict=stations,
                                           statetemplate=config['statetemplate'],
                                           soilfile=config['soilfile'],
                                           masktemplate=config['masktemplate'],
                                           fluxtemplate=config['fluxtemplate'],
                                           vicflowtemplate=config['vicflowtemplate'],
                                           obsflowtemplate=config['obsflowtemplate'],
                                           polytemplate=config['polytemplate'],
                                           latlondigits=config['latlondigits'])
        #    for station in sorted(stationinfo.iterkeys()):
        #        print '{}:'.format(station)
        #        for k,v in stationinfo[station].iteritems():
        #            print '\t{} : {}'.format(k, v)

if __name__ == "__main__":
    main()

#!/usr/bin/env python
import argparse
import ConfigParser
import os
import vic.assimilate as va
from utils.parse import parse_date


def parse_commandline():
    parser = argparse.ArgumentParser(description='Determine polynomial '
                                     'relating streamflow and soil '
                                     'moisture for VIC streamflow '
                                     'assimilation')
    parser.add_argument('--station', required=True,
                        metavar='<station>',
                        help='station')
    parser.add_argument('--config', required=True,
                        metavar='<configuration file>',
                        help='Configuration file')

    args = parser.parse_args()

    configparser = ConfigParser.ConfigParser()
    configparser.read(os.path.expanduser(args.config))
    config = {}
    for section in configparser.sections():
        for key, value in configparser.items(section):
            config[key] = value

    config['timeofconcentration'] = int(config['timeofconcentration'])
    config['startwindow'] = parse_date(config['startwindow'], '-')
    config['endwindow'] = parse_date(config['endwindow'], '-')
    if config['figuretemplate'] == 'None':
        config['figuretemplate'] = None

    return (args.station, config)


def main():
    (station, config) = parse_commandline()
    p = va.determine_update_function(station=station,
                                     timeofconcentration = config['timeofconcentration'],
                                     startwindow=config['startwindow'],
                                     endwindow=config['endwindow'],
                                     masktemplate=config['masktemplate'],
                                     fluxtemplate=config['fluxtemplate'],
                                     vicflowtemplate=config['vicflowtemplate'],
                                     polytemplate=config['polytemplate'],
                                     figuretemplate=config['figuretemplate'],
                                     latlondigits=config['latlondigits'])


if __name__ == "__main__":
    main()

# This code is to support a particle filter implementation using VIC
from __future__ import division
import collections
import ConfigParser
import datetime as dt
#import matplotlib as mpl
#import matplotlib.pyplot as plt
import multiprocessing
import numpy as np
#import pandas as pd
import re
import subprocess
#import sys
from utils.general import flatten
from utils.parse import parse_date
#import uuid

class Error(Exception):
    """Base class for exceptions in this module"""
    pass

class SimulationError(Error):
    """Exception raised for simulation errors"""
    pass

class Simulation(object):
    """Model simulation"""
    def __init__(self, exe=None, configfile=None, args=None, logfile=None):
        """Initiate a simulation instance"""
        self.exe = exe
        self.args = args
        self.configfile = configfile
        self.logfile = logfile

    def __str__(self):
        bstr = 'exe: {}\n'.format(self.exe)
        bstr += 'arguments:\n'
        if self.args:
            bstr += ', '.join(self.args)
        bstr += '\nconfigfile: {}'.format(self.configfile)
        bstr += '\nlogfile: {}'.format(self.logfile)
        return bstr

    def __repr__(self):
        bstr = '{}'.format(self.exe)
        if self.args:
            bstr += ': '
            bstr += ', '.join(self.args)
        if self.logfile:
            bstr += '> {}'.format(self.logfile)
        return bstr

    def run(self):
        """Run the simulation"""
        if self.exe is None:
            raise SimulationError('No executable defined')
        fargs = [self.exe]
        if self.args:
            fargs.append(self.args)
        fargs = list(flatten(fargs))
        loghandle = None
        if type(self.logfile) is str:
            loghandle = open(self.logfile, 'w')
        elif type(self.logfile) is file:
            loghandle = self.logfile
        elif self.logfile is not None:
            raise TypeError('logfile must be a filename, filehandle, or None')
        returnval = subprocess.call(fargs, stdout=loghandle, stderr=loghandle)
        if returnval != 0:
            raise SimulationError('return value not zero')
        if type(self.logfile) is str:
            loghandle.close()

class VicSimulation(Simulation):
    """VIC simulation class"""
    def writeconfig(self, config):
        repl = re.compile('^OUTVAR\d*$')
        with open(self.configfile, 'w') as f:
            for key, value in config.items():
                f.write('{:20s} {}\n'.format(re.sub(repl, 'OUTVAR', key),
                        value))


class RouteSimulation(Simulation):
    """route simulation class"""
    pass

class ParticleFilter(object):
    """Particle filter"""
    def __init__(self, configurationfile):
        """Initialize a particle filter based on a configurationfile"""
        self.pool = None
        # Parse the configuration file
        configparser = ConfigParser.SafeConfigParser(allow_no_value=True)
        configparser.optionxform = str # preserve case of configuration keys
        print 'Reading {}'.format(configurationfile)
        configparser.read(configurationfile)
        self.config = collections.OrderedDict()
        for section in configparser.sections():
            self.config[section] = collections.OrderedDict()
            for key, value in configparser.items(section):
                self.config[section][key] = value
        if not self.config['COMPUTING']['poolsize']:
            self.poolsize = multiprocessing.cpu_count()
        else:
            self.poolsize = configparser.getint('COMPUTING','poolsize')

        # initialize simulation time
        self.start = parse_date(self.config['SIMULATION']['start'])
        self.now = self.start
        self.end = parse_date(self.config['SIMULATION']['end'])
        self.timeinterval = dt.timedelta(configparser.getint('SIMULATION',
                                                             'interval'))

        # Initialize particles
        self.np = configparser.getint('FILTER', 'np')
        self.particles = [VicParticle(i+1, 1./self.np) for i in range(self.np)]
        for i in range(self.np):
            p = self.particles[i]
            p.vic['exe'] = self.config['PATHS']['vicexe']
            p.vic['infile'] = self.config['PATHS']['vicinfiletemplate'].format(i+1)
            p.vic['logfile'] = self.config['PATHS']['viclogfiletemplate'].format(i+1)
            p.route['exe'] = self.config['PATHS']['routeexe']
            p.route['infile'] = self.config['PATHS']['routeinfiletemplate'].format(i+1)
            p.route['logfile'] = self.config['PATHS']['routelogfiletemplate'].format(i+1)
        self.normalize_weights()

    def __repr__(self):
        bstr = '{}\n'.format(self.config)
        bstr += '{}: {}\n'.format(self.poolsize, self.pool)
        bstr += '{}: {}\n'.format(self.np, [self.particles[i] for i in range(self.np)])
        bstr += '{} -- {}({}) -- {}\n'.format(self.start, self.now,
                                              self.timeinterval, self.end)
        return bstr

    def advance(self):
        # Advance all the particles
        self.config['VICCONFIG']['STARTYEAR'] = self.now.year
        self.config['VICCONFIG']['STARTMONTH'] = self.now.month
        self.config['VICCONFIG']['STARTDAY'] = self.now.day
        endsimulation = self.now + self.timeinterval
        self.config['VICCONFIG']['ENDYEAR'] = endsimulation.year
        self.config['VICCONFIG']['ENDMONTH'] = endsimulation.month
        self.config['VICCONFIG']['ENDDAY'] = endsimulation.day

        self.pool = multiprocessing.Pool(self.poolsize)
        results = [self.pool.apply_async(advance, [self.particles[i], self.config['VICCONFIG']])
                   for i in range(self.np)]
        self.pool.close()
        self.pool.join()
        # Report errors
        for r in results:
            if r.get() is not None:
                raise(SimulationError(r.get()))
        # Update the simulation time
        self.now += self.timeinterval

    def evaluate(self):
        for p in self.particles:
            p.evaluate()

    def normalize_weights(self):
        total = np.array([self.particles[i].weight for i in range(self.np)]).sum()
        for p in self.particles:
            p.weight /= total

    def resample(self):
        pass

class Particle(object):
    """Single particle in a particular filter"""

    def __init__(self, ID, weight, state=None):
        """Initialize a particle"""
        self.ID = ID
        self.weight = weight
        self.state = state
        print "initialized particle {}".format(self.ID)

    def __repr__(self):
        bstr = '{}({}): {}'.format(self.ID, self.weight, self.state)
        return bstr

    def advance(self):
        """Advance a particle"""
        pass

    def evaluate(self):
        """Evaluate a particle (calculate its weight)"""
        pass

class VicParticle(Particle):
    """Particle for running VIC simulations"""

    def __init__(self, ID, weight, state=None):
        # invoke the parent's initializer
        super(VicParticle, self).__init__(ID, weight, state)
        self.vic = {'exe': None, 'infile': None, 'logfile': None}
        self.route = {'exe': None, 'infile': None, 'logfile': None}

    def advance(self, config):
        """Advance a particle"""
        self.vicsim = VicSimulation(self.vic['exe'], self.vic['infile'],
                                    ['-g', self.vic['infile']],
                                    self.vic['logfile'])
        self.vicsim.writeconfig(config)
        # self.vicsim.run()
        # self.routesim = RouteSimulation(rout['exe'], [route['infile'], route['logfile']).run()
        # self.routesim.run()


def advance(particle, config):
    particle.advance(config)

if __name__ == '__main__':
    pf = ParticleFilter('/d1/nijssen/vicpf/config/gunnison/particlefilter/gunnison_pf.cfg')
    print pf
    print "Advancing filter 1"
    pf.advance()
#    print pf
#    print "Advancing filter 2"
#    pf.advance()
#    print pf

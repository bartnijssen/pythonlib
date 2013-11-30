# This code is to support a particle filter implementation using VIC
from __future__ import division
import collections
import ConfigParser
import copy
import datetime as dt
#import matplotlib as mpl
#import matplotlib.pyplot as plt
import glob
import multiprocessing
import netCDF4
import numpy as np
import os
import pandas as pd
import re
import subprocess
import sys
from utils.general import flatten, make_sure_path_exists
from utils.parse import parse_date
from vic.vicasc2routenc import vicasc2routenc
from vic.vicasc2nc import vicasc2nc, filterPick
#import uuid


def run(particle):
    particle.run()

def convertvicasc2route(vicconfig, routeconfig, vicasc2routencconfig):
    inpath = vicconfig['RESULT_DIR']
    if inpath[-1] != '/':
        inpath += '/'
    filelist = glob.glob(inpath + vicasc2routencconfig['vicfilefilter'] + "*")
    domainfile = routeconfig['SECTION|DOMAIN|FILE_NAME']
    columns = {}
    columns['runoff'] = vicasc2routencconfig['runoff_column']
    columns['baseflow'] = vicasc2routencconfig['baseflow_column']
    period = {}
    period['start'] = dt.datetime(vicconfig['STARTYEAR'],
                                  vicconfig['STARTMONTH'],
                                  vicconfig['STARTDAY'])
    period['end'] = dt.datetime(vicconfig['ENDYEAR'],
                                vicconfig['ENDMONTH'],
                                vicconfig['ENDDAY'])
    period['interval'] = dt.timedelta(1)
    outfiletemplate = '{}/{}'.format(
                        routeconfig['SECTION|INPUT_FORCINGS|DATL_PATH'],
                        routeconfig['SECTION|INPUT_FORCINGS|DATL_FILE'])
    ncvarinfo = {}
    ncvarinfo['name'] = vicasc2routencconfig['ncvar_name']
    ncvarinfo['long_name'] = vicasc2routencconfig['ncvar_longname']
    ncvarinfo['units'] = vicasc2routencconfig['ncvar_units']
    ncvarinfo['divideby'] = float(vicasc2routencconfig['ncvar_divideby'])
    ncformat = vicasc2routencconfig['ncformat']
    print 'Converting output: {}'.format(outfiletemplate)
    vicasc2routenc(filelist, domainfile, columns, period, outfiletemplate,
                   ncvarinfo, ncformat)

def convertvicasc2nc(vicconfig, domainfile, vicasc2ncconfig):
    inpath = vicconfig['RESULT_DIR']
    if inpath[-1] != '/':
        inpath += '/'
    filelist = glob.glob(inpath + vicasc2ncconfig['vicfilefilter'] + "*")
    period = {}
    period['start'] = dt.datetime(vicconfig['STARTYEAR'],
                                  vicconfig['STARTMONTH'],
                                  vicconfig['STARTDAY'])
    period['end'] = dt.datetime(vicconfig['ENDYEAR'],
                                vicconfig['ENDMONTH'],
                                vicconfig['ENDDAY'])
    period['interval'] = dt.timedelta(1)
    outfiletemplate = '{}/{}'.format(vicconfig['RESULT_DIR'],
                                     'vic_{:04d}_{:02d}_{:02d}.nc')
    ncformat = vicasc2ncconfig['ncformat']
    ncvarinfo = vicasc2ncconfig['ncvarinfo']
    print 'Converting output: {}'.format(outfiletemplate)
    vicasc2nc(filelist, domainfile, period, outfiletemplate,
              ncvarinfo, ncformat)

def readflow(ncfile, ncvar, start=None, end=None):
    nc = netCDF4.Dataset(ncfile, 'r')
    jd = [dt.datetime(x.year, x.month, x.day, x.hour, x.minute, x.second)
          for x in netCDF4.num2date(nc.variables['time'][:],
                                    nc.variables['time'].units,
                                    nc.variables['time'].calendar)]
    print '{}: From {} to {}'.format(ncfile, jd[0], jd[-1])
    df = pd.Series(nc.variables[ncvar][:,0], index=jd)
    nc.close()
    if start:
        df = df[start:]
    if end:
        df = df[:end]
    print df.describe()
    return df

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
        else:
            print 'Running: {}'.format(self.exe)
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
        print "Going to run {}".format(fargs)
        returnval = subprocess.call(fargs, stdout=loghandle, stderr=loghandle)
        if returnval != 0:
            raise SimulationError('return value not zero')
        if type(self.logfile) is str:
            loghandle.close()

class RouteSimulation(Simulation):

    def writeconfig(self, config):
        # Parse the configuration file
        configparser = ConfigParser.SafeConfigParser(allow_no_value=True)
        configparser.optionxform = str # preserve case of configuration keys
        configurationfile = config['GLOBALTEMPLATE']
        print 'Reading {}'.format(configurationfile)
        configparser.read(configurationfile)
        configtemplate = collections.OrderedDict() # preserve order of entries
        for section in configparser.sections():
            configtemplate[section] = collections.OrderedDict()
            for key, value in configparser.items(section):
                configtemplate[section][key] = value
        for key, value in config.items():
            fields = key.split('|')
            if len(fields) == 3:
                configtemplate[fields[1]][fields[2]] = value
        print 'Writing {}'.format(self.configfile)
        with open(self.configfile, 'w') as f:
            for section in configtemplate.keys():
                f.write('#-- ====================================== --#\n')
                f.write('[{}]\n'.format(section))
                for key, value in configtemplate[section].items():
                    f.write('{:20s}: {}\n'.format(key, value))

class VicSimulation(Simulation):

    def writeconfig(self, config):
        print 'Reading {}'.format(config['GLOBALTEMPLATE'])
        print 'Writing {}'.format(self.configfile)
        with open(config['GLOBALTEMPLATE'], 'r') as f, open(self.configfile, 'w') as g:
            for line in f:
                for key, value in config.items():
                    if re.match(key, line):
                       line = '{:20s} {}\n'.format(key, value)
                g.write(line)
        # Create the output directory if it does not exist
        make_sure_path_exists(config['RESULT_DIR'])

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
        self.config = collections.OrderedDict() # preserve order of entries
        for section in configparser.sections():
            self.config[section] = collections.OrderedDict()
            for key, value in configparser.items(section):
                self.config[section][key] = value
        if not self.config['COMPUTING']['poolsize']:
            self.poolsize = multiprocessing.cpu_count()
        else:
            self.poolsize = configparser.getint('COMPUTING','poolsize')

        # parse the nc variables info for vicasc2nc
        ncvarinfo = {}
        varlist = filterPick(self.config['VICASC2NC'].keys(), re.compile('^ncvar.*_name$'))
        varlist = [x.split('_name')[0] for x in varlist]
        for var in varlist:
            ncvarinfo[var] = {}
            ncvarinfo[var]['name'] = self.config['VICASC2NC'][var + '_name']
            ncvarinfo[var]['longname'] = self.config['VICASC2NC'][var + '_longname']
            ncvarinfo[var]['column'] = self.config['VICASC2NC'][var + '_column']
            ncvarinfo[var]['units'] = self.config['VICASC2NC'][var + '_units']
            ncvarinfo[var]['divideby'] = float(self.config['VICASC2NC'][var + '_divideby'])
        self.config['VICASC2NC']['ncvarinfo'] = ncvarinfo

        # initialize simulation time
        self.start = parse_date(self.config['SIMULATION']['start'])
        self.now = self.start
        self.end = parse_date(self.config['SIMULATION']['end'])
        self.runlength = dt.timedelta(configparser.getint('SIMULATION',
                                                          'runlength'))
        self.backstep = dt.timedelta(configparser.getint('SIMULATION',
                                                         'backstep'))

        # Initialize particles
        self.np = configparser.getint('FILTER', 'np')
        self.particles = [VicParticle(i+1, 1./self.np) for i in range(self.np)]
        for i in range(self.np):
            p = self.particles[i]
            p.vic['exe'] = self.config['MODELS']['vicexe']
            p.vic['infile'] = self.config['MODELS']['vicinfiletemplate'].format(i+1)
            p.vic['logfile'] = self.config['MODELS']['viclogfiletemplate'].format(i+1)
            p.route['exe'] = self.config['MODELS']['routeexe']
            p.route['infile'] = self.config['MODELS']['routeinfiletemplate'].format(i+1)
            p.route['logfile'] = self.config['MODELS']['routelogfiletemplate'].format(i+1)
            p.refvar = self.config['FILTER']['referencevar']
        self.normalize_weights()

    def __repr__(self):
        bstr = '{}\n'.format(self.config)
        bstr += '{}: {}\n'.format(self.poolsize, self.pool)
        bstr += '{}: {}\n'.format(self.np, [self.particles[i] for i in range(self.np)])
        bstr += '{} -- {}({},{}) -- {}\n'.format(self.start, self.now,
                                                 self.runlength, self.backstep,
                                                 self.end)
        return bstr

    def advance(self):
        self.now = self.runstate

    def run(self):
        # run all the particles
        # determine the time settings
        self.runstart = self.now
        self.runend = self.runstart + self.runlength - dt.timedelta(1)
        self.runstate = self.runend - self.backstep

        # write the configuration setting
        self.config['VICCONFIG']['STARTYEAR'] = self.runstart.year
        self.config['VICCONFIG']['STARTMONTH'] = self.runstart.month
        self.config['VICCONFIG']['STARTDAY'] = self.runstart.day
        self.config['ROUTECONFIG']['SECTION|OPTIONS|RUN_STARTDATE'] = \
            '{:04d}-{:02d}-{:02d}-00'.format(self.runstart.year, self.runstart.month,
                                             self.runstart.day)
        self.config['VICCONFIG']['ENDYEAR'] = self.runend.year
        self.config['VICCONFIG']['ENDMONTH'] = self.runend.month
        self.config['VICCONFIG']['ENDDAY'] = self.runend.day
        self.config['VICCONFIG']['STATEYEAR'] = self.runstate.year
        self.config['VICCONFIG']['STATEMONTH'] = self.runstate.month
        self.config['VICCONFIG']['STATEDAY'] = self.runstate.day
        self.config['ROUTECONFIG']['SECTION|OPTIONS|STOP_N'] = \
            self.runlength.days
        self.config['ROUTECONFIG']['SECTION|OPTIONS|REST_DATE'] = \
            '{:04d}-{:02d}-{:02d}-00'.format(self.runstate.year, self.runstate.month,
                                             self.runstate.day)
        self.config['ROUTECONFIG']['SECTION|HISTORY|RVICHIST_MFILT'] = \
            '{}'.format(self.runlength.days)

        # changes made to particles.self as part of particle.run() are not
        # persistent, because the run() part happens in different processes
        # that die when finished. Split into a run setup part that is
        # performed as part of the main script and let the run part just be
        # non-persistent changes. Needs refactoring and cleanup! It would be
        # nice if writing of the config files were still to happen in parallel
        for p in self.particles:
            p.configrun(self.config['VICCONFIG'], self.config['ROUTECONFIG'],
                        self.config['VICASC2ROUTENC'], self.config['VICASC2NC'])
            print '{:03d}: {}'.format(p.ID, p.resultfile)

        self.pool = multiprocessing.Pool(self.poolsize)
        results = [self.pool.apply_async(run, [self.particles[i]])
                   for i in range(self.np)]
        self.pool.close()
        self.pool.join()
        # Report errors
        for r in results:
            if r.get() is not None:
                raise(SimulationError(r.get()))

        print 'Evaluating particles'
        self.evaluate()
        self.normalize_weights()
        self.resample()
        print '{Particle: Weight ==> resample'
        for p in self.particles:
            print '{:03d}: {} ==> {:03d}'.format(p.ID, p.weight, p.sample)


    def evaluate(self):
        print 'Reading reference file: {}'.format(self.config['FILTER']['referencefile'])
        dfref = readflow(self.config['FILTER']['referencefile'],
                         self.config['FILTER']['referencevar'],
                         self.runstart, self.runend + dt.timedelta(1))
        for p in self.particles:
            p.evaluate(dfref)

    def normalize_weights(self):
        total = np.array([self.particles[i].weight for i in range(self.np)]).sum()
        for p in self.particles:
            p.weight /= total

    def resample(self):
        cdf = np.array([self.particles[i].weight for i in range(self.np)]).cumsum()
        unif = np.random.uniform(size=self.np)
        samples = [np.where(cdf > x)[0][0] for x in unif]
        for p, s in zip(self.particles, samples):
            p.sample = s

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

    def configrun(self):
        """Configure the run for a single particle"""
        pass

    def run(self):
        """Run a particle"""
        pass

    def evaluate(self):
        """Evaluate a particle (calculate its weight)"""
        pass

class VicParticle(Particle):
    """Particle for running VIC simulations (including routing)"""
    """Refactor: config info should be part of the particle class info,
       rather than passed around as part of the function arguments"""

    def __init__(self, ID, weight, state=None):
        # invoke the parent's initializer
        super(VicParticle, self).__init__(ID, weight, state)
        self.vic = {'exe': None, 'infile': None, 'logfile': None}
        self.route = {'exe': None, 'infile': None, 'logfile': None}

    def configrun(self, vicconfig, routeconfig, vicasc2routencconfig, vicasc2ncconfig):
        # Personalize the config info
        self.vicconfig = copy.copy(vicconfig)
        self.routeconfig = copy.copy(routeconfig)
        self.vicasc2routencconfig = vicasc2routencconfig
        self.vicasc2ncconfig = vicasc2ncconfig
        self.vicconfig['RESULT_DIR'] = vicconfig['RESULT_DIR'].format(self.ID)
        self.vicconfig['STATENAME'] = vicconfig['STATENAME'].format(self.ID)
        self.vicconfig['FORCING1'] = vicconfig['FORCING1'].format(self.ID)
        self.routeconfig['SECTION|OPTIONS|CASEID'] = \
            routeconfig['SECTION|OPTIONS|CASEID'].format(self.ID)
        self.routeconfig['SECTION|OPTIONS|CASE_DIR'] = \
            routeconfig['SECTION|OPTIONS|CASE_DIR'].format(self.ID)
        self.routeconfig['SECTION|INPUT_FORCINGS|DATL_PATH'] = \
            routeconfig['SECTION|INPUT_FORCINGS|DATL_PATH'].format(self.ID)
        print 'X-{:03d}: {}'.format(self.ID, self.routeconfig['SECTION|OPTIONS|CASE_DIR'])
        self.resultfile = self.routeconfig['SECTION|OPTIONS|CASE_DIR'] + '/hist/gunnison_{:03d}.rvic.h0a.{:04d}-{:02d}-{:02d}.nc'.format(self.ID, vicconfig['ENDYEAR'], vicconfig['ENDMONTH'], vicconfig['ENDDAY'])
        print '{:03d}: Setting result file {}'.format(self.ID, self.resultfile)
        self.vicsim = VicSimulation(self.vic['exe'], self.vic['infile'],
                                    ['-g', self.vic['infile']],
                                    self.vic['logfile'])
        self.routesim = RouteSimulation(self.route['exe'],
                                        self.route['infile'],
                                        [self.route['infile']],
                                        self.route['logfile'])

    def evaluate(self, dfref):
        print '{:03d}: refvar = {}'.format(self.ID, self.refvar)
        print '{:03d}: Reading result file: {}'.format(self.ID, self.resultfile)
        df = readflow(self.resultfile, self.refvar)
        print '{:03d}: {} {}'.format(self.ID, len(dfref), len(df))
        print '{:03d}: {} {}'.format(self.ID, dfref[0], df[0])
        print '{:03d}: {} {}'.format(self.ID, dfref[-1], df[-1])
        self.objective = np.sqrt(np.mean((df-dfref)**2))
        self.weight = 1./(self.objective+0.001)
        print '{:03d}: {} {}'.format(self.ID, self.objective, self.weight)

    def run(self):
        """Run a particle"""
        # Run VIC
        self.vicsim.writeconfig(self.vicconfig)
        self.vicsim.run()

        # Convert the VIC flux files into netCDF files that can be read by the
        # routing program
        print 'Completed VIC simulation starting conversion to netcdf'
        convertvicasc2route(self.vicconfig, self.routeconfig, self.vicasc2routencconfig)
        convertvicasc2nc(self.vicconfig, self.routeconfig['SECTION|DOMAIN|FILE_NAME'], self.vicasc2ncconfig)
        print 'Completed conversion to netcdf removing VIC flux files'

        # remove VIC files now that they have been converted
        inpath = self.vicconfig['RESULT_DIR']
        if inpath[-1] != '/':
            inpath += '/'
        filelist = glob.glob(inpath + self.vicasc2ncconfig['vicfilefilter'] + "*")
        for f in filelist:
            os.remove(f)

        print 'Going to run the routing part: {}'.format(self.route)
        # Route the flows
        print 'Logfile: {}'.format(self.route['logfile'])
        self.routesim.writeconfig(self.routeconfig)
        self.routesim.run()

if __name__ == '__main__':
    pf = ParticleFilter('/Users/nijssen/Dropbox/vicpf/config/gunnison/particlefilter/gunnison_pf.cfg')
    print pf
    #sys.exit()
    print "Advancing filter 1"
    pf.run()
#    print pf
#    print "Advancing filter 2"
#    pf.run()
#    print pf



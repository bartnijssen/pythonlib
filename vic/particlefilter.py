# This code is to support a particle filter implementation using VIC
from __future__ import division
import collections
import ConfigParser
import copy
import datetime as dt
import glob
import logging
import multiprocessing
import netCDF4
import numpy as np
import os
import pandas as pd
import re
import shutil
import subprocess
import sys
from utils.general import flatten, make_sure_path_exists
from utils.parse import parse_date
from vic.vicasc2routenc import vicasc2routenc
from vic.vicasc2nc import vicasc2nc, filterPick

loglevel_default = 'info'

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
    logging.debug('Converting output: {}'.format(outfiletemplate))
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
    ncfile = vicconfig['vicncfile']
    ncformat = vicasc2ncconfig['ncformat']
    ncvarinfo = vicasc2ncconfig['ncvarinfo']
    vicasc2nc(filelist, domainfile, period, ncfile,
              ncvarinfo, ncformat)

def readflow(ncfile, ncvar, start=None, end=None):
    nc = netCDF4.Dataset(ncfile, 'r')
    jd = [dt.datetime(x.year, x.month, x.day, x.hour, x.minute, x.second)
          for x in netCDF4.num2date(nc.variables['time'][:],
                                    nc.variables['time'].units,
                                    nc.variables['time'].calendar)]
    df = pd.Series(nc.variables[ncvar][:,0], index=jd)
    nc.close()
    if start:
        df = df[start:]
    if end:
        df = df[:end]
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
            logging.debug('Running: {}'.format(self.exe))
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
        logging.debug('Going to run {}'.format(fargs))
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
        logging.debug('Reading {}'.format(configurationfile))
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
        logging.debug('Writing {}'.format(self.configfile))
        with open(self.configfile, 'w') as f:
            for section in configtemplate.keys():
                f.write('#-- ====================================== --#\n')
                f.write('[{}]\n'.format(section))
                for key, value in configtemplate[section].items():
                    f.write('{:20s}: {}\n'.format(key, value))

class VicSimulation(Simulation):

    def writeconfig(self, config):
        logging.debug('Reading {}'.format(config['GLOBALTEMPLATE']))
        logging.debug('Writing {}'.format(self.configfile))
        with open(config['GLOBALTEMPLATE'], 'r') as f, open(self.configfile, 'w') as g:
            for line in f:
                for key, value in config.items():
                    if re.match(key, line):
                       line = '{:20s} {}\n'.format(key, value)
                g.write(line)

class ParticleFilter(object):
    """Particle filter"""
    def __init__(self, configparser, configurationfile):
        """Initialize a particle filter based on a configurationfile"""
        self.pool = None
        # Parse the configuration file
        self.configfile = configurationfile
        self.config = collections.OrderedDict() # preserve order of entries
        for section in configparser.sections():
            self.config[section] = collections.OrderedDict()
            for key, value in configparser.items(section):
                self.config[section][key] = value

        if not self.config['COMPUTING']['poolsize']:
            self.poolsize = multiprocessing.cpu_count()
        else:
            self.poolsize = configparser.getint('COMPUTING','poolsize')

        self.ID = self.config['FILTER']['runid']
        logging.info('Initializing particle filter: {}'.format(self.ID))

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

        # number of timestep to evaluate
        self.neval = configparser.getint('FILTER','neval')

        # Initialize particles
        self.np = configparser.getint('FILTER', 'np')
        self.particles = [VicParticle(i+1, 1./self.np) for i in range(self.np)]
        for i in range(self.np):
            p = self.particles[i]
            p.vic['exe'] = self.config['MODELS']['vicexe']
            p.route['exe'] = self.config['MODELS']['routeexe']
            p.refvar = self.config['FILTER']['referencevar']
            p.initstatefile = {}
            p.initstatefile['vic'] = self.config['VICCONFIG']['INIT_STATE'].format(p.ID)
            if self.config['ROUTECONFIG']['SECTION|OPTIONS|RUN_TYPE'].lower() == 'startup':
                p.initstatefile['route'] = \
                    self.config['ROUTECONFIG']['SECTION|INITIAL_STATE|FILE_NAME'].\
                    format(p.ID)
            p.neval = self.neval
        self.normalize_weights()
        logging.info('Initialized {} particles'.format(self.np))

        # set up run and archive paths for particles
        runpath = self.config['FILTER']['rundir']
        archivepath = self.config['FILTER']['archivedir']
        subdirs = ['log', 'history', 'settings']
        self.paths = {}
        self.paths['run'] = {}
        self.paths['archive'] = {}
        for subdir in subdirs:
            self.paths['run'][subdir] = '{}/{}'.format(runpath, subdir)
            self.paths['archive'][subdir] = '{}/{}'.format(archivepath, subdir)
        logging.debug('Creating paths for filter {}:'.format(self.ID))
        for topdir, dirdict in self.paths.iteritems():
            logging.debug('\t{:8s}:'.format(topdir))
            for key, val in dirdict.iteritems():
                logging.debug('\t\t{:6s}: {}'.format(key, val))
                make_sure_path_exists(val)

        subdirs = ['config', 'data', 'log', 'state', 'history']
        for p in self.particles:
            p.paths = {}
            p.paths['run'] = {}
            p.paths['archive'] = {}
            runpath = '{}/pid_{:03d}'.format(self.config['FILTER']['rundir'], p.ID)
            archivepath = '{}/pid_{:03d}'.format(self.config['FILTER']['archivedir'], p.ID)
            for subdir in subdirs:
                p.paths['run'][subdir] = '{}/{}'.format(runpath, subdir)
                p.paths['archive'][subdir] = '{}/{}'.format(archivepath, subdir)
            logging.debug('Creating paths for particle {:03d}:'.format(p.ID))
            for topdir, dirdict in p.paths.iteritems():
                logging.debug('\t{:8s}:'.format(topdir))
                for key, val in dirdict.iteritems():
                    logging.debug('\t\t{:6s}: {}'.format(key, val))
                    make_sure_path_exists(val)

        # set history file for filter
        filename = '{}.history'.format(self.ID)
        self.hist = open(os.path.join(self.paths['run']['history'], filename), 'w')

        # set history file for each particle
        for p in self.particles:
            filename = '{}.{:03d}.history'.format(self.ID, p.ID)
            p.hist = open(os.path.join(p.paths['run']['history'], filename), 'w')

        logging.info('Initialization complete')

    def __repr__(self):
        bstr = '{}\n'.format(self.config)
        bstr += '{}: {}\n'.format(self.poolsize, self.pool)
        bstr += '{}: {}\n'.format(self.np, [self.particles[i] for i in range(self.np)])
        bstr += '{} -- {}({},{}) -- {}\n'.format(self.start, self.now,
                                                 self.runlength, self.backstep,
                                                 self.end)
        return bstr

    def advance(self):
        self.now = self.runstate + dt.timedelta(1)

    def run(self):
        # run all the particles
        # determine the time settings
        self.runstart = self.now
        self.runend = self.runstart + self.runlength - dt.timedelta(1)
        self.runstate = self.runend - self.backstep
        logging.info('{} - {} - {}'.format(self.runstart.strftime('%Y-%m-%d'),
                                           self.runstate.strftime('%Y-%m-%d'),
                                           self.runend.strftime('%Y-%m-%d')))

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
            '{:04d}-{:02d}-{:02d}-00'.format(self.runstate.year,
                                             self.runstate.month,
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
            p.configrun(self.config['VICCONFIG'],
                        self.config['ROUTECONFIG'],
                        self.config['VICASC2ROUTENC'],
                        self.config['VICASC2NC'])
        # sys.exit('configrun finished')

        self.pool = multiprocessing.Pool(self.poolsize)
        results = [self.pool.apply_async(run, [self.particles[i]])
                   for i in range(self.np)]
        self.pool.close()
        self.pool.join()
        # Report errors
        for r in results:
            if r.get() is not None:
                raise(SimulationError(r.get()))

        logging.debug('Evaluating particles')
        self.evaluate()
        self.normalize_weights()
        self.resample()
        logging.debug('{Particle: Weight ==> resample')
        for p in self.particles:
            logging.debug('{:03d}: {} ==> {:03d}'.format(p.ID, p.weight, p.sample.ID))
            p.writehistory(self.runstart, self.runstate, self.runend)
            p.updateinitstate(self.ID, self.runstart, self.runstate)

    def archiverun(self):
        logging.debug('Archiving {}: {}-{}'.\
                      format(self.ID, self.runstart.strftime('%Y-%m-%d'),
                             self.runend.strftime('%Y-%m-%d')))
        for p in self.particles:
            p.archiverun()

    def archivefinal(self):
        logging.debug('Final archiving {}: {}-{}'.\
                      format(self.ID, self.start.strftime('%Y-%m-%d'),
                             self.end.strftime('%Y-%m-%d')))
        for p in self.particles:
            p.archivefinal()

        # copy the configuration file and the main particlefilter.py script
        shutil.copy(self.configfile, self.paths['archive']['settings'])
        shutil.copy(os.path.abspath(__file__), self.paths['archive']['settings'])

        filedict = {}
        # move the history file
        filelist = glob.glob('{}/*'.format(self.paths['run']['history']))
        for src in filelist:
            dst = os.path.join(self.paths['archive']['history'],
                                   os.path.basename(src))
            filedict[src] = dst
        # move the log file
        filelist = glob.glob('{}/*'.format(self.paths['run']['log']))
        for src in filelist:
            dst = os.path.join(self.paths['archive']['log'],
                                   os.path.basename(src))
            filedict[src] = dst


        # move the files
        logging.debug('\t  Files to be moved:')
        for src, dst in filedict.iteritems():
            logging.debug('\t    {} ===> {}'.format(src, dst))
            shutil.move(src, dst)

        # close history file
        self.hist.close()

    def evaluate(self):
        reffile = self.config['FILTER']['referencefile'].format(self.runend.strftime('%Y-%m-%d'))
        logging.debug('Reading reference file: {}'.format(reffile))
        dfref = readflow(reffile,
                         self.config['FILTER']['referencevar'],
                         self.runstart, self.runend + dt.timedelta(1))
        dfref = dfref[-(self.neval):]
        # dfref = None
        for p in self.particles:
            p.evaluate(dfref)

    def finalize(self):
        self.archivefinal()
        shutil.rmtree(self.config['FILTER']['rundir'])

    def normalize_weights(self):
        total = np.array([self.particles[i].weight for i in range(self.np)]).sum()
        for p in self.particles:
            p.weight /= total

    def resample(self):
        weights = [self.particles[i].weight for i in range(self.np)]
        cdf = np.array(weights).cumsum()
        unif = np.random.uniform(size=self.np)
        self.hist.write('{} {} {}:\n'.format(self.runstart.strftime('%Y-%m-%d'),
                                             self.runstate.strftime('%Y-%m-%d'),
                                             self.runend.strftime('%Y-%m-%d')))
        # self.hist.write('weights ({}): {}\n'.format(len(weights), list(weights)))
        # self.hist.write('cdf     ({}): {}\n'.format(len(cdf), list(cdf)))
        # self.hist.write('uniform ({}): {}\n'.format(len(unif), list(unif)))
        # samples = [np.where(cdf > x)[0][0] for x in unif]
        # samples = np.arange(self.np)
        # np.random.shuffle(samples)
        samples = np.ones(self.np, dtype=int) * (self.np-1);
        self.hist.write('samples ({}): {}\n'.format(len(samples), list(samples)))
        self.hist.write('set     ({}): {}\n'.format(len(set(samples)), set(samples)))
        for p, s in zip(self.particles, samples):
            p.sample = self.particles[s]
            logging.debug('\t{:03d}==>{:03d}'.format(p.ID, p.sample.ID))
        # for p in self.particles:
        #     p.sample = p
        #     logging.debug('\t{:03d}==>{:03d}'.format(p.ID, p.sample.ID))

    def runfilter(self):
        while (self.now + self.runlength) <= (self.end + dt.timedelta(1)):
            self.run()
            self.archiverun()
            self.advance()

class Particle(object):
    """Single particle in a particular filter"""

    def __init__(self, ID, weight, state=None):
        """Initialize a particle"""
        self.ID = ID
        self.weight = weight
        self.state = state
        logging.debug('initialized particle {}'.format(self.ID))

    def __repr__(self):
        bstr = '{}({}): {}'.format(self.ID, self.weight, self.state)
        return bstr

    def archiverun(self):
        """Archive the run for a single particle"""
        pass

    def archivefinal(self):
        pass

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

    def archiverun(self):
        logging.debug('\tArchiving particle {:03d}'.format(self.ID))
        # first the moves
        filedict = {}
        # list config files
        filelist = glob.glob('{}/*'.format(self.paths['run']['config']))
        for src in filelist:
            dst = os.path.join(self.paths['archive']['config'],
                                   os.path.basename(src))
            filedict[src] = dst
        # list data files
        filelist = glob.glob('{}/*.nc'.format(self.paths['run']['data']))
        filelist += glob.glob('{}/hist/*.nc'.format(self.paths['run']['data']))
        for src in filelist:
            dst = os.path.join(self.paths['archive']['data'],
                                   os.path.basename(src))
            filedict[src] = dst
        # list state files
        filelist = glob.glob('{}/restarts/*.nc'.format(self.paths['run']['data']))
        filelist += glob.glob('{}/*'.format(self.paths['run']['state']))
        for src in filelist:
            dst = os.path.join(self.paths['archive']['state'],
                                   os.path.basename(src))
            filedict[src] = dst

        # move the files
        logging.debug('\t  Files to be moved:')
        for src, dst in filedict.iteritems():
            logging.debug('\t    {} ===> {}'.format(src, dst))
            shutil.move(src, dst)

        # then the file to be deleted
        # data logs from routing program
        filelist = glob.glob('{}/logs/*'.format(self.paths['run']['data']))
        # rpointer files from routing program
        filelist += glob.glob('{}/restarts/rpointer'.format(self.paths['run']['data']))
        # data logs from VIC and routing
        filelist += glob.glob('{}/*'.format(self.paths['run']['log']))

        # delete the files
        logging.debug('\t  Files to be removed:')
        for src in filelist:
            logging.debug('\t    {}'.format(src))
            os.remove(src)

    def archivefinal(self):
        logging.debug('\tFinal archiving particle {:03d}'.format(self.ID))
        # close the history file
        self.hist.close()

        filedict = {}
        # move the history file for each particle
        filelist = glob.glob('{}/*'.format(self.paths['run']['history']))
        for src in filelist:
            dst = os.path.join(self.paths['archive']['history'],
                                   os.path.basename(src))
            filedict[src] = dst

        # move the files
        logging.debug('\t  Files to be moved:')
        for src, dst in filedict.iteritems():
            logging.debug('\t    {} ===> {}'.format(src, dst))
            shutil.move(src, dst)

    def configrun(self, vicconfig, routeconfig, vicasc2routencconfig, vicasc2ncconfig):
        # Personalize the config info
        startstr = '{:04d}-{:02d}-{:02d}'.format(vicconfig['STARTYEAR'],
                                                 vicconfig['STARTMONTH'],
                                                 vicconfig['STARTDAY'])
        endstr = '{:04d}-{:02d}-{:02d}'.format(vicconfig['ENDYEAR'],
                                               vicconfig['ENDMONTH'],
                                               vicconfig['ENDDAY'])
        statestr = '{:04d}-{:02d}-{:02d}'.format(vicconfig['STATEYEAR'],
                                                 vicconfig['STATEMONTH'],
                                                 vicconfig['STATEDAY'])
        self.vicconfig = copy.copy(vicconfig)
        self.routeconfig = copy.copy(routeconfig)
        self.vicasc2routencconfig = vicasc2routencconfig
        self.vicasc2ncconfig = vicasc2ncconfig
        self.vicconfig['RESULT_DIR'] = self.paths['run']['data']
        self.vicconfig['STATENAME'] = '{}/state.{:03d}.{}.{}'.\
            format(self.paths['run']['state'], self.ID, startstr, statestr)
        self.vicconfig['FORCING1'] = vicconfig['FORCING1'].format(self.ID)
        self.vicconfig['vicncfile'] = '{}/vic.{:03d}.{}.{}.nc'.\
                                      format(self.paths['run']['data'],
                                             self.ID, startstr, endstr)
        self.vicconfig['INIT_STATE'] = self.initstatefile['vic']

        self.routeconfig['SECTION|OPTIONS|CASE_DIR'] = \
            self.paths['run']['data']
        self.routeconfig['SECTION|OPTIONS|CASEID'] = \
            routeconfig['SECTION|OPTIONS|CASEID'].format(self.ID)
        self.routeconfig['SECTION|OPTIONS|CASE_DIR'] = \
            self.paths['run']['data']
        self.routeconfig['SECTION|INPUT_FORCINGS|DATL_PATH'] = \
            self.paths['run']['data']
        self.routeconfig['SECTION|INPUT_FORCINGS|DATL_FILE'] = \
            'rvic.{:03d}.{}.{}.nc'.format(self.ID, startstr, endstr)
        if self.routeconfig['SECTION|OPTIONS|RUN_TYPE'].lower() == 'startup':
            self.routeconfig['SECTION|INITIAL_STATE|FILE_NAME'] = \
                self.initstatefile['route']

        self.route['resultfile'] = self.routeconfig['SECTION|OPTIONS|CASE_DIR'] + \
                                   '/hist/{}.rvic.h0a.{}.nc'.\
                                   format(self.routeconfig['SECTION|OPTIONS|CASEID'],
                                          endstr)

        self.vic['infile'] = '{}/vic.global.{:03d}.{}.{}'.\
            format(self.paths['run']['config'], self.ID, startstr, endstr)
        self.vic['logfile'] = '{}/vic.log.{:03d}.{}.{}'.\
            format(self.paths['run']['log'], self.ID, startstr, endstr)
        self.vicsim = VicSimulation(self.vic['exe'], self.vic['infile'],
                                    ['-g', self.vic['infile']],
                                    self.vic['logfile'])

        self.route['infile'] = '{}/route.global.{:03d}.{}.{}'.\
            format(self.paths['run']['config'], self.ID, startstr, endstr)
        self.route['logfile'] = '{}/route.log.{:03d}.{}.{}'.\
            format(self.paths['run']['log'], self.ID, startstr, endstr)
        self.routesim = RouteSimulation(self.route['exe'],
                                        self.route['infile'],
                                        [self.route['infile']],
                                        self.route['logfile'])

    def evaluate(self, dfref):
        logging.debug('{:03d}: Reading result file: {}'.\
                      format(self.ID, self.route['resultfile']))

        df = readflow(self.route['resultfile'], self.refvar)
        df = df[-(self.neval):]

        # relative root mean squared error (smaller is better)
        # rrmse = np.sqrt(np.mean((df-dfref)**2))/np.mean(dfref)
        # normalized number of sign crossings (smaller is better)
        # ncs = ((np.sign(np.diff(df.values)) !=
        #         np.sign(np.diff(dfref.values))).sum() /
        #        (len(df)-1))
        # self.objective = 0.75*rrmse + 0.25*ncs

        # root mean squared error (smaller is better)
        self.objective = np.sqrt(np.mean((df-dfref)**2))

        self.weight = 1./(self.objective+0.001)
        # self.objective = 1
        # self.weight = 1
        logging.debug('{:03d}: {} {}'.format(self.ID, self.objective, self.weight))

    def run(self):
        """Run a particle"""
        # Run VIC
        self.vicsim.writeconfig(self.vicconfig)
        self.vicsim.run()

        # Convert the VIC flux files into netCDF files that can be read by the
        # routing program
        logging.debug('Completed VIC simulation starting conversion to netcdf')
        convertvicasc2route(self.vicconfig, self.routeconfig, self.vicasc2routencconfig)
        convertvicasc2nc(self.vicconfig, self.routeconfig['SECTION|DOMAIN|FILE_NAME'], self.vicasc2ncconfig)
        logging.debug('Completed conversion to netcdf removing VIC flux files')

        # remove VIC files now that they have been converted
        inpath = self.vicconfig['RESULT_DIR']
        if inpath[-1] != '/':
            inpath += '/'
        filelist = glob.glob(inpath + self.vicasc2ncconfig['vicfilefilter'] + "*")
        for f in filelist:
            os.remove(f)

        logging.debug('Going to run the routing part: {}'.format(self.route))
        # Route the flows
        logging.debug('Logfile: {}'.format(self.route['logfile']))
        self.routesim.writeconfig(self.routeconfig)
        self.routesim.run()

    def updateinitstate(self, runid, runstart, runstate):
        self.initstatefile['vic'] = '{}/state.{:03d}.{}.{}_{}'.\
            format(self.sample.paths['archive']['state'],
                   self.sample.ID,
                   runstart.strftime('%Y-%m-%d'),
                   runstate.strftime('%Y-%m-%d'),
                   runstate.strftime('%Y%m%d'))
        routestatedate = runstate + dt.timedelta(1)
        self.initstatefile['route'] = '{}/{}_{:03d}.r.{}.nc'.\
            format(self.sample.paths['archive']['state'],
                   runid,
                   self.sample.ID,
                   routestatedate.strftime('%Y-%m-%d-%H-%M-%S'))
        logging.debug('{:03d}({:03d}): {:40s} {:40s}'.\
                      format(self.ID, self.sample.ID,
                             self.initstatefile['vic'],
                             self.initstatefile['route']))

    def writehistory(self, runstart, runstate, runend):
        hstr = '{:03d} {} {} {} {} {} {:03d}'.\
            format(self.ID,
                   runstart.strftime('%Y-%m-%d'),
                   runstate.strftime('%Y-%m-%d'),
                   runend.strftime('%Y-%m-%d'),
                   self.objective,
                   self.weight,
                   self.sample.ID)
        self.hist.write('{}\n'.format(hstr))

if __name__ == '__main__':
    # get configuration file from command-line
    try:
        configurationfile = sys.argv[1]
    except:
        sys.exit('Usage: {} <configuration file>'.format(sys.argv[0]))

    # parse configuration file
    configparser = ConfigParser.SafeConfigParser(allow_no_value=True)
    configparser.optionxform = str # preserve case of configuration keys
    configparser.read(configurationfile)

    # initiate logging
    logfile = configparser.get('FILTER', 'logfile')
    try:
        loglevel = configparser.get('FILTER', 'loglevel').upper()
    except:
        loglevel = loglevel_default.upper()
    loglevel = getattr(logging, loglevel)
    logging.basicConfig(filename=logfile, filemode='w', level=loglevel)
    logging.info('Parsed configuration file {}'.format(configurationfile))

    # run filter
    pf = ParticleFilter(configparser, configurationfile)
    pf.runfilter()
    pf.finalize()

    logging.info('Particle filter completed')




# This code is to support a particle filter implementation using VIC
#import datetime as dt
#import matplotlib as mpl
#import matplotlib.pyplot as plt
import multiprocessing
#import numpy as np
#import pandas as pd
import subprocess
#import sys
from utils.general import flatten
#import uuid

class Error(Exception):
    """Base class for exceptions in this module"""
    pass

class SimulationError(Error):
    """Exception raised for simulation errors"""
    pass

class Simulation:
    """Model simulation"""
    def __init__(self, exe=None, args=None, logfile=None):
        """Initiate a simulation instance"""
        self.exe = exe
        self.args = args
        self.logfile = logfile

    def __str__(self):
        bstr = 'exe: {}\n'.format(self.exe)
        bstr += 'arguments:\n'
        if self.args:
            bstr += ', '.join(self.args)
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

class Particle():
    """Single particle in a particular filter"""

    def __init__(self, ID=None, weight=None, state=None):
        """Initiate a particle"""
#        self.ID = uuid.uuid1()
        self.ID = ID
        self.weight = weight
        self.state = state

#    def __str__(self):
#        bstr = '{}({}): {}'.format(str(self.ID)[:8], self.weight, self.state)
#        return bstr

    def __repr__(self):
        bstr = '{}({}): {}'.format(self.ID, self.weight, self.state)
        return bstr

    def advance(self, exe, infile, logfile):
        """Advance a particle"""
        self.simulation = Simulation(exe, ['-g', infile], logfile)
        self.simulation.run()

    def evaluate(self):
        """Evaluate a particle (calculate its weight)"""
        pass

def advance(particle, exe, infile, logfile):
    particle.advance(exe, infile, logfile)

if __name__ == '__main__':

    pool = multiprocessing.Pool()

    exe = '/d1/nijssen/vicpf/bin/vicNl'
    infiletemplate = '/d1/nijssen/vicpf/config/gunnison/vic/gunnison.vic.global.test_{:03d}'
    logfiletemplate = '/d1/nijssen/junk/run_{:03d}'

    Np = 10
    particles = [Particle(i+1, 1./Np) for i in range(Np)]
    results = [pool.apply_async(advance, [particles[i], exe, infiletemplate.format(i+1), logfiletemplate.format(i+1)]) for i in range(Np)]
    pool.close()
    pool.join()

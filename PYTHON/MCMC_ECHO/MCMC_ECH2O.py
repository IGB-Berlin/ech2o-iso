#/usr/bin/env python
# -*- coding: utf-8 -*-
# *************************************************
#
# Markov Chain Monte Carlo calibration algorithm for ECH2O
#
# -------
# Routine: Main program
# -------
# Author: S. Kuppel
# Created on 11/2016
# -------------------------------------------------

from pcraster import *
import subprocess
import numpy as np
import os, time, sys, glob, copy
from optparse import OptionParser

# --- Subroutine(s)
import Tools

# ---------
#  OPTIONS
# ---------
parser = OptionParser()

# --- Definition ---
# Configuration file
parser.add_option("--file",dest="file",metavar="FILE",
                  help="Name of the file that defines the calibration to perform")
# Output directory
parser.add_option("--pathout",dest="pathout",metavar="PATHOUT",
                  help="Name of output directory (current by default)")
# Name of the ECH2O executable
parser.add_option("--exe",dest="exe",metavar="exe",help="Name of the ECH2O exec file")
# Name of the reference ECH2O config file
parser.add_option("--cfg",dest="cfg",metavar="cfg",help="Name of the ECH2O config file")
# Output subdirectory
parser.add_option("--outdir",dest="outdir",metavar="outdir",help="Output directory")
# Number of CPUs used
parser.add_option("--ncpu",dest="ncpu",metavar="ncpu",help="Number of CPUs used")
# Only generating parameters ?
parser.add_option("--tlimit",dest="tlimit",metavar="tlimit",help="Time limit of one ECH2O run (in seconds)")

(options, args) = parser.parse_args()

# set working directory
class Config:
    PATH_MAIN = os.getcwd()+'/'

print
print '*********************************************************************************'
print 'MONTE CARLO CALIBRATION WITH ECH2O: '
print #'- construction of the trajectories'
#print '- forward runs'
#print '- storage of outputs and info for posterior analysis : elementary effects, etc.'

#############################################################################################
# Initialization of the calibration
# ---------------------------------

# -- General config file (where parameters info is taken from)
if options.file != None:
    [file_py,ext]=os.path.splitext(options.file)
    if len(glob.glob(options.file)) ==0 : 
        sys.exit('# STOP. The file that define the assimilation characteristics does not exist : \n   ...'+file_py+'.py')
    file = options.file
    frep = os.path.dirname(file)
    if frep == '':  file = os.path.join(Config.PATH_MAIN,file)
    print 'The user provided definition file is: '+options.file

# -- Number of CPUs used in parallel
if options.ncpu == None:
    options.ncpu = 1
Config.ncpu = copy.copy(options.ncpu)

# -- Output path
if options.outdir == None:
    options.outdir = os.path.splitext(options.file)[0]
#print options.outdir

# -- Time wall for ECH2O execution
if options.tlimit == None:
    Config.tlimit = ''
else:
    Config.tlimit = options.tlimit
    Config.tcmd = 'ulimit -t '+str(int(options.tlimit)*int(Config.ncpu))+' ;'

# -- ECH2O executable
if options.exe != None:
    if len(glob.glob(options.exe))==0:
        sys.exit('The user provided EXEC file was not found: '+options.exe)
    exe_ech2o = options.exe
    print 'The user provided EXEC file is: '+options.exe

# -- Config file
if options.cfg != None:
    if len(glob.glob(options.cfg))==0:
        sys.exit('The user provided CFG file was not found: '+options.cfg)
    cfg_ech2o = options.cfg
    print 'The user provided CFG file is: '+options.cfg

# -- Execution command: time limit and OMP use
Config.cmde_ech2o = ' '.join([Config.tcmd,
                              'OMP_NUM_THREADS='+str(Config.ncpu),
                              '../'+exe_ech2o,'../'+cfg_ech2o])
#print Config.cmde_ech2o


# --- Defining the various PATHs
Config.PATH_OUT  = os.path.abspath(os.path.join(Config.PATH_MAIN,options.outdir))
Config.PATH_EXEC = os.path.abspath(os.path.join(Config.PATH_OUT,'tmp'))
Config.PATH_SPA = os.path.abspath(os.path.join(Config.PATH_OUT,'Spatial'))
Config.PATH_CLIM = os.path.abspath(os.path.join(Config.PATH_MAIN,'..','Climate_30m'))
# where the base maps are first copied from
Config.PATH_SPA_REF = os.path.abspath(os.path.join(Config.PATH_MAIN,'Spatial')) 

# -- Creation of output directory
if len(glob.glob(Config.PATH_OUT)) == 0: os.system('mkdir '+ Config.PATH_OUT)
# if len(glob.glob(Config.PATH_EXEC)) == 0: os.system('mkdir '+ Config.PATH_EXEC)

# -- Copy def file in the output
os.system('cp '+file+' '+Config.PATH_OUT)

# -- Some verbose
print
print 'Maps & parameters:', Config.PATH_SPA
print 'Climatic forcing: ', Config.PATH_CLIM
print 'Outputs:          ', Config.PATH_OUT
print

# -- Copy of cfg and executable's symbolic link
if len(glob.glob(os.path.join(Config.PATH_OUT,exe_ech2o))) != 0:
    os.system('rm '+os.path.join(Config.PATH_OUT,exe_ech2o))
os.symlink(os.path.join(Config.PATH_MAIN, exe_ech2o) , os.path.join(Config.PATH_OUT,exe_ech2o) )
#if len(glob.glob(os.path.join(PATH_EXEC,cfg_ech2o))) != 0:
#    os.system('rm '+os.path.join(PATH_EXEC,cfg_ech2o))
os.system('cp '+os.path.join(Config.PATH_MAIN, cfg_ech2o)+' '+os.path.join(Config.PATH_OUT,cfg_ech2o))

# -- Creation of inputs directory (which will reflect parameter's sampling)
if len(glob.glob(Config.PATH_SPA)) == 0: os.system('mkdir '+ Config.PATH_SPA)
# Copy of reference data
os.system('cp -pr '+Config.PATH_SPA_REF+'/* '+Config.PATH_SPA)
    
#############################################################################################
# Initialization
# ------------------------------

#-- Import classes and setup from the def file
exec('from ' + file_py + ' import *')

# Initialize parameters and optimization variables
Tools.InitOpti(Opti, Paras, Site)
# Shuffle the first guess (to be implemented) ?
#if mstart==1:
#    Tools.RandomFG(Opti)

# Initialized maps and other run inputs
Tools.InitInputs(Opti, Paras, Site, Config)

# Get the observations used for optimisation
Data.names = Data.obs.keys()
Data.names.sort()

# Some output print
print '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
print ' Calibrated parameters ('+str(Opti.nvar)+'):'
print ' - '.join(Opti.names)
print
print ' Constraining data fluxes:'
for oname in Data.names:
    print '- '+oname
print
print

#############################################################################################
# Calibration loop
# -------------------

# -- First guess (= prior if mstart=0)
print '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
print 'Prior simulation'
# Create the inputs for ECH2O
Tools.CreateInputs(Opti, Paras, Site, Config, 0)
# Run ECHO
if len(glob.glob(Config.PATH_EXEC))==0:
    os.system('mkdir '+Config.PATH_EXEC)
Tools.runECH2O(Config)
if Tools.runOK(Data, Config) == 0:
    sys.exit('Could not run the prior ! Will stop here then.')
# Get the initial cost function value
Tools.CostFunction(Opti, Data, Config, 0)
# Output initial quantities
#Tools.ParasFG(Opti, Config)
Tools.ParasOutputs(Opti, Config, 0)
Tools.ObsOutputs(Data, Opti, Config, 0)
Tools.OptiOutputs(Opti, Data, Config, 0)
print

# -- Metropolis random-walk
print '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
print ' Entering the metropolis random-walk loop'
print
for it in xrange(1,Opti.nit+1):
    
    itout = '%03i' % int(it)
    print '- - - - - - - - - - - - - - - - - -'
    print 'Iteration',it, 'of',Opti.nit

    # Random walk
    Tools.ParasWalk(Opti, it)
    # Create the inputs for ECH2O
    Tools.CreateInputs(Opti, Paras, Site, Config, it)
    # Run ECH2O
    Tools.runECH2O(Config)
    # Check if it ran properly
    if Tools.runOK(Data, Config)==1:
        # Get the cost function
        Tools.CostFunction(Opti, Data, Config, it)
        # Check if it is accepted
        Tools.CriterionWalk(Opti, it)
        # Write parameters values for this sequence
        Tools.ParasOutputs(Opti, Config, it)
        Tools.ObsOutputs(Data, Opti, Config, it)
        Tools.OptiOutputs(Opti, Data, Config, it)
        # Clean
        os.system('rm -f '+Config.PATH_EXEC+'/*')
    else:
        # Sotrage the run log to check what's wrong
        os.system('mv '+Config.PATH_EXEC+'/ech2o.log '+Config.PATH_OUT+'/ech2o_'+itout+'.log')


# -- After optimisation ends
# Output final parameter set
#Tools.ParasFinal(Opti, Config)
#Tools.ObsFinal(Data, Opti, Config)

        
###############################################################################################
###############################################################################################

#/usr/bin/env python
# -*- coding: utf-8 -*-
# *************************************************
#
# Monte Carlo calibration algorithm for ECH2O
#
# -------
# Routine: Main program
# -------
# Author: S. Kuppel
# Created on 10/2016
# -------------------------------------------------

from pcraster import *
from pcraster.framework import *
from subprocess import *
import numpy as np
import os, time, sys, glob, copy
from optparse import OptionParser
from itertools import chain

# --- Subroutine(s)
import MC_Tools

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
    if len(glob.glob(options.file)) ==0 : sys.exit('# STOP. The file that define the assimilation characteristics does not exist : \n   ...'+file_py+'.py')
    file = options.file
    frep = os.path.dirname(file)
    if frep == '':  file = os.path.join(Config.PATH_MAIN,file)
    print 'The user provided definition file is: '+options.file

# -- Output path
if options.outdir == None:
    options.outdir = os.getcwd()+'/'+os.path.splitext(options.file)[0]

# -- Output path
if options.ncpu == None:
    options.ncpu = 1
Config.ncpu = copy.copy(options.ncpu)

# -- Execution command
if options.exe != None:
    if len(glob.glob(options.exe))==0:
        sys.exit('The user provided EXEC file was not found: '+options.exe)
    exe_ech2o = options.exe
    print 'The user provided EXEC file is: '+options.exe

if options.cfg != None:
    if len(glob.glob(options.cfg))==0:
        sys.exit('The user provided CFG file was not found: '+options.cfg)
    cfg_ech2o = options.cfg
    print 'The user provided CFG file is: '+options.cfg

cmde_ech2o = ' '.join(['../'+exe_ech2o,'../'+cfg_ech2o])

# --- Defining the various PATHs
Config.PATH_OUT  = os.path.abspath(os.path.join(Config.PATH_MAIN,options.outdir))
#Config.PATH_EXEC = os.path.abspath(os.path.join(Config.PATH_OUT,'tmp'))
Config.PATH_SPA = os.path.abspath(os.path.join(Config.PATH_OUT,'Spatial'))
Config.PATH_CLIM = os.path.abspath(os.path.join(Config.PATH_MAIN,'..','Climate'))
# where the base maps are first copied from
Config.PATH_SPA_REF = os.path.abspath(os.path.join(Config.PATH_MAIN,'Spatial')) 

# -- Creation of output directory
if len(glob.glob(Config.PATH_OUT)) == 0: os.system('mkdir '+ Config.PATH_OUT)
#if len(glob.glob(Config.PATH_EXEC)) == 0: os.system('mkdir '+ Config.PATH_EXEC)

# -- Copy def file in the output
os.system('cp '+file+' '+Config.PATH_OUT)

# -- Some verbose
print
print 'Maps & parameters:', Config.PATH_SPA
print 'Climatic forcing: ', Config.PATH_CLIM
print 'Outputs:          ', Config.PATH_OUT
#print 'Runs:             ', Config.PATH_EXEC
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
os.system('cp -p '+Config.PATH_SPA_REF+'/*.map '+Config.PATH_SPA)

#############################################################################################
# Initialization
# ------------------------------

#-- Import classes and setup from the def file
exec('from ' + file_py + ' import *')

#-- Parameters: from definition to all values
Paras.names = Paras.ref.keys()
Paras.n = len(Paras.names)
# Remove the default map files of calibrated param in the inputs directory
# --> helps checking early on if there is an improper map update
for pname in Paras.names:
    if Paras.ref[pname]['veg']==0:
        os.system('rm -f '+Config.PATH_SPA+'/'+Paras.ref[pname]['file']+'.map')
# Read dictionary to get all param setup
Opti.min = []
Opti.max = []
Opti.log = []
Opti.names = []
Opti.ind = []
Paras.ind = {}
ipar=0
ipar2=0
Paras.isveg = 0
for par in Paras.names:
    # Dimensions (soil or veg or 1)
    nr = Paras.ref[par]['soil']*(Site.ns-1)+Paras.ref[par]['veg']*(Site.nv-1)+1
    # Build vectors used in the optimisation
    Opti.min.append(np.repeat(Paras.ref[par]['min'],nr))
    Opti.max.append(np.repeat(Paras.ref[par]['max'],nr))
    Opti.log.append(np.repeat(Paras.ref[par]['log'],nr))
    # Link betwen params and all variables
    Opti.ind.append(np.repeat(ipar,nr))
    ipar+=1
    if nr>1: Paras.ind[par] = list(np.arange(ipar2,ipar2+nr,1))
    if nr==1: Paras.ind[par] = ipar2
    ipar2+=nr
    # For outputs
    if Paras.ref[par]['soil']==1:
        Opti.names.append([par + '_' + s for s in Site.soils])
    elif Paras.ref[par]['veg']==1:
        Opti.names.append([par + '_' + s for s in Site.vegs])
        Paras.isveg += 1
    else:
        Opti.names.append([par])
    
Opti.min = list(chain(*Opti.min))
Opti.max = list(chain(*Opti.max))
Opti.log = list(chain(*Opti.log))
Opti.names = list(chain(*Opti.names))
Opti.ind = list(chain(*Opti.ind))

# Total number of variables
Opti.nvar = len(Opti.min)

# Soils / units maps
setclone(Config.PATH_SPA+'/unit.map')

Site.bmaps = {}
for im in range(Site.ns):
    Site.bmaps[Site.soils[im]] = readmap(Config.PATH_SPA+'/'+Site.sfiles[im])
Site.bmaps['unit'] = readmap(Config.PATH_SPA+'/unit.map')
Site.bmaps['chanmask'] = readmap(Config.PATH_SPA+'/chanmask.map')
Site.bmaps['chanmask_NaN'] = readmap(Config.PATH_SPA+'/chanmask_NaN.map')
    
# Reference dictionary for vegetation inputs file
Opti.vref = {}
with open(Config.PATH_SPA_REF+'/'+Site.vfile,'r') as csvfile:
    paramread = list(csv.reader(csvfile, delimiter='\t'))
exit
# "Head": number of species and of params
Opti.vref['header'] = paramread[0][0:len(paramread[0])]
# All parameters values (keep strings!)
for iv in range(Site.nv):
    Opti.vref[iv] = paramread[iv+1][0:len(paramread[iv+1])]
# "Footers" : name of head1 and of parameters
Opti.vref['footer'] = paramread[Site.nv+1][0:len(paramread[Site.nv+1])]
Opti.vref['name'] = paramread[Site.nv+2][0:len(paramread[Site.nv+2])]

print 'Total number of variables :', Opti.nvar
print 'Number of iterations      :', Opti.nit
print
print '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
print ' Entering the optimisation loop...'
print

#############################################################################################
# Calibration loop
# -------------------
for it in range(Opti.nit):
    
    #exec("itout = '%0"+str(len(str(Opti.nit)))+"i' % int(it+1)")
    itout = '%03i' % int(it+1)
    Config.PATH_EXEC = Config.PATH_OUT+'/'+itout
    if len(glob.glob(Config.PATH_EXEC))==0:
	os.system('mkdir '+Config.PATH_EXEC)
    print 'Iteration',itout, 'of',Opti.nit

    # Parameters values
    MC_Tools.gen_paras(Opti, it+1)
    # Write parameters values
    MC_Tools.output_par(Opti, Config, it)
    # Create the inputs for ECH2O
    MC_Tools.create_inputs(Opti, Paras, Site, Config, it)
    # Run ECH2O
    #os.system('echo $OMP_NUM_THREADS')
    #os.system(cmde_ech2o)
    os.chdir(Config.PATH_EXEC)
    os.system('OMP_NUM_THREADS='+str(Config.ncpu)+' '+cmde_ech2o+' > '+Config.PATH_EXEC+'/ech2o.log')
    # Compute and save optimisation metrics (misfit etc.)
    # MC_Tools.manage_outputs(Data, Config, itout)

###############################################################################################
###############################################################################################

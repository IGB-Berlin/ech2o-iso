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
import subprocess
import numpy as np
import os, time, sys, glob, copy
from optparse import OptionParser
from itertools import chain

# --- Subroutine(s)
import ECH2O_Tools_SAS as ECH2O_Tools

# ---------
#  OPTIONS
# ---------
parser = OptionParser()

# ==============================================================================================
# --- Definition ---

# What do you want to do ?
# 'calib_sampling': generate ensemble of parameters sets for calibration
# 'calib_runs': runs the model for all sampled parameters
# 'forward_runs': runs the model for an ensemble of runs, usually the best configurations from the calibration. Allows to look at observations not used in calibration, e.g. maps

# 'sensi_morris': performs a Morris sensitivity analysis
parser.add_option("--mode",dest="mode",metavar="mode",
                  help="Switch ('calib_sampling','calib_runs','forwards_runs','sensi_morris')")

# Configuration file
parser.add_option("--file",dest="file",metavar="FILE",
                  help="Name of the file that defines the calibration to perform")
# Output directory
parser.add_option("--outdir",dest="outdir",metavar="outdir",help="Output directory")
# Number of CPUs used
parser.add_option("--ncpu",dest="ncpu",metavar="ncpu",help="Number of CPUs used")
# Use scratch ? (saves EcH2O tmp outputs on scratch, then saves on users after each run)
parser.add_option("--scratch",dest="scratch",metavar="scratch",help="Uses of localscratch (1: fastest storage) or shared (2: fast)")
# Resolution of EcH2O
parser.add_option("--Resol",dest="Resol",metavar="Resol",help="resolution (in m) for simulation, by default 30m)")
# Report BasinSummary.txt ? (by default = 0)
parser.add_option("--BSum",dest="BSum",metavar="BSum",help="report BasinSummary.txt files for each simulation")

# == Options for specific routines modes ==

# -- If mode == 'calib_sampling'
# Sampling method
parser.add_option("--sampling",dest="sampling",metavar="sampling",help="Sampling method (LHS, LHS_m = maximim, uniform)")

# -- If mode != 'calib_sampling' (i.e., EcH2O is actually run)
# Name of the ECH2O config tracking file (if mode!='calib_sampling', and tracking activated in simulation)
parser.add_option("--cfgTrck",dest="cfgTrck",metavar="cfgTrck",help="Name of the ECH2O configTrck file")
# Name of the ECH2O executable
parser.add_option("--exe",dest="exe",metavar="exe",help="Name of the ECH2O exec file")
# Name of the reference ECH2O config file
parser.add_option("--cfg",dest="cfg",metavar="cfg",help="Name of the ECH2O config file")
# Flux tracking activated ?
parser.add_option("--isTrck",dest="isTrck",metavar="isTrck",help="Switch of water tracking")
# Time limit for the runs
parser.add_option("--tlimit",dest="tlimit",metavar="tlimit",
                  help="Time limit of one ECH2O run (in seconds)")
# Trim outputs ? If so, the beginning and/or ending of outputs should not be saved (e.g. transient state)
parser.add_option("--trimB",dest="trimB",metavar="trimB",help="Drop the beginning of outputs: 0 if not, length otherwise")
parser.add_option("--trimL",dest="trimL",metavar="trimL",help="Length of the trim (if trimB>0): full trim by default, integer otherwise (has to be larger than total length - trimB)")
# Leap year ? (if spinup is multiple of ONE year of actual forcings) ?
parser.add_option("--leap",dest="leap",metavar="leap",help="Leap year (1) or not (0) for spinup multiple")

# -- If mode == 'calib_runs'
# Spinup ? (if post optim, mode = 2 ) ?
parser.add_option("--spinup",dest="spinup",metavar="spinup",help="Spinup switch (0) or length (days)")
# Restart ? if calibration stopped for some reasons
# --> 0 if not (by default), or 1: the code will take the antelast that worked
parser.add_option("--restart",dest="restart",metavar="restart",help="Restart (1) or not (0)")

# -- If mode == 'forward_runs'
# Input dir (if post optim, mode = 2) ?
parser.add_option("--inEns",dest="inEns",metavar="inEns",help="Ensenble input prefix")
# Size of the ensemble (if post optim, mode = 2 ) ?
parser.add_option("--nEns",dest="nEns",metavar="nEns",help="Size of the ensemble")
# Average maps ? (saves a looooooads of disk space)
parser.add_option("--MapAv",dest="MapAv",metavar="MapAv",
                  help="0 or 1, to temporally average spatial outputs)")
parser.add_option("--MapAvT",dest="MapAvT",metavar="MapAvT",
                  help="Averaging period (week, month, season)")

# If mode == 'sensi_morris'
# Only generating trajectories ?
parser.add_option("--MSinit",dest="MSinit",metavar="MSinit",
                  help="Switch for Morris sensitivity: 0=trajectories generation, 1=runs")
parser.add_option("--MSspace",dest="MSspace",metavar="MSspace",
                  help="Walk in the parameter space for Morris: 'trajectory' or 'radial'")

# Read the options
(options, args) = parser.parse_args()

# =====================================================================================================

# set working directory
class Config:
    PATH_MAIN = os.getcwd()+'/'
    cfgdir = PATH_MAIN+'Configs/'

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

# -- Number of CPUs used in parallel
if options.ncpu == None:
    options.ncpu = 1
Config.ncpu = copy.copy(options.ncpu)

# -- MC init: many things don't happen if it is the case
if options.mode == None:
    sys.exit("Choose which mode you're on mate!!")
else:
    Config.mode = copy.copy(options.mode)

# -- Restart: do not start from first iteration
if options.restart == None:
    Config.restart = 0
else:
    Config.restart = int(options.restart)
    if Config.restart > 1:
        sys.exit('Wrong value for restart')

# -- Output path
#print options.nEns
if options.outdir == None:
    if Config.mode != 'forward_runs':
        options.outdir = os.path.splitext(options.file)[0]
    else:
        tmp = options.cfg.split('_')[1].split('.')[0].split('-')
        if(len(tmp)==2):
            [pref1,pref2] = tmp
            options.outdir = 'Res.ens'+options.nEns+'_'+pref2+'.'+pref1+'.'+options.inEns
        if(len(tmp)==1):
            options.outdir = 'Res.ens'+options.nEns+'.'+tmp[0]+'.'+options.inEns
        if(len(tmp)>2):
            sys.exit('Error: incorrect config file format. Check the code my friend!')

# scratch?
if options.scratch != None:
    if int(options.scratch) == 1:
        Config.scratch = 1
    else:
        Config.scratch = 0
else:
    Config.scratch = 0

# -- Time wall for ECH2O execution
if options.tlimit == None:
    Config.tlimit = '100000'
else:
    Config.tlimit = options.tlimit
Config.tcmd = 'ulimit -t '+str(int(Config.tlimit)*int(Config.ncpu))+' ;'

# -- Spinup 
if options.spinup == None:
    Config.spinup = 0
else:
    Config.spinup = int(options.spinup)

# -- Maps averaging
if options.MapAv != None:
    Config.MapAv = int(options.MapAv)
    if Config.MapAv == 1:
        if options.MapAvT in ['week','month','season']:
            Config.MapAvT = options.MapAvT
        else:
            sys.exit('Wrong maps averaging option!')
else:
    Config.MapAv = 0

# -- Resolution (default 30m, used with reference Spatial folder is picked up)
if options.Resol != None:
    Config.Resol = copy.copy(options.Resol)
else:
    Config.Resol = '30m'

# -- Report BasinSummary.txt
if options.BSum != None:
    Config.repBS = int(options.BSum)
else:
    Config.repBS = 0

# -- MS init: many things don't happen if it is the case
if Config.mode == 'sensi_morris':
    if options.MSinit == None:
        Config.MSinit = 1
        sys.exit("Choose if you're initializating the MS sampling (or not)!")
    else:
        Config.MSinit = int(options.MSinit)
        
    if options.MSspace == None:
        Config.MSspace = 'trajectory'
        #sys.exit("Choose how you walk the Morris space (trajectory or radial)!")
    else:
        Config.MSspace = copy.copy(options.MSspace)

# -- Run ECH2O?
Config.runECH2O = 1
if Config.mode == 'calib_sampling' or (Config.mode=='sensi_morris' and Config.MSinit==1):
    Config.runECH2O = 0

#print options.outdir

print
print '*********************************************************************************'
if Config.mode.split('_')[0] == 'calib':
    print 'CALIBRATION with EcH2O: '
elif Config.mode == 'forward_runs':
    print 'ENSEMBLE RUNS with EcH2O'
elif Config.mode == 'sensi_morris':
    print 'MORRIS SENSITIVITY with EcH2O: '
    print '- construction of the trajectories'
    print '- forward runs'
    print '- storage of outputs and info for posterior analysis : elementary effects, etc.' 
    print '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
    print

# -- Calibration: all parameter path
if Config.mode.split('_')[0] == 'calib':
    Config.PATH_PAR  = os.path.abspath(os.path.join(Config.PATH_MAIN,'Parameters_samples'))
    Config.FILE_PAR = Config.PATH_PAR+'/'+options.outdir.split('.')[0]+'_parameters.'
    # -- Creation of output directory
    if len(glob.glob(Config.PATH_PAR)) == 0: os.system('mkdir '+ Config.PATH_PAR)
    # -- Some verbose
    print
    print "Parameter samples' directory:          ", Config.PATH_PAR
    print

# -- Sensitivity: all parameter path
if Config.mode == 'sensi_morris':

    print
    if(Config.MSspace=='trajectory'):
        Config.PATH_TRAJ  = os.path.abspath(os.path.join(Config.PATH_MAIN,'Trajectories'))
        print "Trajectories directory:          ", Config.PATH_TRAJ
    if(Config.MSspace=='radial'):
        Config.PATH_TRAJ  = os.path.abspath(os.path.join(Config.PATH_MAIN,'RadialPoints'))
        print "Radial points directory:         ", Config.PATH_TRAJ
    Config.FILE_TRAJ = Config.PATH_TRAJ+'/'+options.outdir.split('.')[0]
    # -- Creation of output directory
    if len(glob.glob(Config.PATH_TRAJ)) == 0: os.system('mkdir '+ Config.PATH_TRAJ)
    # -- Output of elementary effects
    if(Config.MSinit==0):
        Config.PATH_EE  = os.path.abspath(os.path.join(Config.PATH_MAIN,'ElementaryEffects'))
        print "Elementary effects directory:          ", Config.PATH_EE
        Config.FILE_EE = Config.PATH_EE+'/'+options.outdir.split('.')[0]
    print

# -- For runs
#print Config.runECH2O

if Config.runECH2O == 1:
    # -- Execution command
    if options.exe != None:
        if len(glob.glob(options.exe))==0:
            sys.exit('The user provided EXEC file was not found: '+options.exe)
        exe_ech2o = options.exe
        print 'The user provided EXEC file is: '+options.exe

    if options.cfg != None:
        cfg_ech2o = options.cfg+'.ini'
        if Config.mode == 1 and len(glob.glob(Config.cfgdir+cfg_ech2o))==0:
            sys.exit('The user provided CFG file was not found: '+cfg_ech2o)
        if Config.mode == 2 and len(glob.glob(Config.cfgdir+cfg_ech2o))==0:
            sys.exit('The user provided CFG file was not found: '+cfg_ech2o)
        print 'The user provided CFG file is: '+cfg_ech2o

        # Spinup?
        if Config.mode == 'forward_runs' and Config.spinup > 0:
            print '~~~~ Spinup activated ! ~~~~~~~'
            if options.leap == None:
                Config.leap = 0
            else:
                Config.leap = int(options.leap)
            cfgSpin_ech2o = options.cfg.split('.')[0].split('-')[0]+'-spin'+str(Config.spinup)+'Y.ini'
            if len(glob.glob(Config.cfgdir+cfgSpin_ech2o))==0:
                sys.exit('The user provided spinup CFG file was not found: '+cfgSpin_ech2o)
            else:
                print 'The user provided spinup CFG file is: '+cfgSpin_ech2o
        
    # We want to run with isotopes...
    if options.isTrck != None:
        Config.isTrck = int(options.isTrck)
    else:
        Config.isTrck = 0

    if Config.isTrck == 1:
        if options.cfgTrck != None:
            cfgTrck_ech2o = options.cfgTrck+'.ini'
        else:
            cfgTrck_ech2o = options.cfg.split('_')[0]+'Trck_'+options.cfg.split('_')[1]+'.ini'

        if len(glob.glob(Config.cfgdir+cfgTrck_ech2o))==0:
            sys.exit('The user provided CFGtrck file was not found: '+cfgTrck_ech2o)
        else:
            print 'The user provided CFGtrck file is: '+cfgTrck_ech2o

        # Spinup (depreciated)?
        if Config.mode == 'forward_runs' and Config.spinup > 0:
            [pref1,pref2] = options.cfg.split('-')[0].split('_')
            cfgTrckSpin_ech2o = pref1+'Trck_'+pref2+'-spin.ini'
            if len(glob.glob(Config.cfgdir+cfgTrckSpin_ech2o))==0:
                sys.exit('The user provided spinup CFGtrck file was not found: '+cfgTrckSpin_ech2o)
            else:
                print 'The user provided spinup CFGtrck file is: '+cfgTrckSpin_ech2o


    # -- Execution command: time limit and OMP use
    Config.cmde_ech2o = ' '.join([Config.tcmd,'OMP_NUM_THREADS='+str(Config.ncpu),
                                  './'+exe_ech2o,cfg_ech2o])
    # Spinup? (depreciated)
    if Config.mode == 'forward_runs' and Config.spinup > 0:
        Config.spin_ech2o = ' '.join([Config.tcmd,'OMP_NUM_THREADS='+str(Config.ncpu),
                                      './'+exe_ech2o,cfgSpin_ech2o])

    # -- Get the parallel job number, based on the output dir name
    Config.numsim = options.outdir.split('.')[::-1][0]

    # --- Defining the various PATHs
    Config.PATH_OUT = os.path.abspath(os.path.join(Config.PATH_MAIN,options.outdir))
    Config.PATH_SPA = os.path.abspath(os.path.join(Config.PATH_OUT,'Spatial'))
        
    # -- Define execution directoy, depends on scratch options
    if Config.scratch > 0:
        if Config.scratch == 1:
            Config.PATH_EXEC = '/scratch/users/s07as8/'+options.outdir
        if Config.scratch == 2:
            Config.PATH_EXEC = '/nobackup/users/s07as8/'+options.outdir
    else:
        Config.PATH_EXEC = os.path.abspath(os.path.join(Config.PATH_OUT,'tmp'))    

if Config.mode == 'forward_runs':

    if options.inEns != None:
        Config.nEns = int(options.nEns)
        Config.FILE_PAR = Config.PATH_MAIN+'Input_ensemble/'+options.inEns+'.'+options.nEns+'bestParams.txt'
        if len(glob.glob(Config.FILE_PAR))==0:
            sys.exit('The param file (based on ensemble set) was not found: '+Config.FILE_PAR)
        print
        print 'The ensemble param file is : '+Config.FILE_PAR
    else:
        sys.exit('The size of ensemble simulations needs to be specified (--nEns)')

if Config.runECH2O == 1:
    # where the base maps are first copied from
    Config.PATH_SPA_REF = os.path.abspath(os.path.join(Config.PATH_MAIN,'Spatial_'+Config.Resol+'m'))
    #Config.PATH_SPA_REF = os.path.abspath(os.path.join(Config.PATH_MAIN,'Spatial_'+Config.Resol+'_old'))
    #Config.PATH_SPA_REF = os.path.abspath(os.path.join(Config.PATH_MAIN,'Spatial_'+Config.Resol+'_noveg'))
    #Config.PATH_CLIM_REF = os.path.abspath(os.path.join(Config.PATH_MAIN,'../Climate_'+str(Config.Resol)+'m')) 
    Config.PATH_CLIM_REF = os.path.abspath(os.path.join(Config.PATH_MAIN,'../Climate'))
    Config.PATH_CLIM_REF_P = os.path.abspath(os.path.join(Config.PATH_MAIN,'../Climate/Steps'))    

    # -- Creation of output directory
    if len(glob.glob(Config.PATH_OUT)) == 0: 
        os.system('mkdir '+ Config.PATH_OUT)

    # -- Creation of inputs directory (spatial will reflect parameter's sampling)
    if len(glob.glob(Config.PATH_SPA)) == 0: 
        os.system('mkdir '+ Config.PATH_SPA)
    # Copy of reference data
    os.system('cp -p '+Config.PATH_SPA_REF+'/*.map '+Config.PATH_SPA)
    os.system('cp -p '+Config.PATH_SPA_REF+'/SpeciesParams*.tab '+Config.PATH_SPA)

    # Creating the execution directory
    if len(glob.glob(Config.PATH_EXEC))==0:
        os.system('mkdir '+Config.PATH_EXEC)
    else:
        os.system('rm -f '+Config.PATH_EXEC+'/*')

    # -- Some verbose
    print
    print 'Original/template maps & parameters:', Config.PATH_SPA_REF
    print 'Original climate data:', Config.PATH_CLIM_REF
    print
    if Config.scratch > 0:
        print '-----------------------------------------'
        print "Scratch storage activated! Hopefully that's gonna speed up things..."
        print 'Temporary maps & parameters:', Config.PATH_SPA
        #print 'Temporary climate data:', Config.PATH_CLIMS
        print 'Temporary outputs:', Config.PATH_EXEC
        print '-----------------------------------------'
    else:
        print 'Maps & parameters:', Config.PATH_SPA
    #print 'Climatic forcing: ', Config.PATH_CLIM
    print 'Final outputs:          ', Config.PATH_OUT

    print

    # -- Symbolic link to executable
    if len(glob.glob(os.path.join(Config.PATH_OUT,exe_ech2o))) != 0:
        os.system('rm '+os.path.join(Config.PATH_OUT,exe_ech2o))
    os.symlink(os.path.join(Config.PATH_MAIN, exe_ech2o) , os.path.join(Config.PATH_OUT,exe_ech2o))

    # ===
    # -- Copy config files (on users even if there's scratch, for debugging)
    os.system('cp '+os.path.join(Config.PATH_MAIN, Config.cfgdir, cfg_ech2o)+' '+os.path.join(Config.PATH_OUT, cfg_ech2o))
    if Config.isTrck == 1:
        os.system('cp -p '+os.path.join(Config.PATH_MAIN, Config.cfgdir, cfgTrck_ech2o)+' '+os.path.join(Config.PATH_OUT,cfgTrck_ech2o))
    # -- Modify template config file with config-specific info
    with open(os.path.join(Config.PATH_OUT, cfg_ech2o),'a') as fw:
        fw.write('\n\n\n#Simulation-specific folder section\n#\n\n')
        fw.write('Maps_Folder = '+Config.PATH_SPA+'\n')
        fw.write('Clim_Maps_Folder = '+Config.PATH_CLIM_REF+'\n')
        fw.write('Output_Folder = '+Config.PATH_EXEC+'\n')
        fw.write('ClimateZones = ClimZones_'+Config.Resol+'m.map\n')
        fw.write('Isohyet_map = isohyet_'+Config.Resol+'m.map\n')
        if Config.isTrck == 1:
            fw.write('Tracking = 1\n')
            fw.write('TrackingConfig = '+os.path.join(Config.PATH_OUT,cfgTrck_ech2o)+'\n')
        else:
            fw.write('Tracking = 0')
        fw.write('\n\n')

    if Config.mode == 'forward_runs' and Config.spinup > 1:
        os.system('cp -p '+os.path.join(Config.PATH_MAIN, Config.cfgdir, cfgSpin_ech2o)+' '+os.path.join(Config.PATH_OUT,cfgSpin_ech2o))
        if Config.isTrck == 1:
            os.system('cp -p '+os.path.join(Config.PATH_MAIN, Config.cfgdir, cfgTrckSpin_ech2o)+' '+os.path.join(Config.PATH_OUT,cfgTrckSpin_ech2o))

        with open(os.path.join(Config.PATH_OUT, cfgSpin_ech2o),'a') as fw:
            fw.write('\n\n\n#Simulation-specific folder section\n#\n\n')
            fw.write('Maps_Folder = '+Config.PATH_SPA+'\n')
            fw.write('Clim_Maps_Folder = '+Config.PATH_CLIM_REF+'\n')
            fw.write('Output_Folder = '+Config.PATH_EXEC+'\n')
            if Config.isTrck == 1:
                fw.write('Tracking = 1\n')
                fw.write('TrackingConfig = '+os.path.join(Config.PATH_OUT,cfgTrckSpin_ech2o)+'\n')
            else:
                fw.write('Tracking = 0')
            fw.write('\n\n')

    # -- Copy def file in the output
    os.system('cp '+file+' '+Config.PATH_OUT)

    
#############################################################################################
# Initialization
# ------------------------------

#-- Import classes and setup from the def file
sys.path.insert(0,Config.PATH_MAIN)
exec('from ' + file_py + ' import *')

# -- Trim: only saves the time steps within the trim
if options.trimB == None:
    Config.trimB = 1
else:
    Config.trimB = int(options.trimB)
if options.trimL == None:
    Config.trimL = Data.lsim - Config.trimB + 1
else:
    if int(options.trimL) <= Data.lsim - Config.trimB+1:
        Config.trimL = int(options.trimL)
    else:
        sys.exit('Error: the specified output slicing start+length goes beyond simulation time!')

# Bare rock simulation
if Opti.simRock == None:
    Opti.simRock = 0

# -- Used for simulations
Config.treal = [Data.simbeg+timedelta(days=x) for x in range(Config.trimL)]
#sys.exit()

# -- Parameters: from definition to all values
Paras.names = Paras.ref.keys()
Paras.names.sort()
Paras.n = len(Paras.names)
# Read dictionary to get all param setup
Opti.min = []
Opti.max = []
Opti.log = []
Opti.names = []
Opti.ind = []
Opti.comp = []
Paras.ind = {}
ipar=0
ipar2=0
Paras.isveg = 0
for par in Paras.names:
        
    # Dimensions (soil or veg or 1)
    nr = Paras.ref[par]['soil']*(Site.ns-1)+Paras.ref[par]['veg']*(Site.nv-1)+1

    # Build vectors used in the optimisation
    if Config.mode != 'forward_runs':
        if type(Paras.ref[par]['min'])==float or type(Paras.ref[par]['min'])==int:
            Opti.min = Opti.min+ [Paras.ref[par]['min']]
            Opti.max = Opti.max+ [Paras.ref[par]['max']]
        else:
            Opti.min = Opti.min+ Paras.ref[par]['min']
            Opti.max = Opti.max+ Paras.ref[par]['max']

        Opti.log = Opti.log+list(np.repeat(Paras.ref[par]['log'],nr))

    # Link betwen params and all variables
    Opti.ind = Opti.ind+list(np.repeat(ipar,nr))
    ipar+=1
    if nr>1: Paras.ind[par] = list(np.arange(ipar2,ipar2+nr,1))
    if nr==1: Paras.ind[par] = ipar2
    ipar2+=nr
    # For outputs
    if Paras.ref[par]['soil']==1:
        Opti.names = Opti.names + [par + '_' + s for s in Site.soils]
        Opti.comp = Opti.comp + [i for i in range(Site.ns)]
    elif Paras.ref[par]['veg']==1:
        Opti.names = Opti.names + [par + '_' + s for s in Site.vegs]
        Opti.comp = Opti.comp + [i for i in range(Site.nv)]
        Paras.isveg += 1
    else:
        Opti.names = Opti.names + [par]
        Opti.comp = Opti.comp + [0]
    
#Opti.min = list(chain(*Opti.min))
#Opti.max = list(chain(*Opti.max))
#Opti.log = list(chain(*Opti.log))
#Opti.names = list(chain(*Opti.names))
#Opti.ind = list(chain(*Opti.ind))

#print Opti.names
#print Opti.min
#print Opti.max
#print Opti.log
#print Opti.ind

# Total number of variables
Opti.nvar = len(Opti.names)

# -- Generate calibration parameters samples
if Config.mode == 'calib_sampling':

    if options.sampling == None:
        Config.sampling = 'uniform'
    else:
        Config.sampling = options.sampling

    # Total number of samples
    Opti.nsamptot = Opti.nit*int(Config.ncpu)
    print
    print str(Opti.nsamptot)+' sets of '+str(Opti.nvar)+' parameters will be generated...'

    # Generate enough parameters sets and save them
    ECH2O_Tools.gen_paras(Opti, Config)
    # Done
    print 'Parameters sample generation done.'

# -- Retrieve previously generated samples
if Config.mode == 'calib_runs':

    print 'Get parameters samples for this job...'
    # Read parameters sample for this parallel run
    ECH2O_Tools.get_par(Opti, Config)

# IN post-optim mode, read directly the params from the file
if Config.mode == 'forward_runs':

    print 'Get parameter set for these ensemble runs...'

    tmp = list(np.genfromtxt(Config.FILE_PAR,delimiter=',',dtype= '|S20',unpack=True)[0])

    #print Opti.names
    #print tmp
    if tmp != Opti.names:
        sys.exit("The definition file and input parameter file ain't matching!")

    #pnames = ','.join([str(tmp[id.x]) for idx in range(Opti.nvar)])
    #print tmp#pnames
    #print Opti.names
    Opti.xpar = np.genfromtxt(Config.FILE_PAR,delimiter=',',unpack=True)[1::]
    #Opti.xpar = np.genfromtxt(Config.FILE_PAR,delimiter=',',skip_header=1)
    #Opti.x = [tmp[int(Config.numsim)-1][x] for x in range(Opti.nvar)]

# -- Generate morris trajectories for sensitivity analysis
if Config.mode == 'sensi_morris':

    # Normalized step: plus-minus 0.5 of the normalized range of each parameter
    Opti.stepN = np.zeros((Opti.nvar),np.float64) + 0.5

    if Config.MSinit == 1: 
        ECH2O_Tools.morris_trajs(Config, Opti)
        # Done
        print 'Parameters trajectory generation done.'

    else:
        # Get the trajectory
        f_in = Config.PATH_TRAJ+'/'+options.outdir.split('.')[0]+'.Bstar_traj'+Config.numsim+'.txt'
        # print f_in
        Opti.xpar = np.genfromtxt(f_in,delimiter=',',skip_header=1)
        # print Opti.xpar.shape
        # Reconstruct step (+- 0.5)
        if(Config.MSspace=='trajectory'):
            Opti.dx = np.diff(Opti.xpar,axis=0)
        elif(Config.MSspace=='radial'):
            Opti.dx = Opti.xpar[1::,:]-Opti.xpar[0,:]
        Opti.dx[Opti.dx!=0] = Opti.dx[Opti.dx!=0] / np.abs(Opti.dx[Opti.dx!=0]) * 0.5
        if(np.ptp(Opti.dx) != 1.0 or np.min(Opti.dx)!= -0.5 or np.max(Opti.dx)!= 0.5):
            sys.exit('Error: THe fetched Bnorm has a problem...')

    # Total number of runs
    Opti.nruns = (Opti.nvar+1) * Opti.nr

# Calibration or SA runs: Check that the same parameters are used
if Config.mode == 'calib_runs' or (Config.mode == 'sensi_morris' and Config.MSinit==0):
    if(len(Opti.xpar[0])!=len(Opti.names)):
        sys.exit("The definition file and input parameter file ain't matching!")
        
# -- Whenever there's EcH2O runs: a few others initializations
if Config.runECH2O == 1:

    # Get the observations used for optimisation
    Data.names = Data.obs.keys()
    Data.names.sort()

    # Remove the default map files of calibrated param in the inputs directory
    # --> helps checking early on if there is an improper map update
    for pname in Paras.names:
        if Paras.ref[pname]['veg']==0:
            os.system('rm -f '+Config.PATH_SPA+'/'+Paras.ref[pname]['file']+'.map')

    # Soils / units maps
    Config.cloneMap = boolean(readmap(Config.PATH_SPA+'/base.map'))
    setclone(Config.PATH_SPA+'/base.map')

    Site.bmaps = {}
    for im in range(Site.ns):
        Site.bmaps[Site.soils[im]] = readmap(Config.PATH_SPA+'/'+Site.sfiles[im])
    Site.bmaps['unit'] = readmap(Config.PATH_SPA+'/unit.map')
    if Opti.simRock == 1:
        # Site.bmaps['nolowK'] = readmap(Config.PATH_SPA+'/unit.nolowK.map')    
        Site.bmaps['rock'] = readmap(Config.PATH_SPA+'/unit.rock.map')    
    
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

    # In ensemble runs, check if spinup: if not, copy reference maps to files
    # with right name to match the "generic" config file (spinup or not)
    #if Config.spinup == 0 :
#    os.system('cp '+Config.PATH_SPA+'/streamflow.map '+Config.PATH_SPA+'/Q.map') # Streamflow
#    os.system('rm -f '+Config.PATH_SPA+'/streamflow.map')
#    os.system('cp '+Config.PATH_SPA+'/swe.map '+Config.PATH_SPA+'/SWE.map')      # Snow water equiv
#    os.system('rm -f '+Config.PATH_SPA+'/swe.map')
#    os.system('cp '+Config.PATH_SPA+'/SWC1.map '+Config.PATH_SPA+'/SWC.L1.map')  # Soil moisture L1
#    os.system('rm -f '+Config.PATH_SPA+'/SWC1.map')
#    os.system('cp '+Config.PATH_SPA+'/SWC2.map '+Config.PATH_SPA+'/SWC.L2.map')  # Soil moisture L2
#    os.system('rm -f '+Config.PATH_SPA+'/SWC2.map')
#    os.system('cp '+Config.PATH_SPA+'/SWC3.map '+Config.PATH_SPA+'/SWC.L3.map')  # Soil moisture L3
#    os.system('rm -f '+Config.PATH_SPA+'/SWC3.map')
#    os.system('cp '+Config.PATH_SPA+'/soiltemp.map '+Config.PATH_SPA+'/Ts.map')  # Soil temp
#    os.system('rm -f '+Config.PATH_SPA+'/soiltemp.map')
    
#    if Config.isTrck == 1:
#        os.system('cp '+Config.PATH_SPA+'/dD_snowpack.map '+Config.PATH_SPA+'/dD.snowpack.map') 
#        os.system('cp '+Config.PATH_SPA+'/dD_surface.map '+Config.PATH_SPA+'/dD.surface.map')
#        os.system('cp '+Config.PATH_SPA+'/dD_soil1.map '+Config.PATH_SPA+'/dD.L1.map')
#        os.system('cp '+Config.PATH_SPA+'/dD_soil2.map '+Config.PATH_SPA+'/dD.L2.map')
#        os.system('cp '+Config.PATH_SPA+'/dD_soil3.map '+Config.PATH_SPA+'/dD.L3.map')
#        os.system('cp '+Config.PATH_SPA+'/dD_groundwater.map '+Config.PATH_SPA+'/dD.GW.map')

#        os.system('cp '+Config.PATH_SPA+'/d18O_snowpack.map '+Config.PATH_SPA+'/d18O.snowpack.map') 
#        os.system('cp '+Config.PATH_SPA+'/d18O_surface.map '+Config.PATH_SPA+'/d18O.surface.map')
#        os.system('cp '+Config.PATH_SPA+'/d18O_soil1.map '+Config.PATH_SPA+'/d18O.L1.map')
#        os.system('cp '+Config.PATH_SPA+'/d18O_soil2.map '+Config.PATH_SPA+'/d18O.L2.map')
#        os.system('cp '+Config.PATH_SPA+'/d18O_soil3.map '+Config.PATH_SPA+'/d18O.L3.map')
#        os.system('cp '+Config.PATH_SPA+'/d18O_groundwater.map '+Config.PATH_SPA+'/d18O.GW.map')
        
#        os.system('cp '+Config.PATH_SPA+'/Age_snowpack.map '+Config.PATH_SPA+'/Age.snowpack.map') 
#        os.system('cp '+Config.PATH_SPA+'/Age_surface.map '+Config.PATH_SPA+'/Age.surface.map')
#        os.system('cp '+Config.PATH_SPA+'/Age_soil1.map '+Config.PATH_SPA+'/Age.L1.map')
#        os.system('cp '+Config.PATH_SPA+'/Age_soil2.map '+Config.PATH_SPA+'/Age.L2.map')
#        os.system('cp '+Config.PATH_SPA+'/Age_soil3.map '+Config.PATH_SPA+'/Age.L3.map')
#        os.system('cp '+Config.PATH_SPA+'/Age_groundwater.map '+Config.PATH_SPA+'/Age.GW.map')
            
print 'Total number of variables :', Opti.nvar

#############################################################################################
# Calibration loop
# -------------------

if Config.mode == 'calib_runs':

    print 'Number of iterations      :', Opti.nit
    if Config.restart == 1:
        ECH2O_Tools.restart(Config, Opti, Data)
        print '...but directly restarting from iter. ',Config.itres
    print
    print '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
    print ' Entering the optimisation loop...'
    print


    # Output initialization
    Config.initpar = 0
    Config.initobs = 0
    Opti.begfail = 0

    if Config.restart == 1:
        it0 = Config.itres-1
    else:
        it0 = 0

    for it in range(it0,Opti.nit):
    
        Opti.itout = '%03i' % int(it+1)

        if len(glob.glob(Config.PATH_EXEC))==0:
            os.system('mkdir '+Config.PATH_EXEC)
        else:
            os.system('rm -f '+Config.PATH_EXEC+'/*')
        print 'Iteration',Opti.itout, 'of',Opti.nit

        # Create the inputs for ECH2O
        os.system('cp -f '+Config.PATH_CLIM_REF_P+'/d2H_'+str(it+1)+'.txt '+Config.PATH_CLIM_REF+'/d2H.txt')
        os.system(Config.PATH_CLIM_REF+'/asc2c '+ Config.PATH_CLIM_REF+'/d2H.txt ' + Config.PATH_CLIM_REF + '/d2H.bin')
        ECH2O_Tools.create_inputs(Opti, Paras, Site, Config, it)
        # Run ECH2O
        os.chdir(Config.PATH_OUT)
        print '--> running ECH2O'
        start = time.time()
        os.system(Config.cmde_ech2o+' > ech2o.log')
        print '    run time:',time.time() - start,'seconds (limit at '+Config.tlimit+')'
        # Check if it ran properly
        os.chdir(Config.PATH_EXEC)
        if ECH2O_Tools.runOK(Data, Opti, Config) == 0:
            # If the run fails, let's give it one more chance!
            os.chdir(Config.PATH_OUT)
            os.system('rm -f '+Config.PATH_EXEC+'/*')
            print '--> running ECH2O'
            start = time.time()
            os.system(Config.cmde_ech2o+' > '+Config.PATH_EXEC+'/ech2o.log')
            print '    run time:',time.time() - start,'seconds (limit at '+Config.tlimit+')'
            os.chdir(Config.PATH_EXEC)
            # Still not running properly? Report
            if ECH2O_Tools.runOK(Data, Opti, Config) == 0:
                f_failpar = Config.PATH_OUT+'/Parameters_fail.txt'
                if len(glob.glob(f_failpar))==0:
                    with open(f_failpar,'w') as f_in:
                        f_in.write('Sample,'+','.join(Opti.names)+'\n')
                        f_in.write(str(it+1)+','+','.join([str(x) for x in Opti.x])+'\n')
                else:
                    with open(f_failpar,'a') as f_in:
                        f_in.write(str(it+1)+','+','.join([str(x) for x in Opti.x])+'\n')
                os.system('mv '+Config.PATH_EXEC+'/ech2o.log '+Config.PATH_OUT+'/ech2o_'+Opti.itout+'.log')
                # If it is the very first iteration, record it for mana_outputs later
                if it==0:
                    Opti.begfail = 1

        # If it worked...
        if ECH2O_Tools.runOK(Data, Opti, Config) == 1:
            # Write parameters values for this sequence
            ECH2O_Tools.output_par(Opti, Config, it)
            # Group sampling outputs
            ECH2O_Tools.manage_outputs(Data, Opti, Config, it)

        os.chdir(Config.PATH_OUT)
        #sys.exit()
        # Clean up
        os.system('rm -f '+Config.PATH_EXEC+'/*')

#############################################################################################
# Single run with one parameter set
# -------------------

if Config.mode == 'forward_runs':

    print '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
    print ' Post-optim runs...'
    print

    for it in range(int(options.nEns)):

        print 'Iteration',str(it+1), 'of ',options.nEns

        # Create the inputs for ECH2O
        ECH2O_Tools.create_inputs(Opti, Paras, Site, Config, it)
        # Run ECH2O
        os.chdir(Config.PATH_OUT)
        print '--> running ECH2O'

        # Spinup ?
        if Config.spinup > 0:
            print '...spinup first...'
            ECH2O_Tools.spinup(Config)
            print "...now the 'real' run..."
        #sys.exit()

        start = time.time()
        os.system(Config.cmde_ech2o+' > '+Config.PATH_EXEC+'/ech2o.log')
        print '    run time:',time.time() - start,'seconds (limit at '+Config.tlimit+')'
        # Check if it ran properly
        os.chdir(Config.PATH_EXEC)
        if ECH2O_Tools.runOK(Data, Opti, Config) == 1:
            # Group outputs
            ECH2O_Tools.manage_outputs(Data, Opti, Config, it)
            # -- Report the full BasinSummary.txt files?
            if Config.repBS == 1:
                os.system('mv '+Config.PATH_EXEC+'/BasinSummary.txt '+Config.PATH_OUT+'/BasinSummary_run'+str(it+1)+'.txt')
                # os.system('rm -f *.tab')
            # Clean up
            os.system('rm -f '+Config.PATH_EXEC+'/*')
        
###############################################################################################
# Simulations when varying the parameters, Morris's one-at-a-time
# ---------------------------------------------------------------
if Config.mode == 'sensi_morris' and Config.MSinit == 0:

    print 'Total Number of iterations :', Opti.nruns
    print
    print '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
    print

    Config.initobs = 0
    Opti.begfail = 0

    print 
    print '======================================'
    print '## Runs along trajectory #'+Config.numsim
    print '--------------------------------------'

    # There are npara+1 runs for each trajectory
    # for irun in range(Opti.nvar+1):

        #runnb = '%02i' % int(irun+1)
    print 'Run '+str(irun+1)+' out of '+str('%02i' % int(Opti.nvar+1))

        #print
        #print '|- Creating parameter maps / table for this run...'

        # Create the inputs for ECH2O
    ECH2O_Tools.create_inputs(Opti, Paras, Site, Config, irun)

        # Run ECH2O
    os.chdir(Config.PATH_OUT)
    print '--> running ECH2O'
    start = time.time()
    os.system(Config.cmde_ech2o+' > '+Config.PATH_EXEC+'/ech2o.log')
    print '    run time:',time.time() - start,'seconds (limit at '+Config.tlimit+')'

        # Check if it ran properly
    os.chdir(Config.PATH_EXEC)
    if ECH2O_Tools.runOK(Data, Opti, Config) == 1:
            # Group outputs
        ECH2O_Tools.manage_outputs(Data, Opti, Config, irun)
            # -- Report the full BasinSummary.txt files?
        if Config.repBS == 1:
            os.system('mv '+Config.PATH_EXEC+'/BasinSummary.txt '+Config.PATH_OUT+'/BasinSummary_run'+str(irun+1)+'.txt')
                # os.system('rm -f *.tab')
        # Not running properly? Report
    else:
        f_failpar = Config.PATH_OUT+'/Parameters_fail.txt'
        if len(glob.glob(f_failpar))==0:
           with open(f_failpar,'w') as f_in:
                f_in.write('Sample,'+','.join(Opti.names)+'\n')
                f_in.write(str(irun+1)+','+','.join([str(x) for x in Opti.x])+'\n')
        else:
            with open(f_failpar,'a') as f_in:
                f_in.write(str(irun+1)+','+','.join([str(x) for x in Opti.x])+'\n')
           # If it is the very first iteration, record it for mana_outputs later
        if irun==0:
            Opti.begfail = 1

        os.system('mv '+Config.PATH_EXEC+'/ech2o.log '+Config.PATH_OUT+'/ech2o_'+str(irun+1)+'.log')
            
        # Clean up
    os.system('rm -f '+Config.PATH_EXEC+'/*')

    # print Data.obs
    # Only for debugging ---------------------------------------------------------
#for oname in Data.names:
#    if Data.obs[oname]['type']!='map':
#        Data.obs[oname]['sim_hist'] = Config.PATH_OUT+'/'+oname+'_all.tab'
    # ----------------------------------------------------------------------------

    # Calculate and output the elementary effects
#ECH2O_Tools.morris_ee(Config, Data, Opti)

###############################################################################################
# END PROGRAM MAIN

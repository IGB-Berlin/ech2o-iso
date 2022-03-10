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
# Only generating parameters ?
parser.add_option("--MCmode",dest="MCmode",metavar="MCmode",help="Switch (0=init, 1=optim, 2=post optim runs)")
# Only generating parameters ?
parser.add_option("--tlimit",dest="tlimit",metavar="tlimit",help="Time limit of one ECH2O run (in seconds)")
# Input dir (if post optim) ?
parser.add_option("--inprf",dest="inprf",metavar="inprf",help="Input prefix")
# Input dir (if post optim) ?
parser.add_option("--nbest",dest="nbest",metavar="nbest",help="Number of best runs")


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

# -- Number of CPUs used in parallel
if options.ncpu == None:
    options.ncpu = 1
Config.ncpu = copy.copy(options.ncpu)

# -- MC init: many things don't happen if it is the case
if options.MCmode == None:
    sys.exit("Choose if you're initializating the MC sampling (or not)!")
else:
    Config.MCmode = int(options.MCmode)

# -- Output path
if options.outdir == None:
    options.outdir = os.path.splitext(options.file)[0]

# -- Time wall for ECH2O execution
if options.tlimit == None:
    Config.tlimit = ''
else:
    Config.tlimit = options.tlimit
    Config.tcmd = 'ulimit -t '+str(int(options.tlimit)*int(Config.ncpu))+' ;'

#print options.outdir
##
# -- all parameter path
if Config.MCmode<2:
    Config.PATH_PAR  = os.path.abspath(os.path.join(Config.PATH_MAIN,'Parameters_samples'))
    Config.FILE_PAR = Config.PATH_PAR+'/'+options.outdir.split('.')[0]+'_parameters.'
    # -- Creation of output directory
    if len(glob.glob(Config.PATH_PAR)) == 0: os.system('mkdir '+ Config.PATH_PAR)
    # -- Some verbose
    print
    print "Parameter samples' directory:          ", Config.PATH_PAR
    print

##
if Config.MCmode>0:
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

    # -- Execution command: time limit and OMP use
    Config.cmde_ech2o = ' '.join([Config.tcmd,
                                  'OMP_NUM_THREADS='+str(Config.ncpu),
                                  '../'+exe_ech2o,'../'+cfg_ech2o])

    # -- Get the parallel job number, based on the output dir name
    Config.numsim = options.outdir.split('.')[::-1][0]

    # --- Defining the various PATHs
    Config.PATH_OUT  = os.path.abspath(os.path.join(Config.PATH_MAIN,options.outdir))
    Config.PATH_SPA = os.path.abspath(os.path.join(Config.PATH_OUT,'Spatial'))

if Config.MCmode == 1:

    Config.PATH_EXEC = os.path.abspath(os.path.join(Config.PATH_OUT,'tmp'))    

if Config.MCmode == 2:

    Config.PATH_PAR  = os.path.abspath(os.path.join(Config.PATH_MAIN,'Outputs'))
    Config.FILE_PAR = Config.PATH_PAR+'/'+options.inprf+'.'+options.nbest+'bestParams.txt'

if Config.MCmode > 0 :

    Config.PATH_CLIM = os.path.abspath(os.path.join(Config.PATH_MAIN,'..','Climate'))
    # where the base maps are first copied from
    Config.PATH_SPA_REF = os.path.abspath(os.path.join(Config.PATH_MAIN,'Spatial')) 

    # -- Creation of output directory
    if len(glob.glob(Config.PATH_OUT)) == 0: os.system('mkdir '+ Config.PATH_OUT)
    # if len(glob.glob(Config.PATH_EXEC)) == 0: os.system('mkdir '+ Config.PATH_EXEC)

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
    os.system('cp -p '+Config.PATH_SPA_REF+'/*.map '+Config.PATH_SPA)

    # -- Copy def file in the output
    os.system('cp '+file+' '+Config.PATH_OUT)
    
#############################################################################################
# Initialization
# ------------------------------

#-- Import classes and setup from the def file
sys.path.insert(0,Config.PATH_MAIN)
exec('from ' + file_py + ' import *')

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
Paras.ind = {}
ipar=0
ipar2=0
Paras.isveg = 0
for par in Paras.names:
        
    # Dimensions (soil or veg or 1)
    nr = Paras.ref[par]['soil']*(Site.ns-1)+Paras.ref[par]['veg']*(Site.nv-1)+1
    # Build vectors used in the optimisation
    if type(Paras.ref[par]['min'])== float:
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
    elif Paras.ref[par]['veg']==1:
        Opti.names = Opti.names + [par + '_' + s for s in Site.vegs]
        Paras.isveg += 1
    else:
        Opti.names = Opti.names + [par]
    
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
Opti.nvar = len(Opti.min)

if Config.MCmode == 0:

    # Total number of samples
    Opti.nsamptot = Opti.nit*int(Config.ncpu)
    # Generate enough parameters sets and save them
    MC_Tools.gen_paras(Opti,Config, 'uniform')
    # Done
    print 'Parameters sample generation done.'

if Config.MCmode == 1:

    print 'Get parameters samples for this job...'
    # Read parameters sample for this parallel run
    MC_Tools.get_par(Opti, Config)

if Config.MCmode > 0:

    # Get the observations used for optimisation
    Data.names = Data.obs.keys()
    Data.names.sort()

    # Remove the default map files of calibrated param in the inputs directory
    # --> helps checking early on if there is an improper map update
    for pname in Paras.names:
        if Paras.ref[pname]['veg']==0:
            os.system('rm -f '+Config.PATH_SPA+'/'+Paras.ref[pname]['file']+'.map')

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

# IN post-optim mode, read directly the params from the file
if Config.MCmode == 2:

    print 'Get parameter set for this post-optim run...'

    tmp = list(np.genfromtxt(Config.FILE_PAR,delimiter=',',dtype= '|S20')[0])
    if tmp != Opti.names:
        sys.exit("The definition file and input parameter file ain't matching!")

    #pnames = ','.join([str(tmp[idx]) for idx in range(Opti.nvar)])
    #print tmp#pnames
    #print Opti.names
    Opti.xpar = np.genfromtxt(Config.FILE_PAR,delimiter=',',skip_header=1)
    #Opti.x = [tmp[int(Config.numsim)-1][x] for x in range(Opti.nvar)]

print 'Total number of variables :', Opti.nvar

#############################################################################################
# Calibration loop
# -------------------

if Config.MCmode == 1:

    print 'Number of iterations      :', Opti.nit
    print
    print '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
    print ' Entering the optimisation loop...'
    print


    # Output initialization
    Config.initpar = 0
    Config.initobs = 0

    for it in range(Opti.nit):
    
        Opti.itout = '%03i' % int(it+1)

        if len(glob.glob(Config.PATH_EXEC))==0:
            os.system('mkdir '+Config.PATH_EXEC)
        else:
            os.system('rm -f '+Config.PATH_EXEC+'/*')
        print 'Iteration',Opti.itout, 'of',Opti.nit

        # Create the inputs for ECH2O
        MC_Tools.create_inputs(Opti, Paras, Site, Config, it)
        # Run ECH2O
        os.chdir(Config.PATH_EXEC)
        print '--> running ECH2O'
        start = time.time()
        os.system(Config.cmde_ech2o+' > ech2o.log')
        print '    run time:',time.time() - start,'seconds (limit at '+Config.tlimit+')'
        # Check if it ran properly
        if MC_Tools.runOK(Data, Opti, Config) == 1:
            # Write parameters values for this sequence
            MC_Tools.output_par(Opti, Config, it)
            # Group sampling outputs
            MC_Tools.manage_outputs(Data, Opti, Config, it)
        else:
            f_failpar = Config.PATH_OUT+'/Parameters_fail.txt'
            if len(glob.glob(f_failpar))==0:
                with open(f_failpar,'w') as f_in:
                    f_in.write('Sample,'+','.join(Opti.names)+'\n')
                    f_in.write(str(it+1)+','+','.join([str(x) for x in Opti.x])+'\n')
            else:
                with open(f_failpar,'a') as f_in:
                    f_in.write(str(it+1)+','+','.join([str(x) for x in Opti.x])+'\n')
            os.system('mv '+Config.PATH_EXEC+'/ech2o.log '+Config.PATH_OUT+'/ech2o_'+Opti.itout+'.log')
            #os.system('mv '+Config.PATH_EXEC+' '+Config.PATH_EXEC+'_'+Opti.itout)
        os.chdir(Config.PATH_OUT)
        #sys.exit()
        os.system('rm -fr '+Config.PATH_EXEC)
            #os.system('rm -fr '+Config.PATH_EXEC+'/*')

#############################################################################################
# Single run with one parameter set
# -------------------

if Config.MCmode == 2:

    print '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
    print ' Post-optim runs...'
    print

    for it in range(int(options.nbest)):

        print 'Iteration',str(it+1), 'of ',options.nbest

        Config.PATH_EXEC = os.path.abspath(os.path.join(Config.PATH_OUT,str(it+1)))
        if len(glob.glob(Config.PATH_EXEC))==0:
            os.system('mkdir '+Config.PATH_EXEC)
        else:
            os.system('rm -f '+Config.PATH_EXEC+'/*')

        # Create the inputs for ECH2O
        MC_Tools.create_inputs(Opti, Paras, Site, Config, it)
        # Run ECH2O
        os.chdir(Config.PATH_EXEC)
        print '--> running ECH2O'
        start = time.time()
        os.system(Config.cmde_ech2o+' > ech2o.log')
        print '    run time:',time.time() - start,'seconds (limit at '+Config.tlimit+')'
        # Check if it ran properly
        if MC_Tools.runOK(Data, Opti, Config) == 1:
            # Group outputs
            MC_Tools.manage_outputs2(Data, Opti, Config, it)
        
###############################################################################################
###############################################################################################

#    !/usr/bin/env python3
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
# Number of iterations
parser.add_option("--nit",dest="nit",metavar="nit",help="Number of iterations")

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
if options.pathout == None:
    options.pathout = os.getcwd()+'/'+os.path.splitext(options.file)[0]

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
Config.PATH_OUT  = os.path.abspath(options.pathout)
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

#############################################################################################
# Calibration "model" (loop)
# Actually we a 1-step dynamic model calling ECH2O...
# -------------------
class IsoModel(DynamicModel, MonteCarloModel):

    # Initialize the sampling
    def __init__(self,Opti,Paras,Config,Site):

        DynamicModel.__init__(self)
        MonteCarloModel.__init__(self)

        # Copy variables, dictionaries, etc.
        self.ref    = copy.copy(Paras.ref)
        self.obs    = copy.copy(Data.obs)

        # Paths
        self.PATH_MAIN    = Config.PATH_MAIN
        self.PATH_SPA     = Config.PATH_SPA
        self.PATH_SPA_REF = Config.PATH_SPA_REF
        self.PATH_OUT     = Config.PATH_OUT

        # Catchement properties
        self.soils = Site.soils
        self.ns    = len(self.soils)
        self.sfiles= Site.sfiles
        self.vegs  = Site.vegs
        self.nv    = len(Site.vegs)
        self.vfile = Site.vfile

        # Number of iterations
        self.nit = Opti.nit

        # Set maps with coord etc.
        setclone(self.PATH_SPA+'/unit.map')
         
    # Parameters, inputs or variables which are constant
    def premcloop(self):

        # Parameters: from definition to all values
        self.pnames = self.ref.keys()
        self.pn     = len(self.pnames)

        # Remove the default map files of calibrated param in the inputs directory
        # --> helps checking early on if there is an improper map update
        for pname in self.pnames:
            if self.ref[pname]['veg']==0:
                os.system('rm -f '+self.PATH_SPA+'/'+self.ref[pname]['file']+'.map')

        # Read dictionary to get all param setup
        self.xmin = []
        self.xmax = []
        self.xlog = []
        self.xnames = []
        self.xind = []
        self.pind = {}
        ipar=0
        ipar2=0
        self.isveg = 0
        for par in self.pnames:
            # Dimensions (soil or veg or 1)
            nr = self.ref[par]['soil']*(self.ns-1)+self.ref[par]['veg']*(self.nv-1)+1
            # Build vectors used in the optimisation
            self.xmin.append(np.repeat(self.ref[par]['min'],nr))
            self.xmax.append(np.repeat(self.ref[par]['max'],nr))
            self.xlog.append(np.repeat(self.ref[par]['log'],nr))
            # Link betwen params and all variables
            self.xind.append(np.repeat(ipar,nr))
            ipar+=1
            if nr>1: self.pind[par] = list(np.arange(ipar2,ipar2+nr,1))
            if nr==1: self.pind[par] = ipar2
            ipar2+=nr
            # For outputs
            if self.ref[par]['soil']==1:
                self.xnames.append([par + '_' + s for s in self.soils])
            elif self.ref[par]['veg']==1:
                self.xnames.append([par + '_' + s for s in self.vegs])
                self.isveg += 1
            else:
                self.xnames.append([par])
    
        self.xmin   = list(chain(*self.xmin))
        self.xmax   = list(chain(*self.xmax))
        self.xlog   = list(chain(*self.xlog))
        self.xnames = list(chain(*self.xnames))
        self.xind   = list(chain(*self.xind))

        # Total number of variables
        self.xn = len(self.xmin)

        # Soils / units maps
        self.bmaps = {}
        for im in range(self.ns):
            self.bmaps[self.soils[im]] = readmap(self.PATH_SPA+'/'+self.sfiles[im])
        self.bmaps['unit'] = readmap(self.PATH_SPA+'/unit.map')
        self.bmaps['chanmask'] = readmap(self.PATH_SPA+'/chanmask.map')
        self.bmaps['chanmask_NaN'] = readmap(self.PATH_SPA+'/chanmask_NaN.map')

        # Reference dictionary for vegetation inputs file
        self.vref = {}
        with open(self.PATH_SPA_REF+'/'+self.vfile,'r') as csvfile:
            paramread = list(csv.reader(csvfile, delimiter='\t'))
        exit
        # "Head": number of species and of params
        self.vref['header'] = paramread[0][0:len(paramread[0])]
        # All parameters values (keep strings!)
        for iv in range(self.nv):
            self.vref[iv] = paramread[iv+1][0:len(paramread[iv+1])]
        # "Footers" : name of head1 and of parameters
        self.vref['footer'] = paramread[self.nv+1][0:len(paramread[self.nv+1])]
        self.vref['name'] = paramread[self.nv+2][0:len(paramread[self.nv+2])]

        # -----------------------------------------------------------
        print 'Total number of variables :', self.xn
        print 'Number of iterations      :', self.nit
        print
        print '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
        print ' Entering the optimisation loop...'
        print

    ###########################################    
    # What you change between sampling
    def initial(self):

        # --------------------------------------------------------------------------
        # Create a temporary + Spatial directory so that parallel runs don't mess up
        self.itfmt = '%06i' % int(self.currentSampleNumber())
        #print self.itfmt

        self.PATH_TMP = self.PATH_OUT+'/'+self.itfmt
        if len(glob.glob(self.PATH_TMP)) != 0:
            os.system('rm -fr '+self.PATH_TMP)
        os.system('mkdir '+self.PATH_TMP)

        PATH_SPA_TMP = self.PATH_TMP+'/Spatial'
        if len(glob.glob(PATH_SPA_TMP)) != 0:
            os.system('rm -fr '+PATH_SPA_TMP)
        os.system('cp -pr '+self.PATH_SPA+' '+PATH_SPA_TMP)

        # -----------------------
        os.chdir(self.PATH_TMP)

        # --------------------------------------------------------------------------
        # Parameters values
        self.x = np.zeros(self.xn)
        for i in range(self.xn):
            # log transform
            if self.xlog[i]==1:
                self.x[i] = np.log10(np.random.uniform(10**self.xmin[i],10**self.xmax[i]))
            else:
                self.x[i] = np.random.uniform(self.xmin[i],self.xmax[i])
            # print i, self.xnames[i], self.x[i]
                
        # Output parameters record
        f_par = open('parameters.txt','w')
        #print self.xnames
        #print ','.join(self.xnames)
        #print ','.join(str(self.x))+'\n'
        #f_par.write(','.join(self.xnames)+'\n')
        #f_par.write(','.join(str(list(self.x)))+'\n')
        csv.writer(f_par,self.xnames)
        csv.writer(f_par,self.x)
        f_par.close()
       
        # -------------------------------------------------------------------------------
        # Write inputs files for model runs
        for pname in self.pnames:

            # - Mapped parameters
            if self.ref[pname]['veg']==0:

                # Soil unit dependence
                if self.ref[pname]['soil']==1 :
                    # rint 'Soil dependent !!'
                    outmap = self.bmaps['unit']*0
                    # Read each soil map unit and apply param value
                    for im in range(self.ns):
                        outmap+= self.bmaps[self.soils[im]]*self.x[self.pind[pname][im]]

                # No spatial/veg dependence, but channel stuff
                else:
                    # print 'Not dependent !!'
                    if self.ref[pname]['file'] in ['chanwidth','chanmanningn']:
                        outmap = self.bmaps['chanmask']*self.x[self.pind[pname]]
                    elif Paras.ref[pname]['file'] == 'chanparam':
                        outmap = self.bmaps['chanmask_NaN']*self.x[self.pind[pname]]
                    else:
                        outmap = self.bmaps['unit']*self.x[self.pind[pname]]
    
                # Outputs updated map
                report(outmap,PATH_SPA_TMP+'/'+self.ref[pname]['file']+'.map')
        
                # Check for initial condition/other parameter dependence
                # Initial soil water content
                if self.ref[pname]['file']=='poros':
                    report(outmap*0.8,PATH_SPA_TMP+'/SWC1.map')
                    report(outmap*0.8,PATH_SPA_TMP+'/SWC2.map')
                    report(outmap*0.8,PATH_SPA_TMP+'/SWC3.map')
                # Rootin layer 2: assuming that all roots are within the first two hydraulics layers
                if self.ref[pname]['file']=='rootfrac1':
                    base = readmap(PATH_SPA_TMP+'/unit.map')
                    report(base - outmap, PATH_SPA_TMP+'/rootfrac2.map')

            # - Vegetation parameters
            else:
                # Change the value based in the dict based on name correspondance
                vegnew = copy.copy(self.vref)
                for iv in range(self.nv):
                    vegnew[iv][vegnew['name'].index(pname)] = str(self.x[self.pind[pname][iv]])
        
        # ------------------------------------------------------------------------------
        # Finalizing the preps....

        # - Finalizing soil parameterization
        # Check that initial soil moisture is not smaller residual soil 
        # tbd... for now just pick the porosity and theta_r ranges wisely enough
            
        # - Finalizing the vegetation parameterization
        if self.isveg > 0:
            # Equalize leaf turnover and additional turnover due to water and/or temperature stree
            for iv in range(self.nv):
                vegnew[iv][vegnew['name'].index('TurnovL_MWS')] = copy.copy(vegnew[iv][vegnew['name'].index('TurnovL')])
                vegnew[iv][vegnew['name'].index('TurnovL_MTS')] = copy.copy(vegnew[iv][vegnew['name'].index('TurnovL')])
            # Write the vegetation params file (if there's at least one veg param)
            vegfile = open(PATH_SPA_TMP+'/'+self.vfile,'w')
            vegfile.write('\t'.join(self.vref['header'])+'\n')
            for iv in range(self.nv):
                vegfile.write('\t'.join(vegnew[iv]) +'\n')
            vegfile.write('\t'.join(self.vref['footer'])+'\n')
            vegfile.write('\t'.join(self.vref['name'])+'\n')
            vegfile.close()
        
    def dynamic(self):

        self.setQuiet(quiet=True)
        # Run ECH2O
        os.chdir(self.PATH_TMP)
        print 'Iteration',str(self.currentSampleNumber()), 'of',self.nit        
        os.system('export OMP_NUM_THREADS=1')
        os.system(cmde_ech2o + ' > ech2o.log')
        os.system('rm -fr '+self.PATH_MAIN+'/'+str(self.currentSampleNumber()))

##############################################################################################
myModel = IsoModel(Opti,Paras,Config,Site)
dynamicModel = DynamicFramework(myModel, 1, 1)
#staticModel = StaticFramework(myModel)
mcModel = MonteCarloFramework(dynamicModel, nrSamples=Opti.nit)
#mcModel = MonteCarloFramework(staticModel, nrSamples=Opti.nit)

mcModel.setForkSamples(fork=True,nrCPUs=30)
mcModel.run(postmc=False)
###############################################################################################
###############################################################################################

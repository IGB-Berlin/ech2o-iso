#/usr/bin/env python
# -*- coding: utf-8 -*-
# *************************************************
#
# Morris sensitivity analysis for ECH2O
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
import csv

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
# Time limit
parser.add_option("--tlimit",dest="tlimit",metavar="tlimit",help="Time limit of one ECH2O run (in seconds)")

(options, args) = parser.parse_args()


# set working directory
class Config:
    PATH_MAIN = os.getcwd()+'/'

print '*********************************************************************************'
print 'MORRIS SENSITIVITY TEST : '
print '- construction of the trajectories'
print '- forward runs'
print '- storage of outputs and info for posterior analysis : elementary effects, etc.' 
print '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
print

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

# -- Output path
if options.outdir == None:
    options.outdir = os.path.splitext(options.file)[0]

# -- Time wall for ECH2O execution
if options.tlimit == None:
    Config.tlimit = ''
else:
    Config.tlimit = options.tlimit
    Config.tcmd = 'ulimit -t '+str(int(options.tlimit)*int(Config.ncpu))+' ;'

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
Config.PATH_EXEC = os.path.abspath(os.path.join(Config.PATH_OUT,'tmp'))
Config.PATH_SPA = os.path.abspath(os.path.join(Config.PATH_OUT,'Spatial'))
Config.PATH_CLIM = os.path.abspath(os.path.join(Config.PATH_MAIN,'..','Climate'))
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
os.system('cp -p '+Config.PATH_SPA_REF+'/*.map '+Config.PATH_SPA)
    
#############################################################################################
# Initialization
# ------------------------------

#-- Import classes and setup from the def file
exec('from ' + file_py + ' import *')

#-- Parameters: from definition to all values
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

    print par

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
        Opti.comp = Opti.comp + [i for i in range(Site.ns)]
    elif Paras.ref[par]['veg']==1:
        Opti.names = Opti.names + [par + '_' + s for s in Site.vegs]
        Opti.comp = Opti.comp + [i for i in range(Site.nv)]
        Paras.isveg += 1
    else:
        Opti.names = Opti.names + [par]
        Opti.comp = Opti.comp + [0]


print Opti.names
print Opti.comp

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

# Add a symbolic link to armadillo
#libarma = 'libarmadillo.so.7'
#if len(glob.glob(os.path.join(PATH_EXEC,libarma))) != 0:
#    os.system('rm '+os.path.join(PATH_EXEC,libarma))
#os.symlink(os.path.join(PATH_MAIN, libarma) , os.path.join(PATH_EXEC,libarma) )

# Total number of runs
nruns = (Opti.nvar+1) * Opti.nr
print
print 'REMINDER --> Total number of runs : ', str(nruns)
#sys.exit()

print 'Total number of variables :', Opti.nvar
print 'Total Number of iterations      :', nruns
print
print '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
print


# -----------------------------------------------------
# --- Creation of the Morris trajectories ---
# -----------------------------------------------------

Opti.vals={}

# Value possible for each parameter
Opti.step = np.zeros((Opti.nvar),np.float64)

for i in range(Opti.nvar):
    if Opti.log[i]==0: # Linear
        Opti.step[i] = (Opti.max[i]-Opti.min[i])/(Opti.nlev-1)
        Opti.vals[i] = np.arange(Opti.min[i],Opti.max[i]+Opti.step[i], Opti.step[i], np.float64)
    if Opti.log[i]==1: # Log
        Opti.step[i] = (np.log10(Opti.max[i])-np.log10(Opti.min[i]))/(Opti.nlev-1)
        Opti.vals[i] = 10**(np.arange(np.log10(Opti.min[i]),np.log10(Opti.max[i])+Opti.step[i], Opti.step[i], np.float64))
    
# Construct B* matrix, for each repetition
Opti.Bstar = np.zeros((Opti.nvar,Opti.nvar+1,Opti.nr),np.float64)
Opti.Xstar = np.zeros((Opti.nvar),np.float64)
Opti.ind_hist = np.zeros((Opti.nvar,Opti.nr)) #history of the index of trajectory construction

for ir in range(Opti.nr):

    # Generate 'base' vector X*
    ind_tmp = np.zeros((Opti.nvar))
    for i in range(Opti.nvar):
        ind_tmp[i] = random.randint(0,Opti.nlev-1)
        Opti.Xstar[i] = Opti.vals[i][ind_tmp[i]]
        #if Opti.log[i]==0: # Linear variation
            #Opti.Xstar[i] = Opti.min[i]+ind_tmp[i]*Opti.step[i]
        if Opti.Xstar[i] < Opti.min[i]-np.abs(Opti.min[i])*0.000001:
            print Opti.names[i],Opti.Xstar[i], Opti.min[i],ind_tmp[i],Opti.step[i]#,mins[i]+ind_tmp[i]*step[i]
            sys.exit()
        #else: # Log variation
        #    Opti.Xstar[i]=10**(np.log10(Opti.min[i])+ind_tmp[i]*Opti.step[i])
            
        #    if Opti.Xstar[i] < Opti.min[i]-np.abs(Opti.min[i])*0.000001:
        #        print Opti.names[i],Opti.Xstar[i], Opti.min[i],ind_tmp[i],Opti.step[i]#,mins[i]+ind_tmp[i]*step[i]
        #        sys.exit()

    # Generate trajectory
    # first: the first point comes from an increase by some steps of a random number of components
    Opti.Bstar[:,0,ir] = Opti.Xstar

    # then: all the next are changed by only one component with + or - one step
    #if ir==0:
    ind2 = range(Opti.nvar)
    random.shuffle(ind2)
    Opti.ind_hist[:,ir]=ind2
    for i in range(Opti.nvar):
        #print varnames[ind2[i]]
        Opti.Bstar[:,i+1,ir] = Opti.Bstar[:,i,ir] #copy previous location
        eps = random.choice([-1,1])

        if Opti.log[ind2[i]]==0: #linear
            tmp = Opti.Bstar[ind2[i],i,ir] + Opti.step[ind2[i]]*eps #add (or substract) one step for a random component (not used)
            if ind_tmp[ind2[i]] == 0: # case where we already are on the lower boundary : add
                tmp = Opti.Bstar[ind2[i],i,ir] + Opti.step[ind2[i]]
            if ind_tmp[ind2[i]] == Opti.nlev-1 : # case where we already are on the upper boundary : sbtrkt
                tmp = Opti.Bstar[ind2[i],i,ir] - Opti.step[ind2[i]]

        if Opti.log[ind2[i]]==1: #log
            tmp = 10**(np.log10(Opti.Bstar[ind2[i],i,ir]) + Opti.step[ind2[i]]*eps) #add (or substract) one step for a random component (not used)
            if ind_tmp[ind2[i]] == 0: # case where we already are on the lower boundary : add
                tmp = 10**(np.log10(Opti.Bstar[ind2[i],i,ir]) + Opti.step[ind2[i]])
            if ind_tmp[ind2[i]] == Opti.nlev-1 : # case where we already are on the upper boundary : sbtrkt
                tmp = 10**(np.log10(Opti.Bstar[ind2[i],i,ir]) - Opti.step[ind2[i]])

                
        Opti.Bstar[ind2[i],i+1,ir] = tmp
        if tmp>Opti.max[ind2[i]]+np.abs(Opti.max[ind2[i]])*0.000001 or tmp<Opti.min[ind2[i]]-np.abs(Opti.min[ind2[i]])*0.000001:
            print 'Error in the incrementation of the parameter '+Opti.names[ind2[i]]
            print Opti.Xstar[ind2[i]],Opti.Bstar[ind2[i],i,ir],Opti.Bstar[ind2[i],i+1,ir],tmp
            print Opti.min[ind2[i]], Opti.max[ind2[i]],Opti.step[ind2[i]]#,ind2_tmp[ind2[i]],gap_tmp[ind2[i]]
            sys.exit()

    #print Bstar[ind2[i],i,ir], Bstar[ind2[i],i+1,ir]
    #print ind2
    #print ind_hist[:,ir]
        

# -------------------------------------
#- Write parameters values file

# -- Change directory
os.chdir(Config.PATH_OUT)
#PATH_EXEC = os.path.abspath('')
#print
#print 'Current directory : ', PATH_OUT

# -- Write trajectories

# Each parameters in a file
for ip in range(Opti.nvar):
    with open('trajectory_'+Opti.names[ip]+'.csv','wb') as csvfile:
        csv_writer = csv.writer(csvfile)
        for ir in range(Opti.nr):
            csv_writer.writerow(Opti.Bstar[ip,:,ir])
    exit

# Altogether (not possible for now with csv, will need to get netCDF...)
#Vars['Bstar'] = {'dim_name': ('npara','nsteps','nrepet',), \
#               'dim_size': (npara,npara+1,nr),\
#               'attr_name' : ['long_name'],\
#               'attr_value' : ['Matrix of Morris parameter trajectories'],\
#               'value' : Bstar}
# Write Bstar for each trajectory
for ir in range(Opti.nr):
    trajnb = '%02i' % int(ir+1)
    #print trajnb
    with open('Bstar_traj'+trajnb+'.csv','wb') as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(Opti.names)
        for irun in range(Opti.nvar+1):
            csv_writer.writerow(Opti.Bstar[:,irun,ir])
    exit


# History of the component changed at each step of the trajectories
with open('trajectory_history.csv','wb') as csvfile:
    csv_writer = csv.writer(csvfile)
    for ir in range(Opti.nr):
        csv_writer.writerow(Opti.ind_hist[:,ir])
exit

# -----------------------------------------------
# --- Simulations when varying the parameters ---
# -----------------------------------------------

# Each trajectory
for ir in range(Opti.nr):
    trajnb = '%02i' % int(ir+1)
    print 
    print '======================================'
    print '## Trajectory '+trajnb+' out of '+str('%02i' % int(Opti.nr))
    print '----------------------------------'

    vrun = {}
    
    # There are npara+1 runs for each trajectory
    for irun in range(Opti.nvar+1):

        # -- Prepare the runs : constructs maps and/or parameter file with
        #                       initialized/updated parameters
        os.chdir(Config.PATH_SPA)
                
        runnb = '%02i' % int(irun+1)
        print 'Run '+str(irun+1)+' out of '+str('%02i' % int(Opti.nvar+1))

        if irun == 0:

            print
            print '|- Initializing parameter maps / table for this trajectory...'
            print

            # Assign the first vector parameters values
            for pname in Paras.names:

                ## - Mapped parameters
                if Paras.ref[pname]['veg']==0:

                    # Soil unit dependence
                    if Paras.ref[pname]['soil']==1 and pname!='Rootp.1' and pname!='Rootp.1+2':
                        # print 'Soil dependent !!'
                        outmap = Site.bmaps['unit']*0
                        # Read each soil map unit and apply param value
                        for im in range(Site.ns):
                            outmap+= Site.bmaps[Site.soils[im]]*Opti.Bstar[Paras.ind[pname][im],irun,ir]

                        report(outmap,Config.PATH_SPA+'/'+Paras.ref[pname]['file']+'.map')

                    # No spatial/veg dependence, but channel stuff
                    elif pname!='Rootp.1+2' and pname!='Rootp.1':
                        # print 'Not dependent !!'
                        if Paras.ref[pname]['file'] in ['chanwidth','chanmanningn']:
                            outmap = Site.bmaps['chanmask']*Opti.Bstar[Paras.ind[pname],irun,ir]
                        elif Paras.ref[pname]['file'] == 'chanparam':
                            outmap = Site.bmaps['chanmask_NaN']*Opti.Bstar[Paras.ind[pname],irun,ir]
                        else:
                            outmap = Site.bmaps['unit']*Opti.Bstar[Paras.ind[pname],irun,ir]

                        report(outmap,Config.PATH_SPA+'/'+Paras.ref[pname]['file']+'.map')

                    # Root layer 2: summed fraction over 1+2 minus fraction over 1
                    elif pname=='Rootp.1':
                        outmap = Site.bmaps['unit']*0
                        outmap2 = Site.bmaps['unit']*0
                        for im in range(Site.ns):
                            tmp = Opti.Bstar[Paras.ind['Rootp.1'][im],irun,ir]
                            outmap+= Site.bmaps[Site.soils[im]]*Opti.Bstar[Paras.ind['Rootp.1+2'][im],irun,ir]*tmp
                            outmap2+= Site.bmaps[Site.soils[im]]*Opti.Bstar[Paras.ind['Rootp.1+2'][im],irun,ir]*(1-tmp)
                        # base = readmap(Config.PATH_SPA+'/unit.map')
                        report(outmap, Config.PATH_SPA+'/rootfrac1.map')
                        report(outmap2, Config.PATH_SPA+'/rootfrac2.map')

                    # Check for initial condition/other parameter dependence
                    # Initial soil water content
                    if pname=='Porosity':
                        report(outmap*0.8,Config.PATH_SPA+'/SWC1.map')
                        report(outmap*0.8,Config.PATH_SPA+'/SWC2.map')
                        report(outmap*0.8,Config.PATH_SPA+'/SWC3.map')

                # - Vegetation parameters
                else:
                    # Change the value based on param name correspondance
                    vegnew = copy.copy(Opti.vref)
                    # print Opti.vref
                    for iv in range(Site.nv):
                        vegnew[iv][vegnew['name'].index(pname)] = str(Opti.Bstar[Paras.ind[pname][iv],irun,ir])

            # ------------------------------------------------------------------------------
            # Finalizing the preps....

            # - Finalizing soil parameterization
            # Check that initial soil moisture is not smaller residual soil
            # tbd... for now just pick the porosity and thetar reange wisely enough

            # - Finalizing the vegetation parameterization
            if Paras.isveg > 0:
                # Equalize leaf turnover and additional turnover due to water and/or temperature stress
                # for iv in range(Site.nv):
                #    vegnew[iv][vegnew['name'].index('TurnovL_MWS')] = copy.copy(vegnew[iv][vegnew['name'].index('TurnovL')])
                #    vegnew[iv][vegnew['name'].index('TurnovL_MTS')] = copy.copy(vegnew[iv][vegnew['name'].index('TurnovL')])
                # Write the vegetation params file (if there's at least one veg param)
                vegfile = open(Config.PATH_SPA+'/'+Site.vfile,'w')
                vegfile.write('\t'.join(Opti.vref['header'])+'\n')
                for iv in range(Site.nv):
                    vegfile.write('\t'.join(vegnew[iv]) +'\n')
                vegfile.write('\t'.join(Opti.vref['footer'])+'\n')
                vegfile.write('\t'.join(Opti.vref['name'])+'\n')
                vegfile.close()


        if irun >= 1:

            # Change only the relevant parameter
            print
            print '|- Updating parameter maps / table within the trajectory...'

            ind_new = Opti.ind_hist[irun-1,ir]
            vold = Opti.Bstar[ind_new,irun-1,ir] # Back up old value
            vnew = Opti.Bstar[ind_new,irun,ir]   # Assign new value
            print '     Modif of '+Opti.names[ind_new]+' : '+str(vold)+' --> '+str(vnew)
            if vnew == vold: 
                sys.exit('ERROR : the wrong parameter was assumed as changed !!!')

            pname = Paras.names[Opti.ind[ind_new]]
        
            ## - Mapped parameters
            if Paras.ref[pname]['veg']==0:

                im = Opti.comp[ind_new]

                # Soil unit dependence
                if Paras.ref[pname]['soil']==1 and pname!='Rootp.1' and pname!='Rootp.1+2':
                    # Update map and re-write
                    inmap = readmap(Config.PATH_SPA+'/'+Paras.ref[pname]['file']+'.map')
                    outmap = inmap + Site.bmaps[Site.soils[im]]*(vnew-vold)
                    report(outmap,Config.PATH_SPA+'/'+Paras.ref[pname]['file']+'.map')

                # No spatial/veg dependence, but channel stuff
                elif pname!='Rootp.1+2' and pname!='Rootp.1':
                    # print 'Not dependent !!'
                    if Paras.ref[pname]['file'] in ['chanwidth','chanmanningn']:
                        outmap = Site.bmaps['chanmask']*vnew
                    elif Paras.ref[pname]['file'] == 'chanparam':
                        outmap = Site.bmaps['chanmask_NaN']*vnew
                    else:
                        outmap = Site.bmaps['unit']*vnew

                    report(outmap,Config.PATH_SPA+'/'+Paras.ref[pname]['file']+'.map')

                # Root layers 2: summed fraction over 1+2 minus fraction over 1
                elif pname=='Rootp.1' or pname=='Rootp.1+2':
                    inmap  = readmap(Config.PATH_SPA+'/rootfrac1.map')
                    inmap2 = readmap(Config.PATH_SPA+'/rootfrac2.map')
                    tmp = Opti.Bstar[Paras.ind['Rootp.1'][im],irun,ir]
                    if pname=='Rootp.1':
                        outmap = inmap + Site.bmaps[Site.soils[im]]*Opti.Bstar[Paras.ind['Rootp.1+2'][im],irun,ir]*(vnew-vold)
                        outmap2= inmap2 + Site.bmaps[Site.soils[im]]*Opti.Bstar[Paras.ind['Rootp.1+2'][im],irun,ir]*((1-vnew)-(1-vold))
                    if pname=='Rootp.1+2':
                        outmap = inmap + Site.bmaps[Site.soils[im]]*(vnew-vold)*Opti.Bstar[Paras.ind['Rootp.1'][im],irun,ir]
                        outmap2= inmap2 + Site.bmaps[Site.soils[im]]*(vnew-vold)*(1-Opti.Bstar[Paras.ind['Rootp.1'][im],irun,ir])
                    # base = readmap(Config.PATH_SPA+'/unit.map')
                    report(outmap, Config.PATH_SPA+'/rootfrac1.map')
                    report(outmap2, Config.PATH_SPA+'/rootfrac2.map')

                # Check for initial condition/other parameter dependence
                # Initial soil water content
                if pname=='Porosity':
                    report(outmap*0.8,Config.PATH_SPA+'/SWC1.map')
                    report(outmap*0.8,Config.PATH_SPA+'/SWC2.map')
                    report(outmap*0.8,Config.PATH_SPA+'/SWC3.map')

            # - Vegetation parameters
            else:
                # Reference dictionary for vegetation inputs file
                vegnew = {}
                with open(Config.PATH_SPA+'/'+Site.vfile,'r') as csvfile:
                    paramread = list(csv.reader(csvfile, delimiter='\t'))
                exit
                # Read previous parameters values (keep strings!)
                for iv in range(Site.nv):
                    vegnew[iv] = paramread[iv+1][0:len(paramread[iv+1])]
                # Change the value based on param name correspondance
                vegnew[Opti.comp[ind_new]][vegnew['name'].index(pname)] = str(vnew)
                vegfile = open(Config.PATH_SPA+'/'+Site.vfile,'w')
                vegfile.write('\t'.join(Opti.vref['header'])+'\n')
                for iv in range(Site.nv):
                    vegfile.write('\t'.join(vegnew[iv]) +'\n')
                vegfile.write('\t'.join(Opti.vref['footer'])+'\n')
                vegfile.write('\t'.join(Opti.vref['name'])+'\n')
                vegfile.close()

        #############################################################################################
        # Calibration loop
        # -------------------

        if len(glob.glob(Config.PATH_EXEC))==0:
            os.system('mkdir '+Config.PATH_EXEC)
        else:
            os.system('rm -f '+Config.PATH_EXEC+'/*')

        # -- Move to run directory and create links
        os.chdir(Config.PATH_EXEC)
        
        # -- Execute ECH2O --
        print
        print 'ECH2O simulation in progress (run ',str(irun+1),' of ',str(Opti.nvar+1),')...'
        runall = '%03i' % int(Opti.nr*ir+irun+1)
        #os.system(cmde_orchidee + ' > orchidee.log')

        #print cmde_ech2o
        #os.system(cmde_ech2o)
        os.system(Config.cmde_ech2o + ' > ech2o.'+runall+'.log')
        
        # -- Outputs reorganization
        # Basin summary
        if len(glob.glob(os.path.join(Config.PATH_EXEC,'BasinSummary.txt'))) == 0:
            print 'Error in the execution of ECH2O...'
            sys.exit()
        os.system('mv -f BasinSummary.txt BasinSummary_traj'+trajnb+'_run'+runnb+'.txt')      
        
        # Observation used for sensitivity
        for oname in range(Data.names):
            obs_out = os.path.join(Config.PATH_OUT,oname+'_traj'+trajnb+'_run'+runnb+'.tab')
            print 'Moving sensitivity outputs to : ',obs_out
            os.system('mv -f '+oname+'.tab '+obs_out)
        print
    #sys.exit()

# END PROGRAM MAIN

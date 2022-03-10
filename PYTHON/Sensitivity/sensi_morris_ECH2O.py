#!/usr/bin/env python
#*******************************************************************************
# PROGRAMME	: SENSI_MORRIS
# AUTEUR	: S. KUPPEL
# CREATION	: 01/2012
# COMPILATEUR	: PYTHON
#
# Analyse de sensibilite d'ORCHIDEE avec la methode de Morris
#
# Pour l'instant avec un seul PFT par site
#
#
#*******************************************************************************
"""

"""


# ==============================================================================
#  MAIN PROGRAM
# ==============================================================================

import time, os, glob, sys, copy
import csv
import ConfigParser
from optparse import OptionParser
import socket
from time import localtime, strftime
import random
import numpy as np
#import IO_tools
#import io as IO_tools

#os.system('./ech2o config_bk.ini')
#sys.exit()
# ---------
#  OPTIONS
# ---------
parser = OptionParser()

# --- Definition ---
# fichier de configuration
parser.add_option("--file",dest="file",metavar="FILE",
                  help="Name of the file that defines the simulations to perform")

# chemin de sortie
parser.add_option("--pathout",dest="pathout",metavar="PATHOUT",                       
                  help="Name of output directory (current by default)")

# nom executable
parser.add_option("--exe",dest="exe",metavar="exe",
                  help="Name of the ECH2O exec file")

# nom executable
parser.add_option("--cfg",dest="cfg",metavar="cfg",
                  help="Name of the ECH2O config file")

# nombre d'annees de simulation
# parser.add_option("--nyears",dest="nyears",metavar="nyears",
#                  help="Number of years considered")

(options, args) = parser.parse_args()

# --- Gestion ---
PATH_MAIN = os.getcwd()+'/'

print '*********************************************************************************'
print 'MORRIS SENSITIVITY TEST : '
print '- construction of the trajectories'
print '- forward runs'
print '- storage of outputs and info for posterior analysis : elementary effects, etc.' 
print '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
print

if options.file != None:
    [file,ext]=os.path.splitext(options.file)
    if len(glob.glob(file+'.csv')) ==0 : sys.exit('# STOP. The file that define the assimilation characteristics does not exist : \n   ...'+file+'.csv')
    file = options.file
    frep = os.path.dirname(file)   
    if frep == '':  file = os.path.join(PATH_MAIN,file)

if options.pathout == None:
    options.pathout = os.getcwd()+'/'+os.path.splitext(options.file)[0]

PATH_OUT  = os.path.abspath(options.pathout)
PATH_EXEC = os.path.abspath(os.path.join(PATH_OUT,'tmp'))
PATH_SPA = os.path.abspath(os.path.join(PATH_OUT,'Spatial'))
PATH_SPA_REF = os.path.abspath(os.path.join(PATH_MAIN,'Spatial')) # where the base maps are first copied from
PATH_CLIM = os.path.abspath(os.path.join(PATH_MAIN,'..','Climate'))

print 'Output:', PATH_OUT
print 'Maps & parameters:', PATH_SPA
print 'Climatic forcing:', PATH_CLIM

# -----------------------------------------------
# --- DEFINITION DES CONSTANTES DE SIMULATION ---
# -----------------------------------------------

### Catchment
# Number of soil unit
nsoil = 5
soils = np.asarray(['PG','PP','BR','PR','P'])
# Vegetation types
nveg = 4
vegs = np.asarray(['Pine','Hthr','Sphagnum','Molinia'])

### Simulation data
# ---------------------

# --- Execute command definition ---

cdumpsize = 'ulimit -c 33000'      #'limit coredumpsize 33000'
stacksize = 'ulimit -s unlimited ' #'limit stacksize unlimited'

if options.exe != None:
    if len(glob.glob(options.exe))==0:
        sys.exit('The user provided EXEC file was not found: '+options.exe)
    exe_ech2o = options.exe
    print
    print 'The user provided EXEC file is: '+options.exe

if options.cfg != None:
    if len(glob.glob(options.cfg))==0:
        sys.exit('The user provided CFG file was not found: '+options.cfg)
    cfg_ech2o = options.cfg
    print
    print 'The user provided CFG file is: '+options.cfg

#cmde_ech2o = ' ; '.join([cdumpsize,stacksize,'./'+exe_ech2o])
cmde_ech2o = ' '.join(['./'+exe_ech2o,cfg_ech2o])

# -- Creation of output directory
if len(glob.glob(PATH_OUT)) == 0: os.system('mkdir '+ PATH_OUT)
if len(glob.glob(PATH_EXEC)) == 0: os.system('mkdir '+ PATH_EXEC)

print
print 'Runs: ', PATH_EXEC
print 'Outputs: ', PATH_OUT
#print 'PATH_MAIN : ', PATH_MAIN

# -- Copy of executable and cfg's symbolic link
if len(glob.glob(os.path.join(PATH_EXEC,exe_ech2o))) != 0:
    os.system('rm '+os.path.join(PATH_EXEC,exe_ech2o))
os.symlink(os.path.join(PATH_MAIN, exe_ech2o) , os.path.join(PATH_EXEC,exe_ech2o) )
#if len(glob.glob(os.path.join(PATH_EXEC,cfg_ech2o))) != 0:
#    os.system('rm '+os.path.join(PATH_EXEC,cfg_ech2o))
os.system('cp '+os.path.join(PATH_MAIN, cfg_ech2o)+' '+os.path.join(PATH_EXEC,cfg_ech2o))
# Add a symbolic link to armadillo
libarma = 'libarmadillo.so.7'
if len(glob.glob(os.path.join(PATH_EXEC,libarma))) != 0:
    os.system('rm '+os.path.join(PATH_EXEC,libarma))
os.symlink(os.path.join(PATH_MAIN, libarma) , os.path.join(PATH_EXEC,libarma) )

# -- Creation of inputs directory
if len(glob.glob(PATH_SPA)) == 0: os.system('mkdir '+ PATH_SPA)
# Copy of reference data
os.system('cp -p '+PATH_SPA_REF+'/* '+PATH_SPA)

### Morris data
# -------------------
# Number of levels by parameter
nlevel = 7
# Number of trajectories
nr = 7
# Read parameters names, range, etc. in configuration file (.csv)
vardef = {}
with open(file, 'rb') as csvfile:
    paramread = list(csv.reader(csvfile, delimiter=','))
    vardef['name']  = np.asarray(paramread[0][1:len(paramread[1])],dtype='|S20')
    vardef['prior'] = np.asarray(paramread[1][1:len(paramread[1])],dtype=np.float32)
    vardef['min']   = np.asarray(paramread[2][1:len(paramread[1])],dtype=np.float32)
    vardef['max']   = np.asarray(paramread[3][1:len(paramread[1])],dtype=np.float32)
    vardef['soil']  = np.asarray(paramread[4][1:len(paramread[1])],dtype=np.int64)
    vardef['veg']   = np.asarray(paramread[5][1:len(paramread[1])],dtype=np.int64)
    vardef['log']   = np.asarray(paramread[6][1:len(paramread[1])],dtype=np.int64)
    vardef['inv']   = np.asarray(paramread[7][1:len(paramread[1])],dtype=np.int64)
    vardef['file']   = np.asarray(paramread[8][1:len(paramread[1])],dtype='|S20')    
exit

nparadef = len(vardef['name'])

# ---------------------------------------
# --- "Spatialize" the parameters     ---
# ---------------------------------------

repets = vardef['soil']*(nsoil-1)+vardef['veg']*(nveg-1)+1
varnames = np.repeat(vardef['name'], repets, axis=0)
mins = np.repeat(vardef['min'], repets, axis=0)
maxs = np.repeat(vardef['max'], repets, axis=0)
priors = np.repeat(vardef['prior'], repets, axis=0)
vlog = np.repeat(vardef['log'], repets, axis=0)
vinv = np.repeat(vardef['inv'], repets, axis=0)
vfile = np.repeat(vardef['file'], repets, axis=0)

varsoil = np.repeat(vardef['soil'], repets, axis=0)
varveg = np.repeat(vardef['veg'], repets, axis=0)
npara = len(varnames)

# Rename params
i=0
tmp = copy.copy(varnames)
varnames2 = copy.copy(varnames)
while i < npara:
    if varsoil[i] == 1:
        for j in range(nsoil):
            varnames[i+j]='_'.join([tmp[i+j],soils[j]])
            varnames2[i+j]='_'.join([tmp[i+j],str(j+1)])
        i=i+nsoil
    elif varveg[i] == 1:
        for j in range(nveg):
            varnames[i+j]='_'.join([tmp[i+j],vegs[j]])
            varnames2[i+j]='_'.join([tmp[i+j],str(j+1)])
        i=i+nveg
    else :
        i=i+1

# Total number of runs
nruns = (npara+1) * nr
print
print 'REMINDER --> Total number of runs : ', str(nruns)
#sys.exit()

# -------------------------------------------
# Output(s) used for the elementary effects -
# -------------------------------------------

obsname = ['Streamflow','Evap','SoilMoistureL1','SoilMoistureL2','SoilMoistureL3','SoilMoistureAv']
nobs = len(obsname)

# -----------------------------------------------------
# --- Creation of the Morris trajectories ---
# -----------------------------------------------------

varvals={}

# Value possible for each parameter
step = np.zeros((npara),np.float64)
i=0
for pname in varnames:
    if vlog[i]==0: # Linear
        step[i] = (maxs[i]-mins[i])/(nlevel-1)
        varvals[pname] = np.arange(mins[i],maxs[i], step[i], np.float64)
    if vlog[i]==1: # Log
        step[i] = (np.log(maxs[i])-np.log(mins[i]))/(nlevel-1)
        varvals[pname] = np.exp(np.arange(np.log(mins[i]),np.log(maxs[i]), step[i], np.float64))
    i=i+1
    
# Construct B* matrix, for each repetition
Bstar = np.zeros((npara,npara+1,nr),np.float64)
Xstar = np.zeros((npara),np.float64)
ind_hist = np.zeros((npara,nr)) #history of the index of trajectory construction

for ir in range(nr):

    # Generate 'base' vector X*
    ind_tmp = np.zeros((npara))
    for i in range(npara):
        ind_tmp[i] = random.randint(0,nlevel-1)
        if vlog[i]==0: # Linear variation
            Xstar[i]=mins[i]+ind_tmp[i]*step[i]
            if Xstar[i] < mins[i]-np.abs(mins[i])*0.000001:
                print varnames[i],Xstar[i], mins[i],ind_tmp[i],step[i]#,mins[i]+ind_tmp[i]*step[i]
                sys.exit()
        else: # Log variation
            Xstar[i]=np.exp(np.log(mins[i])+ind_tmp[i]*step[i])
            if Xstar[i] < mins[i]-np.abs(mins[i])*0.000001:
                print varnames[i],Xstar[i], mins[i],ind_tmp[i],step[i]#,mins[i]+ind_tmp[i]*step[i]
                sys.exit()

    # Generate trajectory
    # first: the first point comes from an increase by some steps of a random number of components
    Bstar[:,0,ir]=Xstar
    #if ir==0:
    # ind1 = range(npara)
    # random.shuffle(ind1)

    # #maxind = random.randint(1,npara)
    # gap_tmp = N.zeros((npara),savespace=1)
    # for i in range(maxind):
    #     seq = range(ind_tmp[ind1[i]],nniveaux)
    #     gap_tmp[ind1[i]] = random.choice(seq)
    #     Bstar[ind1[i],0,ir] = mins[ind1[i]]+gap_tmp[ind1[i]]*step[ind1[i]]
    #     if Bstar[ind1[i],0,ir] < Xstar[ind1[i]]: 
    #         #print mins[ind1[i]], maxs[ind1[i]],step[ind1[i]],ind1_tmp[ind1[i]]
    #         #print Bstar[ind1[i],0,ir], Xstar[ind1[i]]
    #         sys.exit('Error in the definition of X(1) !!')
    #print (Bstar[:,0,ir]-Xstar)/step
    # then: all the next are changed by only one component with + or - one step
    #if ir==0:
    ind2 = range(npara)
    random.shuffle(ind2)
    ind_hist[:,ir]=ind2
    for i in range(npara):
        #print varnames[ind2[i]]
        Bstar[:,i+1,ir] = Bstar[:,i,ir] #copy previous location
        eps = random.choice([-1,1])

        if vlog[ind2[i]]==0: #linear
            tmp = Bstar[ind2[i],i,ir] + step[ind2[i]]*eps #add (or substract) one step for a random component (not used)
            if ind_tmp[ind2[i]] == 0: # case where we already are on the lower boundary : add
                tmp = Bstar[ind2[i],i,ir] + step[ind2[i]]
            if ind_tmp[ind2[i]] == nlevel-1 : # case where we already are on the upper boundary : sbtrkt
                tmp = Bstar[ind2[i],i,ir] - step[ind2[i]]

        if vlog[ind2[i]]==1: #log
            tmp = np.exp(np.log(Bstar[ind2[i],i,ir]) + step[ind2[i]]*eps) #add (or substract) one step for a random component (not used)
            if ind_tmp[ind2[i]] == 0: # case where we already are on the lower boundary : add
                tmp = np.exp(np.log(Bstar[ind2[i],i,ir]) + step[ind2[i]])
            if ind_tmp[ind2[i]] == nlevel-1 : # case where we already are on the upper boundary : sbtrkt
                tmp = np.exp(np.log(Bstar[ind2[i],i,ir]) - step[ind2[i]])

                
        Bstar[ind2[i],i+1,ir] = tmp
        if tmp>maxs[ind2[i]]+np.abs(maxs[ind2[i]])*0.000001 or tmp<mins[ind2[i]]-np.abs(mins[ind2[i]])*0.000001:
            print 'Error in the incrementation of the parameter '+varnames[ind2[i]]
            print Xstar[ind2[i]],Bstar[ind2[i],i,ir],Bstar[ind2[i],i+1,ir],tmp
            print mins[ind2[i]], maxs[ind2[i]],step[ind2[i]],ind2_tmp[ind2[i]],gap_tmp[ind2[i]]
            sys.exit()

    #print Bstar[ind2[i],i,ir], Bstar[ind2[i],i+1,ir]
    #print ind2
    #print ind_hist[:,ir]
        

# -------------------------------------
#- Write parameters values file

# -- Change directory
os.chdir(PATH_OUT)
#PATH_EXEC = os.path.abspath('')
#print
#print 'Current directory : ', PATH_OUT

# -- Write trajectories

# Each parameters in a file
for ip in range(npara):
    with open('trajectory_'+varnames[ip]+'.csv','wb') as csvfile:
        csv_writer = csv.writer(csvfile)
        for ir in range(nr):
            csv_writer.writerow(Bstar[ip,:,ir])
    exit

# Altogether (not possible for now with csv, will need to get netCDF...)
#Vars['Bstar'] = {'dim_name': ('npara','nsteps','nrepet',), \
#               'dim_size': (npara,npara+1,nr),\
#               'attr_name' : ['long_name'],\
#               'attr_value' : ['Matrix of Morris parameter trajectories'],\
#               'value' : Bstar}
# Write Bstar for each trajectory
for ir in range(nr):
    trajnb = '%02i' % int(ir+1)
    #print trajnb
    with open('Bstar_traj'+trajnb+'.csv','wb') as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(varnames)
        for irun in range(npara+1):
            csv_writer.writerow(Bstar[:,irun,ir])
    exit


# History of the component changed at each step of the trajectories
with open('trajectory_history.csv','wb') as csvfile:
    csv_writer = csv.writer(csvfile)
    for ir in range(nr):
        csv_writer.writerow(ind_hist[:,ir])
exit

# -----------------------------------------------
# --- Simulations when varying the parameters ---
# -----------------------------------------------

# A few functions...

pcrcalc = 'pcrcalc \''

# Writing a PCRaster map, soil-unit-dependent
def DefPCRMap(name,value,filename):
    if name in ['CHANWIDTH','MANNING']:
        cmde = pcrcalc+filename+'.map = chanmask.map * '+str(value)+'\''
    elif name == 'CHANSEEP':
        cmde = pcrcalc+filename+'.map = chanmask.map/chanmask.map * '+str(value)+'\''
    elif name in ['TDAMP','TSOILI','TDAMP','DAMPDEPTH']:
        cmde = pcrcalc+filename+'.map = unit.map * '+str(value)+'\''
    else:
        print 'Error: this case was not foreseen !'
        sys.exit()
    return cmde


# Remove the default map for the tested parameters (just in case)
for ip in range(nparadef):
    if vardef['veg'][ip]==0:
        os.system('rm -f '+os.path.join(PATH_SPA,vardef['file'][ip]+'.map'))

# Each trajectory
for ir in range(nr):
    trajnb = '%02i' % int(ir+1)
    print 
    print '======================================'
    print '## Trajectory '+trajnb+' out of '+str('%02i' % int(nr))
    print '----------------------------------'

    vrun = {}
    
    # There are npara+1 runs for each trajectory
    for irun in range(npara+1):

        # -- Prepare the runs : constructs maps and/or parameter file with
        #                       initialized/updated parameters
        os.chdir(PATH_SPA)
                
        runnb = '%02i' % int(irun+1)
        print 'Run '+str(irun+1)+' out of '+str('%02i' % int(npara+1))
        if irun == 0:

            print
            print '|- Initializing parameter maps / table for the runs...'
            print

            # Assign the first vector parameters values
            for ip in range(npara):
                vrun[varnames[ip]] = Bstar[ip,irun,ir]
                #print '     '+varnames[ip]+' : '+str(vrun[varnames[ip]])
                if vinv[ip]==1:
                    vrun[varnames[ip]] = 1/vrun[varnames[ip]]
                    
                # Create the input file for this parameter
                # 1. if soil-dependent
                if varsoil[ip]==1:
                    cmd = pcrcalc + vfile[ip]+'.map = '
                    isl = int(varnames2[ip].split('_')[1])
                    if(isl == 1):
                        cmde = cmd + 'unit.soil1.map *'+str(vrun[varnames[ip]])+'\''
                        print varnames[ip]+' :',cmde
                        os.system(cmde)
                    else:
                        if vfile[ip]!=vfile[ip-1] or isl>nsoil:
                           print 'Error: you messed up with the soil dependent asignation !!'
                           sys.exit()
                        cmde = cmd+vfile[ip]+'.map + unit.soil'+str(isl)+'.map *'+str(vrun[varnames[ip]])+'\''
                        print varnames[ip]+' :',cmde
                        os.system(cmde)

                # 2. if vegetation dependent : TBD
                elif varveg[ip]==1:
                    print 'haha veg!'
                # 3. If not dependent
                else:
                    cmde = DefPCRMap(varnames[ip],vrun[varnames[ip]],vfile[ip])
                    print varnames[ip]+' :',cmde
                    os.system(cmde)

                # Update parameters that depend on others !
                #    here assuming that all roots are within the first two hydraulics layers,
                #    so that fracroot2 = 1 - fracroot1
                if vfile[ip]=='rootfrac1':
                    cmde = pcrcalc+'rootfrac2.map = unit.map - rootfrac1.map\''
                    print 'Update second layer roots :',cmde
                    os.system(cmde)

                # 
        if irun >= 1:

            # Change only the relevant parameter
            print
            print '|- Updating parameter maps / table within the trajectory...'

            ind_new = ind_hist[irun-1,ir]
            vold = vrun[varnames[ind_new]]                   # Back up old value
            vrun[varnames[ind_new]] = Bstar[ind_new,irun,ir] # Assign new value
            print '     Modif of '+varnames[ind_new]+' : '+str(Bstar[ind_new,irun-1,ir])+' --> '+str(Bstar[ind_new,irun,ir])
            if Bstar[ind_new,irun,ir] == Bstar[ind_new,irun-1,ir]: 
                sys.exit('ERROR : the wrong parameter was assumed as changed !!!')
        
            if vinv[ind_new]==1:
                print 'REMINDER : this parameter value is inverted for the simulation, so --'
                vrun[varnames[ind_new]] = 1/vrun[varnames[ind_new]]
                    
            # Create the input file for this parameter
            # 1. if soil-dependent
            if varsoil[ind_new]==1:
                cmd = pcrcalc + vfile[ind_new]+'.map = '+vfile[ind_new]+'.map + unit.soil' 
                isl = int(varnames2[ind_new].split('_')[1])
                # "update" vfile, new - old
                cmde = cmd+str(isl)+'.map * ('+str(vrun[varnames[ind_new]])+'-'+str(vold)+')\''
                print cmde
                os.system(cmde)
            # 2. if vegetation dependent : TBD
            elif varveg[ind_new]==1:
                print 'haha veg!'
            # 3. If not dependent
            else:
                cmde = DefPCRMap(varnames[ind_new],vrun[varnames[ind_new]],vfile[ind_new])
                print cmde
                os.system(cmde)

            # Update parameters that depend on others !
            #    here assuming that all roots are within the first two hydraulics layers,
            #    so that fracroot2 = 1 - fracroot1
            if vfile[ind_new]=='rootfrac1':
                cmde = pcrcalc+'rootfrac2.map = unit.map - rootfrac1.map\''
                print 'Update second layer roots :',cmde
                os.system(cmde)
                    
        # -- Update initial conditions that depend on parameterization !
        #    e.g., initial soil moisture (typically depends on porosity)
        if irun==0 or (irun>=1 and vfile[ind_new]=='poros'):
            print 'Update initial conditions !!'
            cmde = pcrcalc+'SWC1.map = poros.map * 0.8\''
            os.system(cmde)
            cmde = pcrcalc+'SWC2.map = poros.map * 0.8\''
            os.system(cmde)
            cmde = pcrcalc+'SWC3.map = poros.map * 0.8\''
            os.system(cmde)
            
        # -- Move to run directory and create links
        os.chdir(PATH_EXEC)
        
        # -- Execute ECH2O --
        print
        print 'ECH2O simulation in progress...'
        #os.system(cmde_orchidee + ' > orchidee.log')
        runall = '%03i' % int(nr*ir+irun+1)
        #print cmde_ech2o
        #os.system(cmde_ech2o)
        os.system('time ' + cmde_ech2o + ' > ech2o.'+runall+'.log')
        print
        
        # -- Outputs reorganization
        # Basin summary
        if len(glob.glob(os.path.join(PATH_EXEC,'BasinSummary.txt'))) == 0:
            print 'Error in the execution of ECH2O...'
            sys.exit()
        os.system('mv -f BasinSummary.txt BasinSummary_traj'+trajnb+'_run'+runnb+'.txt')      
        
        # Observation used for sensitivity
        for io in range(nobs):
            obs_out = os.path.join(PATH_OUT,obsname[io]+'_traj'+trajnb+'_run'+runnb+'.tab')
            print 'Moving sensitivity outputs to : ',obs_out
            os.system('mv -f '+obsname[io]+'.tab '+obs_out)
        print
    #sys.exit()

# END PROGRAM MAIN

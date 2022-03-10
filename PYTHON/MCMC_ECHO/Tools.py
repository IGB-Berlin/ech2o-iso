#!/usr/bin/env python
# -*- coding: utf-8 -*-
# *************************************************
#
# MCMC calibration algorithm for ECH2O
#
# -------
# Routine: Subroutines for parameters & obs manip & outputs
# -------
# Author: S. Kuppel
# Created on 11/2016
# -------------------------------------------------

from pcraster import *
import time, os, glob, sys, copy
import csv
from time import localtime, strftime
import random
import numpy as np
from datetime import datetime, timedelta

import eval_measures as metrics

# ----------------------------------------------------------------------------
# -- Initialize parameters and optimization variables

def InitOpti(Opti, Paras, Site):

    # -- Parameters: from definition to all values
    Paras.names = Paras.ref.keys()
    Paras.names.sort()
    Paras.n = len(Paras.names)
    # Read dictionary to get all param setup
    Opti.xprior = [] # First guess
    Opti.sig = []    # Prior standard deviation
    Opti.min = []    # Lower boundary of sampling range
    Opti.max = []    # Upper boundary of sampling range
    Opti.log = []    # Logarithmic sampling ? (0=no, 1=yes)
    Opti.names = []  # Upper boundary of sampling range
    Opti.ind = []    # Parameter number (common for all its components)
    Paras.ind = {}   # Index of each parameter components in Opti vectors
    ipar=0
    ipar2=0
    Paras.isveg = 0
    for par in Paras.names:
        # Dimensions (soil or veg or 1)
        nr = Paras.ref[par]['soil']*(Site.ns-1)+Paras.ref[par]['veg']*(Site.nv-1)+1
        # Build vectors used in the optimisation
        if type(Paras.ref[par]['min'])== float:
            Opti.xprior += [Paras.ref[par]['xfg']]
            Opti.min += [Paras.ref[par]['min']]
            Opti.max += [Paras.ref[par]['max']]
        else:
            Opti.xprior += Paras.ref[par]['xfg']
            Opti.min += Paras.ref[par]['min']
            Opti.max += Paras.ref[par]['max']

        Opti.log += list(np.repeat(Paras.ref[par]['log'],nr))
        # Link betwen params and all variables
        Opti.ind += list(np.repeat(ipar,nr))
        ipar+=1
        if nr>1: Paras.ind[par] = list(np.arange(ipar2,ipar2+nr,1))
        if nr==1: Paras.ind[par] = ipar2
        ipar2+=nr
        # For outputs
        if Paras.ref[par]['soil']==1:
            Opti.names += [par + '_' + s for s in Site.soils]
        elif Paras.ref[par]['veg']==1:
            Opti.names += [par + '_' + s for s in Site.vegs]
            Paras.isveg += 1
        else:
            Opti.names += [par]
    
    # Total number of variables
    Opti.nvar = len(Opti.min)

    # Prior standard deviation (used for random-walk)
    Opti.sig = copy.copy(Opti.xprior)
    for i in range(Opti.nvar):
        # If log sampling, the standard deviation will in log unit (see ParasWalk) !
        if Opti.log[i]==1:
            Opti.sig[i] = (10**Opti.max[i]-10**Opti.min[i])*Opti.sigp
        else:
            Opti.sig[i] = (Opti.max[i]-Opti.min[i])*Opti.sigp


# ----------------------------------------------------------------------------
# -- Initialize maps and other un inputs

def InitInputs(Opti, Paras, Site, Config):

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
    # "Header": number of species and of params
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

# ----------------------------------------------------------------------------
# -- Generate random set of parameters values of write it

def ParasWalk(Opti, it):

    if it==1:
        Opti.nsamp = 0
        Opti.x = copy.copy(Opti.xprior)

    # Dummy copy just to get the shape of Opti.xnew
    Opti.xnew = copy.copy(Opti.x)

    #-- Normally-distributed value centered on current sample
    for i in range(Opti.nvar):
        # Log transform where needed
        if Opti.log[i]==1:
            Opti.xnew[i] = 10**(np.random.normal(np.log10(Opti.x[i]),Opti.sig[i]))
            # Check boundaries...if outside, resample
            while Opti.xnew[i]>Opti.max[i] or Opti.xnew[i]<Opti.min[i]:
                Opti.xnew[i] = 10**(np.random.normal(np.log10(Opti.x[i]),Opti.sig[i]))
        else:
            Opti.xnew[i] = np.random.normal(Opti.x[i],Opti.sig[i])
            # Check boundaries...if outside, resample
            while Opti.xnew[i]>Opti.max[i] or Opti.xnew[i]<Opti.min[i]:
                Opti.xnew[i] = np.random.normal(Opti.x[i],Opti.sig[i])

    # -- Total number of sampling (to compute acceptance rate)
    Opti.nsamp += 1

# ----------------------------------------------------------------------------
# -- Updating inputs for ECH2O

def CreateInputs(Opti, Paras, Site, Config, it):

    # Small switch not to read the vegetation params file every time
    readveg=1
    
    # -- Get the parameter sample prior or current iteration
    if it==0:
        Opti.x = copy.copy(Opti.xprior)
    else:
        Opti.x = copy.copy(Opti.xnew)

    # --
    for pname in Paras.names:

        #print pname
        
        ## - Mapped parameters
        if Paras.ref[pname]['veg']==0:

            # Soil unit dependence
            if Paras.ref[pname]['soil']==1 and pname!='Rootp.1+2':
                #print 'Soil dependent !!'
                outmap = Site.bmaps['unit']*0
                # Read each soil map unit and apply param value
                for im in range(Site.ns):
                    outmap+= Site.bmaps[Site.soils[im]]*Opti.x[Paras.ind[pname][im]]

            # No spatial/veg dependence, but channel stuff
            elif pname!='Rootp.1+2':
                #print 'Not dependent !!'
                if Paras.ref[pname]['file'] in ['chanwidth','chanmanningn']:
                    outmap = Site.bmaps['chanmask']*Opti.x[Paras.ind[pname]]
                elif Paras.ref[pname]['file'] == 'chanparam':
                    outmap = Site.bmaps['chanmask_NaN']*Opti.x[Paras.ind[pname]]
                else:
                    outmap = Site.bmaps['unit']*Opti.x[Paras.ind[pname]]
            
            # Outputs updated map
            report(outmap,Config.PATH_SPA+'/'+Paras.ref[pname]['file']+'.map')
        
            # Check for initial condition/other parameter dependence
            # Initial soil water content
            if pname=='Porosity':
                report(outmap*0.8,Config.PATH_SPA+'/SWC1.map')
                report(outmap*0.8,Config.PATH_SPA+'/SWC2.map')
                report(outmap*0.8,Config.PATH_SPA+'/SWC3.map')
            # Root layer 2: summed fraction over 1+2 minus fraction over 1
            if pname=='Rootp.1':
                outmap2 = Site.bmaps['unit']*0
                for im in range(Site.ns):
                    outmap2+= Site.bmaps[Site.soils[im]]*Opti.x[Paras.ind['Rootp.1+2'][im]]
                #base = readmap(Config.PATH_SPA+'/unit.map')
                report(outmap2 - outmap, Config.PATH_SPA+'/rootfrac2.map')
        
        ## - Vegetation parameters
        else:
            # Change the value based on param name correspondance
            vegnew = copy.copy(Opti.vref)
            #print Opti.vref
            for iv in range(Site.nv):
                vegnew[iv][vegnew['name'].index(pname)] = str(Opti.x[Paras.ind[pname][iv]])

    # ------------------------------------------------------------------------------
    # Finalizing the preps....

        
    ## - Finalizing soil parameterization
    # Check that initial soil moisture is not smaller residual soil 
    # tbd... for now just pick the porosity and thetar reange wisely enough

    ## - Finalizing the vegetation parameterization
    if Paras.isveg > 0:
    # Equalize leaf turnover and additional turnover due to water and/or temperature stress
        #for iv in range(Site.nv):
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


# ----------------------------------------------------------------------------
# -- Run ECH2O

def runECH2O(Config):

    os.chdir(Config.PATH_EXEC)
    print '--> running ECH2O on '+Config.ncpu+' threads...'
    start = time.time()
    os.system(Config.cmde_ech2o+' > ech2o.log')
    print '    run time:',time.time() - start,'seconds (limit at '+Config.tlimit+')'

# ----------------------------------------------------------------------------
# -- Check is ECH2O ran properly

def runOK(Data, Config):

    isOK = 0

    # 1. Check it ran
    f_test1 = Config.PATH_EXEC+'/BasinSummary.txt'
    #f_test2 = Config.PATH_EXEC+'/'+Data.obs[Data.names[0]]['sim_file']
    f_test2 = Config.PATH_EXEC+'/'+Data.obs['Streamflow']['sim_file']
    if len(glob.glob(f_test1))*len(glob.glob(f_test2))!=0:
        #print 'Outputs are there...'

        # 2. Check it ran until the end
        #tmp1 = np.genfromtxt(f_test1,skip_header=1,unpack=True)[0]
        tmp2 = np.genfromtxt(f_test2,skip_header=Data.nts+3,unpack=True)[0]
        if len(tmp2)==Data.lsim: #& len(tmp2)==Data.lsim:
            #print 'It ran until the end !'
            isOK = 1

    return isOK

# ----------------------------------------------------------------------------
# -- Cost function, with a chosen metric

def CostFunction(Opti, Data, Config, it):

    if it==0:
        Opti.sigpr = {}
        Opti.Jobs = {}
        Opti.Jprior = 0
    else:
        Opti.Jnew = 0

    ### For each obs -------------------------------------------------------------------------
    for oname in Data.names:
        #print oname
        #print '           for '+obsnames[iobs]+': ',

        # -- Get the obs
        f_obs = Data.obs[oname]['obs_file']
        tmp = np.genfromtxt(f_obs,delimiter=',',skip_header=1,unpack=True, dtype= '|S10')[0]
        lobs = len(tmp)
        obst = [datetime.strptime(a, '%d/%m/%Y') for a in tmp]
        obs = np.genfromtxt(f_obs,delimiter=',',skip_header=1,unpack=True)[Data.obs[oname]['obs_col']]
        # -- Get simulation outputs
        f_sim = Data.obs[oname]['sim_file']

        # Get the right index (because ECH2O messes up point order...) using Data.sim_order
        # ...not forgetting that point 1 is index 0 (thus -1) and that the first column is time (thus +1)
        # careful with Saturation it is BasinSummary !!
        if oname == 'SaturationArea':
            hskip = 1
            idx = Data.obs[oname]['sim_pts']-1
        else:
            idx = np.argsort(np.array(Data.sim_order))[Data.obs[oname]['sim_pts']-1]+1
            hskip = Data.nts+3

        #print lobs, Data.lsim, oname, idx,
        # ...and apply conversion factor (for Transpiration mostly)
        sim = np.genfromtxt(f_sim,delimiter='\t',skip_header=hskip,unpack=True)[idx]*Data.obs[oname]['conv']

        # Crop between desired common time frame (obs is limiting)
        obs = [obs[idx] for idx in range(lobs) if obst[idx]>=Data.obs[oname]['fit_beg'] and obst[idx]<=Data.obs[oname]['fit_end']]
        sim = [sim[idx] for idx in range(Data.lsim) if Data.simt[idx]>=Data.obs[oname]['fit_beg'] and Data.simt[idx]<=Data.obs[oname]['fit_end']]
        lfit = len(sim)
        # Check the time frames
        if len(obs)!=len(sim):
            sys.exit('Mismatch between sim et obs time frames! (sim='+str(len(sim))+', obs='+str(len(obs))+')')

        # -- Cost function for this observation
        # 
        if it==0:
            # Prior Opti.Jobspr is the prior standard deviation in observation space,
            # i.e. here the L2 misfit divdided by the number of date points
            # thus assuming a diagonal error covariance matrix without cross-correlations...
            Opti.sigpr[oname] = metrics.mse(sim,obs)
            # Prior cost function = L2 misfit divided by standard deviation = number of date points
            Opti.Jobs[oname] = lfit
            # Total cost function
            Opti.Jprior += Opti.Jobs[oname]
        else:
            # L2 misfit divided by standard deviation
            Opti.Jobs[oname] = lfit*metrics.mse(sim,obs) / Opti.sigpr[oname]
            # Total cost function
            Opti.Jnew += Opti.Jobs[oname]

# ----------------------------------------------------------------------------
# -- Accept or reject new sample based on cost function

def CriterionWalk(Opti, it):

    # Initialization
    Opti.walkOK = 0
    if it == 1 :
        Opti.walks = 0
        Opti.walkhist = []
        if Opti.Jnew < Opti.Jprior:
            Opti.x = copy.copy(Opti.xnew)
            Opti.J = copy.copy(Opti.Jnew)
            # Acceptance occurence
            Opti.walkOK = 1
            Opti.walks += 1
            Opti.walkhist += [1.]
            # Some output..
            print '    Random-walk accepted...moving on !'
        else:
            Opti.J = copy.copy(Opti.Jprior)
            # Acceptance occurence
            Opti.walkhist += [0.]
            print '    Random-walk rejected...try again !'
        
    # Compare old and new cost functions
    elif Opti.Jnew < Opti.J:
        Opti.x = copy.copy(Opti.xnew)
        Opti.J = copy.copy(Opti.Jnew)
        # Acceptance occurence
        Opti.walkOK = 1
        Opti.walks += 1
        Opti.walkhist += [1.]
        # Some output..
        print '    Random-walk accepted...moving on !'
    else:
        Opti.walkhist += [0.]
        print '    Random-walk rejected...try again !'

        
# ----------------------------------------------------------------------------
# -- Write parameters historic

def ParasOutputs(Opti,Config, it):
  
    # Open one file for all samples
    if it==0:
        # Whole historic
        Config.f_par = Config.PATH_OUT+'/Parameters_hist.txt'
        with open(Config.f_par,'w') as f_in:
            f_in.write('Iteration,'+','.join(Opti.names)+'\n')
            f_in.write(str(it)+','+','.join([str(x) for x in Opti.xprior])+'\n')
        # Convergence historic
        Config.f_parOK = Config.PATH_OUT+'/Parameters_histconv.txt'
        with open(Config.f_parOK,'w') as f_in:
            f_in.write('Iteration,'+','.join(Opti.names)+'\n')
            f_in.write(str(it)+','+','.join([str(x) for x in Opti.xprior])+'\n')

    else:            
        # Whole historic
        with open(Config.f_par,'a') as f_in:
            f_in.write(str(it)+','+','.join([str(x) for x in Opti.x])+'\n')
        # Convergence historic
        if it==0 or Opti.walkOK == 1:
            with open(Config.f_parOK,'a') as f_in:
                f_in.write(str(it)+','+','.join([str(x) for x in Opti.x])+'\n')

# ----------------------------------------------------------------------------
# -- Write parameters historic

def OptiOutputs(Opti, Data, Config, it):
  
    # Open one file for all samples
    if it==0:
        # Whole historic
        Opti.f_hist = Config.PATH_OUT+'/Opti_hist.txt'
        with open(Opti.f_hist,'w') as f_in:
            tmp = ['J_'+x for x in Opti.names]
            f_in.write('Iteration,'+','.join(tmp)+',Jtotal,OKRate,OKrate_last20\n')
            f_in.write(str(it)+','+','.join([str(Opti.Jobs[x]) for x in Data.names])+','+
                       str(Opti.Jprior)+',NA,NA\n')
        #print strout
        # Convergence historic
        Opti.f_histOK = Config.PATH_OUT+'/Opti_histconv.txt'
        with open(Opti.f_histOK,'w') as f_in:
            tmp = ['J_'+x for x in Opti.names]
            f_in.write('Iteration,'+','.join(tmp)+',Jtotal,OKRate,OKRate_last20\n')
            f_in.write(str(it)+','+','.join([str(Opti.Jobs[x]) for x in Data.names])+','+
                       str(Opti.Jprior)+',NA,NA\n')

    else:            
        if len(Opti.walkhist)>=20:
            rate20 = str(100*sum(Opti.walkhist[::-1][i]) for i in range(20))/20.)
        else:
            rate20 = 'NA'

        with open(Opti.f_hist,'a') as f_in:
            f_in.write(str(it)+','+','.join([str(Opti.Jobs[x]) for x in Data.names])+','+
                       str(Opti.Jnew)+','+str(100.*Opti.walks/Opti.nsamp)+','+rate20+'\n')        
        # Convergence only
        if it==0 or Opti.walkOK == 1:
            with open(Opti.f_histOK,'a') as f_in:
                f_in.write(str(it)+','+','.join([str(Opti.Jobs[x]) for x in Data.names])+','+
                           str(Opti.J)+','+str(100.*Opti.walks/Opti.nsamp)+','+rate20+'\n')        

        
# ----------------------------------------------------------------------------
# -- Write simulation histoic

def ObsOutputs(Data, Opti, Config, it):

    # -- Group the output files in one across simulations, 
    #    separating by observations points and veg type where it applies

    # Initialization
    if it == 0:
        # Basin summary
        Config.f_bs = Config.PATH_EXEC+'/BasinSummary.txt'
        Opti.f_bs_hist = Config.PATH_OUT+'/BasinSummary_hist.txt'
        Opti.f_bs_histOK = Config.PATH_OUT+'/BasinSummary_histconv.txt'
        # Header of files
        with open(Config.f_bs,'r') as f_in:
            tmp = f_in.readline()
            with open(Opti.f_bs_hist,'w') as f_out:
                f_out.write('Iteration,'+','.join(tmp.split('\t')))
            with open(Opti.f_bs_histOK,'w') as f_out:
                f_out.write('Iteration,'+','.join(tmp.split('\t')))

        for oname in Data.names:
            # Historic time series file names
            Data.obs[oname]['sim_hist'] = Config.PATH_OUT+'/'+oname+'_hist.tab'
            Data.obs[oname]['sim_hist2'] = Config.PATH_OUT+'/'+oname+'_histconv.tab'
            # Header of files
            with open(Data.obs[oname]['sim_hist'],'w') as f_out:
                f_out.write('Iteration,'+','.join([str(i+1) for i in range(Data.lsim)])+'\n')
            with open(Data.obs[oname]['sim_hist2'],'w') as f_out:
                f_out.write('Iteration,'+','.join([str(i+1) for i in range(Data.lsim)])+'\n')

    # Save current run outputs (and delete the file to relieve Maxwell...)
    # Summary: last line
    tmp = np.genfromtxt(Config.f_bs,skip_header=Data.lsim)
    with open(Opti.f_bs_hist,'a') as f_out:
        f_out.write(str(it+1)+','+','.join([str(i+1) for i in list(tmp)])+'\n')
    if it==0 or Opti.walkOK == 1:
        with open(Opti.f_bs_histOK,'a') as f_out:
            f_out.write(str(it+1)+','+','.join([str(i+1) for i in list(tmp)])+'\n')
        
        
    # Time series
    for oname in Data.names:
        if oname == 'SaturationArea':
            hskip = 1
            idx = Data.obs[oname]['sim_pts']-1
        else:
            hskip = Data.nts+3
            idx = np.argsort(np.array(Data.sim_order))[Data.obs[oname]['sim_pts']-1]+1

        tmp = np.genfromtxt(Data.obs[oname]['sim_file'],delimiter='\t',
                            skip_header=hskip,unpack=True)[idx]

        with open(Data.obs[oname]['sim_hist'],'a') as f_out:
                f_out.write(str(it)+','+','.join([str(j) for j in list(tmp)])+'\n')
        if it==0 or Opti.walkOK == 1:
            with open(Data.obs[oname]['sim_hist2'],'a') as f_out:
                f_out.write(str(it)+','+','.join([str(j) for j in list(tmp)])+'\n')
            


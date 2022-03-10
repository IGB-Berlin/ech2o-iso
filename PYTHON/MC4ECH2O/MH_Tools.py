#!/usr/bin/env python
# -*- coding: utf-8 -*-
# *************************************************
#
# Monte Carlo calibration algorithm for ECH2O
#
# -------
# Routine: Subroutines for parameters & obs manip & outputs
# -------
# Author: S. Kuppel
# Created on 10/2016
# -------------------------------------------------

#from pcraster import *
import time, os, glob, sys, copy
import csv
from time import localtime, strftime
import random
import numpy as np

# ----------------------------------------------------------------------------
# -- Generate random set of parameters values

def gen_paras(Opti, it):

    #-- Random value for each parameters
    Opti.xnew = []
    for i in range(Opti.nvar):
        # Generate random-walk
        tmp = list(np.random.normal(Opti.xcurTr,Opti.sigT))
        # Log transform whre needed
        
        # log transform
        if Opti.log[i]==1:
            Opti.x.append(np.log10(np.random.uniform(10**Opti.xb[i],10**Opti.sigb[i])))
        else:
            Opti.x.append(np.random.uniform(Opti.min[i],Opti.max[i]))
        #print i, Opti.names[i], Opti.x[i]

    #-- Keep values taken over the iterations
    #if it==1:
    #    Opti.hist_x = np.zeros((Opti.nit,Opti.nvar))*np.nan
    #Opti.hist_x[it-1,:] = copy.copy(Opti.x)


# ----------------------------------------------------------------------------
# -- Write parameters values file

def output_par(Opti,Config,it):
  
    # Open one file for all samples
    if it==0:
        f_par = open(Config.PATH_OUT+'/parameters.txt','w')
        f_par.write(','.join(Opti.names)+'\n')
        f_par.close()

    f_par = open(Config.PATH_OUT+'/parameters.txt','a')
    f_par.write(','.join([str(x) for x in Opti.x])+'\n')
    #print Opti.names
    #print str(Opti.x)
    f_par.close()

    # All in the same file
    #with open(Config.PATH_OUT+'/VarHist.all.csv','w') as csvfile:
    #    csv_writer = csv.writer(csvfile)
    #    csv_writer.writerow(Opti.names)
    #    for ir in range(it+1):
    #        csv_writer.writerow(list(Opti.hist_x[it]))
    #exit

# ----------------------------------------------------------------------------
# -- Updating inputs for ECH2O

def create_inputs(Opti, Paras, Site, Config, it):

    # Small switch not to read the vegetation params file every time
    readveg=1

    for pname in Paras.names:

        ## - Mapped parameters
        if Paras.ref[pname]['veg']==0:

            # Soil unit dependence
            if Paras.ref[pname]['soil']==1 :
                #print 'Soil dependent !!'
                outmap = Site.bmaps['unit']*0
                # Read each soil map unit and apply param value
                for im in range(Site.ns):
                    outmap+= Site.bmaps[Site.soils[im]]*Opti.x[Paras.ind[pname][im]]

            # No spatial/veg dependence, but channel stuff
            else:
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
            if Paras.ref[pname]['file']=='poros':
                report(outmap*0.8,Config.PATH_SPA+'/SWC1.map')
                report(outmap*0.8,Config.PATH_SPA+'/SWC2.map')
                report(outmap*0.8,Config.PATH_SPA+'/SWC3.map')
            # Rootin layer 2: assuming that all roots are within the first two hydraulics layers
            if Paras.ref[pname]['file']=='rootfrac1':
                base = readmap(Config.PATH_SPA+'/unit.map')
                report(base - outmap, Config.PATH_SPA+'/rootfrac2.map')
        

        ## - Vegetation parameters
        else:
            # Change the value based on param name correspondance
            vegnew = copy.copy(Opti.vref)
            for iv in range(Site.nv):
                vegnew[iv][vegnew['name'].index(pname)] = str(Opti.x[Paras.ind[pname][iv]])

    # ------------------------------------------------------------------------------
    # Finalizing the preps....

        
    ## - Finalizing soil parameterization
    # Check that initial soil moisture is not smaller residual soil 
    # tbd... for now just pick the porosity and thetar reange wisely enough

    ## - Finalizing the vegetation parameterization
    if Paras.isveg > 0:
    # Equalize leaf turnover and additional turnover due to water and/or temperature stree
        for iv in range(Site.nv):
            vegnew[iv][vegnew['name'].index('TurnovL_MWS')] = copy.copy(vegnew[iv][vegnew['name'].index('TurnovL')])
            vegnew[iv][vegnew['name'].index('TurnovL_MTS')] = copy.copy(vegnew[iv][vegnew['name'].index('TurnovL')])
    # Write the vegetation params file (if there's at least one veg param)
        vegfile = open(Config.PATH_SPA+'/'+Site.vfile,'w')
        vegfile.write('\t'.join(Opti.vref['header'])+'\n')
        for iv in range(Site.nv):
            vegfile.write('\t'.join(vegnew[iv]) +'\n')
        vegfile.write('\t'.join(Opti.vref['footer'])+'\n')
        vegfile.write('\t'.join(Opti.vref['name'])+'\n')
        vegfile.close()
        
# ----------------------------------------------------------------------------
# -- Manage outputs        

def manage_outputs(Data, Opti, Config, it):

        # Check that the model ran OK...
        #f_bs = Config.PATH_EXEC+'/BasinSummary.txt'
        #if len(glob.glob(f_bs)) == 0:
        #    print 'Error in the execution of ECH2O...'
            #sys.exit()
        
        # -- Group the output files in one across simulations, 
        #    separating by observations points and veg type where it applies
        
        # Initialization
        if it == 0:
            # Grouped output file names
            Opti.fobs = {}
            for obs in Data.obs:
                # Streamflow: only one observation point
                if obs == 'Streamflow':
                    Opti.fobs[obs] = ['Streamflow.txt']
                else:
                    Opti.fobs[obs] = [obs+'.pt'+str(j+1)+'.tab' for j in range(Data.nts)]
                #elif Site.obsveg[i]==0:
                    #for j in range(Data.nts):
                        #Opti.fobs = Opti.fobs + [Site.obs[i]+'.pt'+str(j+1)+'.tab']
        #elif Site.obsveg[i]==1:       

            # Header of files
            # Summary
            #with open(Config.PATH_EXEC+'/BasinSummary.txt','r') as f_in:
            #    with open(Config.PATH_OUT+'/BasinSummary.txt','w') as f_out:
            #        f_out.write('Sample,'+','.join(f_in.readline().split('\t')))
            # Time series
            for obs in Data.obs:
                #print obs
                f_in = Config.PATH_EXEC+'/'+obs+'.tab'
                tmp = np.genfromtxt(f_in,skip_header=Data.nts+3,unpack=True)[0]
                if obs == Data.obs[0]:
                    # Read length of time series
                    Data.lts = len(tmp)
                for fobs in Opti.fobs[obs]:
                    #print fobs
                    with open(Config.PATH_OUT+'/'+fobs,'w') as f_out:
                        f_out.write('Sample,'+','.join([str(i+1) for i in range(Data.lts)])+'\n')

        # Save current run outputs (and delete the file to relieve Maxwell...)
        # Summary: last line
        #tmp = np.genfromtxt(Config.PATH_EXEC+'/BasinSummary.txt',skip_header=Data.nts)
        #with open(Config.PATH_OUT+'/BasinSummary.txt','a') as f_out:
        #    f_out.write(str(it+1)+','+','.join([str(i+1) for i in list(tmp)])+'\n')
        # Time series
        for obs in Data.obs:
            f_in = Config.PATH_EXEC+'/'+obs+'.tab'
            for i in range(len(Opti.fobs[obs])):
                tmp = np.genfromtxt(f_in,delimiter='\t',skip_header=Data.nts+3,unpack=True)[Data.order[i]]
                with open(Config.PATH_OUT+'/'+Opti.fobs[obs][i],'a') as f_out:
                    f_out.write(str(it+1)+','+','.join([str(j) for j in list(tmp)])+'\n')
        
            # Remove original time series
            #os.system('rm -f '+f_in)



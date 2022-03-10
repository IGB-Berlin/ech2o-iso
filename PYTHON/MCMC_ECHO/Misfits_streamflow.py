#    !/usr/bin/env python3
# -*- coding: utf-8 -*-
# *************************************************
#
# Model-data misfit after Monte-Carlo sampling +
# assemble parameters samples
#
# -------
# Author: S. Kuppel
# Created on 10/2016
# -------------------------------------------------

import numpy as np
from datetime import datetime, timedelta
import os, time, sys, glob, copy
import json

# --- Subroutine(s)
import eval_measures as metrics

# Switch : fetching obs metrics or corresponding parameters / best runs ?
# 0 = obs metrics, 1 = parameters / best runs
switch = 1
nbest = 200

# --- Definition ---
# Path to data
obsdir = '/users/s01ks6/Data/BB'
# Output subdirectories
MCname = 'MC1'
outdir = os.getcwd()+'/'+MCname+'_sampling.'
# Observations
obsnames = ['Streamflow']

# Number of samples
nit = 300
# There is a bug in the outputs (for now), so we skip a few runs...
rngit = range(nit)#range(27,nit)
# Number of parallel runs
njob = 96
# Number of locations for simulated time series
nts = 5

# Time frame frame used for the metrics
tb = datetime(2014,7,1)#datetime(2014,7,1)
te = datetime(2016,6,30)#datetime(2016,6,30)
lobs = 1827

#############################################################################################
# Observations
# ------------------------------

if switch == 0:

    print 'Getting the observation time series...'
    print

    obs={}

    # Streamflow
    file = obsdir+'/Q/BB_discharge_daily_01062011-04072016.csv'
    tmp = np.genfromtxt(file,delimiter=',',skip_header=1,unpack=True, dtype= '|S10')[0]
    date = [datetime.strptime(a, '%d/%m/%Y') for a in tmp]
    tmp = np.genfromtxt(file,delimiter=',',skip_header=1,unpack=True)[1]
    # Crop between desired time frame
    obs['Streamflow'] = [tmp[idx] for idx in range(len(tmp)) if date[idx]>=tb and date[idx]<=te]

#############################################################################################
# Simulations outputs : saving the metrics

    print 'Summarizing the metrics....'
    print
    iok = []
    jok = []
    # Streamflow
    # Time frame
    ltf = (te-tb).days+1
    # Open the 2 column (outlet pixel) in each output file
    with open(os.getcwd()+'/'+MCname+'_metrics_'+obsnames[0]+'.txt','w') as f_out:
        f_out.write('job,sample,NSE,KGE,RMSE,bias\n')
        for i in range(njob):
            print 'Job', i+1,'...'
            j = 0
            while j < nit:
                it = '%03i' % int(j+1)
                d_in = outdir+str(i+1)+'/'+it
                if len(glob.glob(d_in))==0:
                    break
                f_in = d_in+'/'+obsnames[0]+'.tab'
                if len(glob.glob(f_in))==0:
                    j+=1
                    continue
                tmp = np.genfromtxt(f_in,delimiter='\t',skip_header=nts+3,unpack=True)[1]
                #print i ,j ,len(tmp)
                if len(tmp) < lobs:
                    j+=1
                    continue
                # Crop between desired time frame (based on the length of this time frame)
                sim = [tmp[idx] for idx in range(len(tmp)-ltf,len(tmp))]
                # Metrics
                md_nse = metrics.nash_sutcliff(sim,obs['Streamflow'])
                md_kge = metrics.kling_gupta(sim,obs['Streamflow'],method='2012')
                md_rmse = metrics.rmse(sim,obs['Streamflow'])
                md_bias = metrics.bias(sim,obs['Streamflow'])
                # Write
                f_out.write(','.join([str(i+1),str(j+1),
                                      str(md_nse),str(md_kge),str(md_rmse),str(md_bias)])+'\n')
                # Save the 'coordinates' of this sample
                iok.append(i)
                jok.append(j)
                j+=1

    nok = len(iok)
    print 'Number of complete samples:',nok
    print

#############################################################################################
# Get parameters / best runs

if switch == 1:

    f_ijok = os.getcwd()+'/'+MCname+'_metrics_Streamflow.txt'
    iok = np.genfromtxt(f_ijok,delimiter=',',skip_header=1,unpack=True,dtype=np.int64)[0]-1
    jok = np.genfromtxt(f_ijok,delimiter=',',skip_header=1,unpack=True,dtype=np.int64)[1]-1
    fit_nse  = np.genfromtxt(f_ijok,delimiter=',',skip_header=1,unpack=True)[2]
    fit_kge  = np.genfromtxt(f_ijok,delimiter=',',skip_header=1,unpack=True)[3]
    fit_rmse = np.genfromtxt(f_ijok,delimiter=',',skip_header=1,unpack=True)[4]
    nok = len(iok)

# ----------------------------------------------------
# Parameters
# ------------------------------
    print 'Gathering the parameter samples...'
    print

    with open(os.getcwd()+'/'+MCname+'_parameters.txt','w') as f_out:
        i = 0
        while i < nok:
            if i==0 or (i>0 and iok[i]!=iok[i-1]):
                f_in = open(outdir+str(iok[i]+1)+'/parameters.txt','r')
                pnames = f_in.readline()#.rstrip().split(',')
                j = 0
                if i==0:
                    f_out.write('job,sample,'+pnames)
        
            data = f_in.readline()
            while j < jok[i]:
                data = f_in.readline()
                j+=1
        
            f_out.write(str(iok[i]+1)+','+str(jok[i]+1)+','+data)
            if i==nok-1 or iok[i]!=iok[i+1]:
                f_in.close()

            j+=1
            i+=1

# ------------------------------------------------------
# Time series of the best runs
# ----------------------------
    print 'Fetch',nbest,'best runs'
    print

    # Best runs selections
    idsort_nse = np.argsort(np.asanyarray(fit_nse))[::-1] # inverse (high = better)
    idsort_kge = np.argsort(np.asanyarray(fit_kge))[::-1] # inverse (high = better)
    idsort_rmse = np.argsort(np.asanyarray(fit_rmse)) # normal (lower = better)

    ibest_nse = [iok[i] for i in idsort_nse[range(nbest)]]
    ibest_kge = [iok[i] for i in idsort_kge[range(nbest)]]
    ibest_rmse = [iok[i] for i in idsort_rmse[range(nbest)]]

    jbest_nse = [jok[i] for i in idsort_nse[range(nbest)]]
    jbest_kge = [jok[i] for i in idsort_kge[range(nbest)]]
    jbest_rmse = [jok[i] for i in idsort_rmse[range(nbest)]]

    print 'NSE-wise...'
    with open(os.getcwd()+'/'+MCname+'_Streamflow_'+str(nbest)+'bestruns_NSE.txt','w') as f_out:
        f_out.write('NSE,'+','.join([str(it+1) for it in range(lobs)])+'\n')
        for i in range(nbest):
            it = '%03i' % int(jbest_nse[i]+1)
            f_in = outdir+str(ibest_nse[i]+1)+'/'+it+'/'+obsnames[0]+'.tab'
            sim = list(np.genfromtxt(f_in,delimiter='\t',skip_header=nts+3,unpack=True)[1])
            tmp = str(fit_nse[idsort_nse[i]])
            f_out.write(tmp+','+','.join([str(val) for val in sim])+'\n')
            
    print 'KGE-wise...'
    with open(os.getcwd()+'/'+MCname+'_Streamflow_'+str(nbest)+'bestruns_KGE.txt','w') as f_out:
        f_out.write('KGE,'+','.join([str(it+1) for it in range(lobs)])+'\n')
        for i in range(nbest):
            it = '%03i' % int(jbest_kge[i]+1)
            f_in = outdir+str(ibest_kge[i]+1)+'/'+it+'/'+obsnames[0]+'.tab'
            sim = list(np.genfromtxt(f_in,delimiter='\t',skip_header=nts+3,unpack=True)[1])
            tmp = str(fit_kge[idsort_kge[i]])
            f_out.write(tmp+','+','.join([str(val) for val in sim])+'\n')

    print 'RMSE-wise...'
    with open(os.getcwd()+'/'+MCname+'_Streamflow_'+str(nbest)+'bestruns_RMSE.txt','w') as f_out:
        f_out.write('RMSE,'+','.join([str(it+1) for it in range(lobs)])+'\n')
        for i in range(nbest):
            it = '%03i' % int(jbest_rmse[i]+1)
            f_in = outdir+str(ibest_rmse[i]+1)+'/'+it+'/'+obsnames[0]+'.tab'
            sim = list(np.genfromtxt(f_in,delimiter='\t',skip_header=nts+3,unpack=True)[1])
            tmp = str(fit_rmse[idsort_rmse[i]])
            f_out.write(tmp+','+','.join([str(val) for val in sim])+'\n')
            
        

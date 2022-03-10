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
nbest = 100

# -- Output directory
outdir = 'Outputs'

# --- Definition ---
# Path to data
obsdir = '/users/s01ks6/Data/BB'
# Simulation subdirectories
MCname = 'MC2'
# Observations
obsnames = ['Streamflow','SWC_1','SWC_2','SWC_3']#,'SWC_forest','T_forest']
nobs = len(obsnames)
# Corresponding simulations outputs
simdir = os.getcwd()+'/'+MCname+'_sampling.'
simfiles = ['Streamflow.txt','SoilMoistureAv.pt2.tab','SoilMoistureAv.pt3.tab',
          'SoilMoistureAv.pt4.tab','SoilMoistureAv.pt5.tab','Transpiration[0].pt5.tab']
# Unit conversion sim --> obs
simfct = [1,1,1,1,1,86400000] # sim T is in m/s --> mm/d
# Obs filenames
obsfiles = ['/Q/BB_discharge_daily_01062011-04072016.csv',
            '/SWC/VSM_Site1_daily_01062011-07092016.csv',
            '/SWC/VSM_Site2_daily_20042011-27092016.csv',
            '/SWC/VSM_Site3_daily_21042011-14102016.csv',
            '/SWC/VSM_Forest_daily_24022015-05102016.csv',
            '/Transpiration/forest_transpiration_08072015-28092015.csv']
obsbeg = [datetime(2011,6,1),datetime(2011,6,1),datetime(2011,4,20),
          datetime(2011,4,21),datetime(2015,2,24),datetime(2015,7,8)]
obsend = [datetime(2016,7,4),datetime(2016,9,9),datetime(2016,9,27),
          datetime(2016,10,14),datetime(2016,10,5),datetime(2015,9,28)]
obscol = [1,1,4,4,4,1]
# Number of samples
nit = 400
# There is a bug in the outputs (for now), so we skip a few runs...
rngit = range(nit)#range(27,nit)
# Number of parallel runs
njob = 96
# Number of locations for simulated time series
nts = 5

# Time frames used for the metrics
fitbeg = [datetime(2014,7,1),datetime(2014,7,1),datetime(2014,7,1),
       datetime(2014,7,1),datetime(2015,2,24),datetime(2015,7,8)]
fitend = [datetime(2016,6,30),datetime(2016,6,30),datetime(2016,6,30),
       datetime(2016,6,30),datetime(2016,6,30),datetime(2015,9,28)]

# Length of the simulated outputs
simlen = 1827
simdate = [datetime(2016,6,30)-timedelta(days=x) for x in range(simlen)[::-1]]

# Metrics
metrics = ['NSE','KGE','RMSE','Corr']#,'rSTD']
idmet = [2,3,4,6]
nmet = len(metrics)

#############################################################################################
# Simulations outputs : saving the metrics
if switch == 0:

    print 'Summarizing the metrics....'
    print

    obs={}
    iok = []
    jok = []

    ### For each obs -------------------------------------------------------------------------
    for iobs in range(nobs):
        print '           for '+obsnames[iobs]+': ',

        # -- Get the obs
        file = obsdir+obsfiles[iobs]
        tmp = np.genfromtxt(file,delimiter=',',skip_header=1,unpack=True, dtype= '|S10')[0]
        date = [datetime.strptime(a, '%d/%m/%Y') for a in tmp]
        tmp = np.genfromtxt(file,delimiter=',',skip_header=1,unpack=True)[obscol[iobs]]
        # Reference simulation date
        lobs = (obsend[iobs]-obsbeg[iobs]).days+1
        date = [obsend[iobs]-timedelta(days=x) for x in range(lobs)[::-1]]
        # Crop obs between desired time frame
        obs = [tmp[idx] for idx in range(len(tmp)) if date[idx]>=fitbeg[iobs] and date[idx]<=fitend[iobs]]

        # -- Metrics outputs
        with open(os.getcwd()+'/'+outdir+'/'+MCname+'_metrics_'+obsnames[iobs]+'.txt','w') as f_out:
            f_out.write('job,sample,'+','.join(metrics)+'\n')
            for i in range(njob):
                print 'Job', i+1,'...',
                # -- Open the simulation outfile, summarized over all the iterations
                f_in = simdir+str(i+1)+'/'+simfiles[iobs]
                sim = np.genfromtxt(f_in,delimiter=',',skip_header=1)
                for j in range(sim.shape[0]):
                    jcomp = int(sim[j][0])
                    # Get sim outpts and apply conversion factor !!
                    tmp = [simfct[iobs]*sim[j][x] for x in range(1,simlen+1)]
                    # Crop between desired time frame
                    sim2 = [tmp[idx] for idx in range(simlen) if simdate[idx]>=fitbeg[iobs] and simdate[idx]<=fitend[iobs]]
                    # Metrics
                    md_nse = metrics.nash_sutcliff(sim2,obs)
                    md_kge = metrics.kling_gupta(sim2,obs,method='2012')
                    md_rmse = metrics.rmse(sim2,obs)
                    #md_bias = metrics.bias(sim2,obs)
                    md_corr = metrics.corr(sim2,obs)
                    md_rstd = metrics.rstd(sim2,obs)
                    # Write
                    f_out.write(','.join([str(i+1),str(jcomp),str(md_nse),str(md_kge),
                                          str(md_rmse),str(md_corr),str(md_rstd)])+'\n')
                    # Save the 'coordinates' of this sample (only the first time)
                    if iobs==0:
                        iok.append(i)
                        jok.append(jcomp)

            print
            print

    nok = len(iok)
    print 'Number of complete samples:',nok
    print

#############################################################################################
# Get parameters / best runs

if switch == 1:

    # Read prior RMSEs
    f_prior=os.getcwd()+'/'+outdir+'/Prior_metrics.csv'
    RMSE_prior = np.genfromtxt(f_prior,delimiter=',',skip_header=2)[2]

    # Locate indexes
    f_ijok = os.getcwd()+'/'+outdir+'/'+MCname+'_metrics_'+obsnames[0]+'.txt'
    iok = np.genfromtxt(f_ijok,delimiter=',',skip_header=1,unpack=True,dtype=np.int64)[0]
    jok = np.asarray((np.genfromtxt(f_ijok,delimiter=',',skip_header=1,unpack=True,dtype=np.int64)[1]),dtype=int)
    nok = len(iok)

# ------------------------------
# Parameters

    print 'Gathering the parameter samples...'
    print

    with open(os.getcwd()+'/'+outdir+'/'+MCname+'_parameters.txt','w') as f_out:
        i = 0
        while i < nok:
            if i==0 or (i>0 and iok[i]!=iok[i-1]):
                f_in = open(simdir+str(iok[i])+'/parameters.txt','r')
                pnames = f_in.readline()#.rstrip().split(',')
                j = 1
                if i==0:
                    f_out.write('job,sample,'+pnames)
        
            data = f_in.readline()
            while j < jok[i]:
                data = f_in.readline()
                j+=1
        
            f_out.write(str(iok[i])+','+str(jok[i])+','+data)
            if i==nok-1 or iok[i]!=iok[i+1]:
                f_in.close()

            j+=1
            i+=1

# ------------------------------------------------------
# Time series of the best runs
# ----------------------------

    print 'Observations: '
    print

    # Metrics combining several obs -------------------------------------------------------------------------
    # Sqrt of sum of mse weighted by prior mse = cost function
    print 'Total RMSE using prior RMSE as normalization '
    print
    for iobs in range(nobs):
        f_in = os.getcwd()+'/'+outdir+'/'+MCname+'_metrics_'+obsnames[iobs]+'.txt'
        fit_rmse = np.genfromtxt(f_in,delimiter=',',skip_header=1,unpack=True)[4]
        if iobs==0:
            costJ = (fit_rmse/RMSE_prior[iobs])**2
            indi = np.genfromtxt(f_in,delimiter=',',skip_header=1,unpack=True,dtype=np.int64)[0]
            indj = np.genfromtxt(f_in,delimiter=',',skip_header=1,unpack=True,dtype=np.int64)[1]
        if iobs>0:
            costJ = costJ + (fit_rmse/RMSE_prior[iobs])**2
    # Output a file with the values of J
    with open(os.getcwd()+'/'+outdir+'/'+MCname+'_metrics_CostFunction.txt','w') as f_out:
        f_out.write('job,sample,J\n')
        for i in range(len(indi)):
            f_out.write(','.join([str(indi[i]),str(indj[i]),str(costJ[i])])+'\n')

    # Best runs selections -----------------------------------------
    print 'Fetch',nbest,' best runs'
    
    # Using the 'cost function'
    print '   -Using the cost function...'

    idsort = np.argsort(np.asanyarray(costJ)) # normal (lower = better)
    ibest = [iok[i] for i in idsort[range(nbest)]]
    jbest = [jok[i] for i in idsort[range(nbest)]]

    for iobs in range(nobs):
        with open(os.getcwd()+'/'+outdir+'/'+MCname+'_CostFunction.'+str(nbest)+'best_'+obsnames[iobs]+'_Ts.txt','w') as f_out:
            f_out.write('J,'+','.join([str(it+1) for it in range(simlen)])+'\n')
            for i in range(nbest):
                f_in = simdir+str(ibest[i])+'/'+simfiles[iobs]
                js = np.genfromtxt(f_in,delimiter=',',skip_header=1,unpack=True)[0]
                tmp = np.genfromtxt(f_in,delimiter=',',skip_header=1)
                fit2 = str(costJ[idsort[i]])
                for k in range(nit):
                    if js[k]==jbest[i]:
                        sim = [simfct[iobs]*tmp[k][x] for x in range(1,simlen+1)]
                        f_out.write(fit2+','+','.join([str(val) for val in sim])+'\n')
                        break


    # Across the metrics
    # for im in range(nmet):
    for im in [2,0,1,3]:
        print
        print '   -'+metrics[im]+'-based metrics:',
        # Using metrics on each type of obs
        for iobs in range(nobs):
            #print '   -'+obsnames[iobs]+'-based metrics:',
            print ' using '+obsnames[iobs]+'...',
            f_in = os.getcwd()+'/'+outdir+'/'+MCname+'_metrics_'+obsnames[iobs]+'.txt'
            fit  = np.genfromtxt(f_in,delimiter=',',skip_header=1,unpack=True)[idmet[im]]
            if metrics[im] in ['NSE','KGE','Corr']:
                idsort = np.argsort(np.asanyarray(fit))[::-1] # inverse (high = better)
            else:
                idsort = np.argsort(np.asanyarray(fit)) # normal (lower = better)

            ibest = [iok[i] for i in idsort[range(nbest)]]
            jbest = [jok[i] for i in idsort[range(nbest)]]

            #print metrics[im]+'...',
            for iobs2 in range(nobs):
                with open(os.getcwd()+'/'+outdir+'/'+MCname+'_'+obsnames[iobs]+'_'+metrics[im]+'.'+str(nbest)+'best_'+obsnames[iobs2]+'_Ts.txt','w') as f_out:
                    f_out.write(metrics[im]+','+','.join([str(it+1) for it in range(simlen)])+'\n')
                    for i in range(nbest):
                        f_in = simdir+str(ibest[i])+'/'+simfiles[iobs2]
                        js = np.genfromtxt(f_in,delimiter=',',skip_header=1,unpack=True)[0]
                        tmp = np.genfromtxt(f_in,delimiter=',',skip_header=1)
                        fit2 = str(fit[idsort[i]])
                        for k in range(nit):
                            if js[k]==jbest[i]:
                                sim = [simfct[iobs2]*tmp[k][x] for x in range(1,simlen+1)]
                                f_out.write(fit2+','+','.join([str(val) for val in sim])+'\n')
                                break

    

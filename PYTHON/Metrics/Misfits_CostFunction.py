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
nbest = 50

# -- Output directory
outdir = 'Outputs'

# --- Definition ---
# Path to data
obsdir = '/users/s01ks6/Data/BB'
# Simulation subdirectories
MCname = 'MC4'
# Observations
obsnames = ['Streamflow','SaturationArea','SWC_Peat','SWC_Gley','SWC_Podzol',
            'SWC_NFforest','SWC_SFforest','T_NFforest']#,'SWC_forest','T_forest']
nobs = len(obsnames)
obsnames2 = obsnames +['NetRad_WS1_bare','NetRad_WS1_heather','NetRad_WS1_molinia',
                       'NetRad_WS1_pine','NetRad_WS1_sphagnum',
                       'NetRad_WS2_bare','NetRad_WS2_heather','NetRad_WS2_molinia',
                       'NetRad_WS2_pine','NetRad_WS2_sphagnum',
                       'NetRad_WS3_bare','NetRad_WS3_heather','NetRad_WS3_molinia',
                       'NetRad_WS3_pine','NetRad_WS3_sphagnum']
nobs2 = len(obsnames2)

# Corresponding simulations outputs
simdir = os.getcwd()+'/'+MCname+'_sampling.'
simfiles = ['Streamflow_hist.tab','SaturationArea_hist.tab',
            'SWC_Peat_hist.tab','SWC_Gley_hist.tab','SWC_Podzol_hist.tab',
            'SWC_NFforest_hist.tab','SWC_SFforest_hist.tab','T_NFforest_hist.tab']
simfiles2 = simfiles + ['NetRad_WS1_bare_hist.tab','NetRad_WS1_heather_hist.tab',
                        'NetRad_WS1_molinia_hist.tab','NetRad_WS1_pine_hist.tab',
                        'NetRad_WS1_sphagnum_hist.tab',
                        'NetRad_WS2_bare_hist.tab','NetRad_WS2_heather_hist.tab',
                        'NetRad_WS2_molinia_hist.tab','NetRad_WS2_pine_hist.tab',
                        'NetRad_WS2_sphagnum_hist.tab',
                        'NetRad_WS3_bare_hist.tab','NetRad_WS3_heather_hist.tab',
                        'NetRad_WS3_molinia_hist.tab','NetRad_WS3_pine_hist.tab',
                        'NetRad_WS2_sphagnum_hist.tab']
sim_order = [6,1,8,2,5,3,7,4,9]
# Length of the simulated outputs
lsim = 1827
simt = [datetime(2011,6,1)+timedelta(days=x) for x in range(lsim)]

prdir = os.getcwd()+'/Results_prior'
prfiles = ['Streamflow.tab','BasinSummary.txt',
            'SoilMoistureAv.tab','SoilMoistureAv.tab','SoilMoistureAv.tab',
            'SoilMoistureAv.tab','SoilMoistureAv.tab','Transpiration[0].tab']
prcol = [1,9,2,3,4,5,6,5]
# Number of locations for simulated time series
nts = 9#len(prcol)
# Unit conversion sim --> obs
simfct = [1,1,1,1,1,1,1,8.64e7,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1] # sim T is in m/s --> mm/d
# Obs filenames
obsfiles = ['/Q/BB_discharge_daily_01062011-20092016.csv',
            '/SatExt/SatExt_daily_01062011-30092014.csv',
            '/SWC/VSM_Site1_daily_01062011-07092016.csv',
            '/SWC/VSM_Site2_daily_20042011-27092016.csv',
            '/SWC/VSM_Site3_daily_21042011-14102016.csv',
            '/SWC/VSM_Forest_daily_24022015-05102016.csv',
            '/SWC/VSM_Plantation_daily_10102012-14102016.csv',
            '/Transpiration/forest_transpiration_08072015-28092015.csv']
#obsbeg = [datetime(2011,6,1),datetime(2011,6,1),datetime(2011,4,20),
#          datetime(2011,4,21),datetime(2015,2,24),datetime(2015,7,8)]
#obsend = [datetime(2016,7,4),datetime(2016,9,9),datetime(2016,9,27),
#          datetime(2016,10,14),datetime(2016,10,5),datetime(2015,9,28)]
obscol = [1,1,1,4,4,4,4,1]
# Number of samples
nit = 200
# There is a bug in the outputs (for now), so we skip a few runs...
rngit = range(nit)#range(27,nit)
# Number of parallel runs
njob = 50

# Time frames used for the metrics
fitbeg = [datetime(2012,11,1),datetime(2012,10,1),datetime(2012,11,1),
          datetime(2012,11,1),datetime(2012,11,1),datetime(2015,2,24),
          datetime(2012,11,1),datetime(2015,7,8)]
fitend = [datetime(2014,10,31),datetime(2014,9,30),datetime(2014,10,31),
          datetime(2014,10,31),datetime(2014,10,31),datetime(2016,2,23),
          datetime(2014,10,31),datetime(2015,9,28)]

# Metrics
#metrics = ['CF_MSE','CF_NSE']
#nmet = len(metrics)

#############################################################################################
# Simulations outputs : saving the metrics
if switch == 0:

    print 'Summarizing the metrics....'
    print

    obs={}
    sigma = {}
    nok = 0

    print 'Get obs and prior misfit...'
    print
    ### For each obs -------------------------------------------------------------------------
    for iobs in range(nobs):

        # -- Get the obs
        f_obs = obsdir+obsfiles[iobs]
        tmp = np.genfromtxt(f_obs,delimiter=',',skip_header=1,unpack=True, dtype= '|S10')[0]
        obst = [datetime.strptime(a, '%d/%m/%Y') for a in tmp]
        lobs = len(obst)
        tmp = np.genfromtxt(f_obs,delimiter=',',skip_header=1,unpack=True)[obscol[iobs]]
        # Crop obs between desired time frame
        obs[obsnames[iobs]] = [tmp[idx] for idx in range(lobs) if obst[idx]>=fitbeg[iobs] and obst[idx]<=fitend[iobs]]
        # print obs[obsnames[iobs]]

        # -- Prior (reference)
        # Some reformatting of the original outputs...
        if obsnames[iobs] == 'SaturationArea':
            hskip=1
            idx = prcol[iobs]-1
        else:
            hskip= nts+3
            idx = np.argsort(np.array(sim_order))[prcol[iobs]-1]+1

        tmp = (np.genfromtxt(prdir+'/'+prfiles[iobs],delimiter='\t',
                            skip_header=hskip,unpack=True)[idx])*simfct[iobs]

        # f_sim = prdir+'/'+simfiles[iobs]
        # tmp = np.genfromtxt(f_sim,delimiter=',',skip_header=1)*simfct[iobs]
        prior = [tmp[idx] for idx in range(1729) if simt[idx]>=fitbeg[iobs] and simt[idx]<=fitend[iobs]]
        # print prior
        sigma[obsnames[iobs]] = metrics.mse(prior,obs[obsnames[iobs]])
    # print sigma
    print

    print 'Get posterior misfits...'
    # Loop over jobs
    with open(os.getcwd()+'/'+outdir+'/'+MCname+'_costfunctions.txt','w') as f_out:        
        f_out.write('job,sample,CF_MSE,CF_NSE\n')

    for i in range(njob):
        
        print 'Job', i+1,'...',

        # Now misfit across obs
        for iobs in range(nobs):
            f_sim = simdir+str(i+1)+'/'+simfiles[iobs]
            # Apply conversion factor
            tmp = np.genfromtxt(f_sim,delimiter=',',skip_header=1)

            if iobs==0:
                # Get the number of runs that worked
                if len(tmp.shape)==1:
                    nj = 1
                else:
                    nj = tmp.shape[0]
                print 'nj =',nj,':'
                # initialize costfunction
                Jmse = [0]*nj
                Jnse = [0]*nj
                iok = []
                jok = []

            for j in range(nj):

                if nj==1:
                    tmp2 = copy.copy(tmp)
                else:
                    tmp2 = tmp[j]

                # Get the run number and save the 'coordinates' of this sample
                if iobs==0:
                    jcomp = int(tmp2[0])
                    iok += [i+1]
                    jok += [jcomp]

                # Get the right sim outpts 
                sim = [tmp2[x]*simfct[iobs] for x in range(1,lsim+1)]
                # Crop between desired time frame
                sim = [sim[idx] for idx in range(lsim) if simt[idx]>=fitbeg[iobs] and simt[idx]<=fitend[iobs]]
                lfit = len(sim)
                # Increment cost function
                Jmse[j] += lfit*metrics.mse(sim,obs[obsnames[iobs]])/sigma[obsnames[iobs]]
                Jnse[j] += (1-metrics.nash_sutcliff(sim,obs[obsnames[iobs]]))**2


        # Write
        with open(os.getcwd()+'/'+outdir+'/'+MCname+'_costfunctions.txt','a') as f_out:
            for j in range(nj):
                f_out.write(','.join([str(iok[j]),str(jok[j]),str(Jmse[j]),str(Jnse[j])])+'\n')

        print
        nok += len(iok)

    print 'Number of complete samples:',nok
    print

#############################################################################################
# Get parameters / best runs

if switch == 1:

    # Locate indexes
    f_ijok = os.getcwd()+'/'+outdir+'/'+MCname+'_costfunctions.txt'
    iok = np.genfromtxt(f_ijok,delimiter=',',skip_header=1,unpack=True,dtype=int)[0]
    jok = np.genfromtxt(f_ijok,delimiter=',',skip_header=1,unpack=True,dtype=int)[1]
    nok = len(iok)
    #print iok
    #print jok
    print nok

# ------------------------------
# Parameters

    print 'Gathering the parameter samples...'
    print

    nok2 = 0
    with open(os.getcwd()+'/'+outdir+'/'+MCname+'_parameters.txt','w') as f_out:
        k = 0
        for i in range(1,njob+1):
            print i,
            f_in = simdir+str(i)+'/Parameters.txt'
            if i==1:
                tmp = list(np.genfromtxt(f_in,delimiter=',',dtype= '|S20')[0])
                pnames = ','.join([str(x) for x in tmp])
                f_out.write('job,sample,'+pnames+'\n')
        
            tmp = np.genfromtxt(f_in,delimiter=',',skip_header=1)
            #print tmp.shape
            if len(tmp.shape)==1:
                nok2+= 1
                print 'Single sample for this job!!'
                tmp2 = list(tmp)
                f_out.write(str(i)+','+','.join([str(x) for x in tmp2])+'\n')
                k+=1

            else:
                nj = tmp.shape[0]
                nok2+=nj
                print nj, nok2
                for j in range(nj):
                    tmp2 = list(tmp[j])
                    f_out.write(str(i)+','+','.join([str(x) for x in tmp2])+'\n')
                    #print i, jok[k], k, j, nj
                    #print k,
                    k+=1

    print nok2
# ------------------------------------------------------
# Time series of the best runs
# ----------------------------

    print 'Observations: '
    print

    # Metrics combining several obs -------------------------------------------------------------------------
    # Sqrt of sum of mse weighted by prior mse = cost function
    f_in = os.getcwd()+'/'+outdir+'/'+MCname+'_costfunctions.txt'
    Jmse = np.genfromtxt(f_in,delimiter=',',skip_header=1,unpack=True)[2]
    Jnse = np.genfromtxt(f_in,delimiter=',',skip_header=1,unpack=True)[3]

    # Best runs selections -----------------------------------------
    print 'Fetch',nbest,' best runs'
    
    # Using the MSE J
    print '   -Using MSE cost function...'

    idsort = np.argsort(np.asanyarray(Jmse)) # normal (lower = better)
    ibest = [iok[i] for i in idsort[range(nbest)]]
    jbest = [jok[i] for i in idsort[range(nbest)]]

    for i in range(nbest):
        print i,
        for iobs in range(nobs2):
            f_Jout = os.getcwd()+'/'+outdir+'/'+MCname+'_CFmse.'+str(nbest)+'best_'+obsnames2[iobs]+'_Ts.txt'
            if i==0:
                with open(f_Jout,'w') as f_out:
                    f_out.write('Jmse,'+','.join([str(it+1) for it in range(lsim)])+'\n')

            # Get simulation
            f_in = simdir+str(ibest[i])+'/'+simfiles2[iobs]
            js = np.genfromtxt(f_in,delimiter=',',skip_header=1,unpack=True)[0]
            if len(js)==1:
                tmp = np.genfromtxt(f_in,delimiter=',',skip_header=1)
                sim = [simfct[iobs]*tmp[x] for x in range(1,lsim+1)]
            else:
                k=0
                while(js[k]!=jbest[i]): k+=1
                tmp = np.genfromtxt(f_in,delimiter=',',skip_header=1)
                sim = [simfct[iobs]*tmp[k][x] for x in range(1,lsim+1)]
            # Add to list
            with open(f_Jout,'a') as f_out:
                f_out.write(str(Jmse[idsort[i]])+','+','.join([str(val) for val in sim])+'\n')

    print

    # Using the NSE J
    print '   -Using NSE cost function...'

    idsort = np.argsort(np.asanyarray(Jnse)) # normal (lower = better)
    ibest = [iok[i] for i in idsort[range(nbest)]]
    jbest = [jok[i] for i in idsort[range(nbest)]]

#    for iobs in range(nobs2):
#        print obsnames2[iobs]
#        with open(os.getcwd()+'/'+outdir+'/'+MCname+'_CFnse.'+str(nbest)+'best_'+obsnames2[iobs]+'_Ts.txt','w') as f_out:
#            f_out.write(','.join([str(it+1) for it in range(lsim)])+'\n')
#            for i in range(nbest):
#                f_in = simdir+str(ibest[i])+'/'+simfiles2[iobs]
#                js = np.genfromtxt(f_in,delimiter=',',skip_header=1,unpack=True)[0]
#               tmp = np.genfromtxt(f_in,delimiter=',',skip_header=1)
#               if len(js)==1:
#                   sim = [simfct[iobs]*tmp[x] for x in range(1,lsim+1)]
#               else:
#                   k=0
#                   while(js[k]!=jbest[i]): k+=1
#                  sim = [simfct[iobs]*tmp[k][x] for x in range(1,lsim+1)]
#               f_out.write(','.join([str(val) for val in sim])+'\n')

    
    for i in range(nbest):
        print i,
        for iobs in range(nobs2):
            f_Jout = os.getcwd()+'/'+outdir+'/'+MCname+'_CFnse.'+str(nbest)+'best_'+obsnames2[iobs]+'_Ts.txt'
            if i==0:
                with open(f_Jout,'w') as f_out:
                    f_out.write('Jnse,'+','.join([str(it+1) for it in range(lsim)])+'\n')

            # Get simulation
            f_in = simdir+str(ibest[i])+'/'+simfiles2[iobs]
            js = np.genfromtxt(f_in,delimiter=',',skip_header=1,unpack=True)[0]
            if len(js)==1:
                tmp = np.genfromtxt(f_in,delimiter=',',skip_header=1)
                sim = [simfct[iobs]*tmp[x] for x in range(1,lsim+1)]
            else:
                k=0
                while(js[k]!=jbest[i]): k+=1
                tmp = np.genfromtxt(f_in,delimiter=',',skip_header=1)
                sim = [simfct[iobs]*tmp[k][x] for x in range(1,lsim+1)]
            # Add to list
            with open(f_Jout,'a') as f_out:
                f_out.write(str(Jmse[idsort[i]])+','+','.join([str(val) for val in sim])+'\n')

    print

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
# 0 = metrics + par, 1 = best runs
switch = 1
nbest = 30

# -- Output directory
outdir = 'Outputs'

# --- Definition ---
# Path to data
obsdir = '/users/s01ks6/Data/BB'
# Simulation subdirectories
MCname = 'MC6'
# Observations
obsnames = ['Streamflow','SaturationArea','SWC_Peat','SWC_Gley','SWC_Podzol',
            'SWC_NFforest','SWC_SFforest','T_NFforest','T_SFforest']#,'SWC_forest','T_forest']
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
            'SWC_NFforest_hist.tab','SWC_SFforest_hist.tab','T_NFforest_hist.tab','T_SFforest_hist.tab']
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
lsim = 1888
simt = [datetime(2011,6,1)+timedelta(days=x) for x in range(lsim)]

prdir = os.getcwd()+'/Results_prior'
prfiles = ['Streamflow.tab','BasinSummary.txt',
            'SoilMoistureAv.tab','SoilMoistureAv.tab','SoilMoistureAv.tab',
            'SoilMoistureAv.tab','SoilMoistureAv.tab','Transpiration[0].tab','Transpiration[0].tab']
prcol = [1,9,2,3,4,5,6,5,6]
# Number of locations for simulated time series
nts = 9#len(prcol)
# Unit conversion sim --> obs
simfct = [1,1,1,1,1,1,1,8.64e7,8.64e7,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1] # sim T is in m/s --> mm/d
# Obs filenames
obsfiles = ['/Q/BB_discharge_daily_01062011-20092016.csv',
            '/SatExt/SatExt_daily_01062011-30092014.csv',
            '/SWC/VSM_Site1_daily_01062011-07092016.csv',
            '/SWC/VSM_Site2_daily_20042011-27092016.csv',
            '/SWC/VSM_Site3_daily_21042011-14102016.csv',
            '/SWC/VSM_Forest_daily_24022015-05102016.csv',
            '/SWC/VSM_Plantation_daily_10102012-14102016.csv',
            '/Transpiration/T_NFforest_08072015-28092015.csv',
            '/Transpiration/T_SFforest_01042016-22092016.csv']
#obsbeg = [datetime(2011,6,1),datetime(2011,6,1),datetime(2011,4,20),
#          datetime(2011,4,21),datetime(2015,2,24),datetime(2015,7,8)]
#obsend = [datetime(2016,7,4),datetime(2016,9,9),datetime(2016,9,27),
#          datetime(2016,10,14),datetime(2016,10,5),datetime(2015,9,28)]
obscol = [1,1,1,4,4,4,4,1,1]
# Number of samples
nit = 400
# There is a bug in the outputs (for now), so we skip a few runs...
rngit = range(nit)#range(27,nit)
# Number of parallel runs
njob = 100

# Time frames used for the metrics
fitbeg = [datetime(2012,11,1),datetime(2012,10,1),datetime(2012,11,1),
          datetime(2012,11,1),datetime(2012,11,1),datetime(2015,2,24),
          datetime(2012,11,1),datetime(2015,7,8),datetime(2016,4,1)]
fitend = [datetime(2014,10,31),datetime(2014,9,30),datetime(2014,10,31),
          datetime(2014,10,31),datetime(2014,10,31),datetime(2016,2,23),
          datetime(2014,10,31),datetime(2015,9,28),datetime(2016,7,31)]
# Metrics
#metrics = ['CF_MSE','CF_NSE']
#nmet = len(metrics)

#############################################################################################
# Simulations outputs : saving the metrics
# if switch == 0:

print 'Summarizing the metrics....'
print

obs={}
sigma = {}
pr_nse = {}
nok = 0
CFMSEprior = 0
CFNSEprior = 0

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
    prior = [tmp[idx] for idx in range(lsim) if simt[idx]>=fitbeg[iobs] and simt[idx]<=fitend[iobs]]
    # print prior
    sigma[obsnames[iobs]] = metrics.mse(prior,obs[obsnames[iobs]])
    pr_nse[obsnames[iobs]] = metrics.nash_sutcliff(prior,obs[obsnames[iobs]])
    CFMSEprior += sigma[obsnames[iobs]]*len(prior)
    CFNSEprior += (1-pr_nse[obsnames[iobs]])**2

    print obsnames[iobs], len(obs[obsnames[iobs]]), len(prior)
print sigma
print pr_nse

tmp = ','.join(['CF_'+x for x in obsnames])
with open(os.getcwd()+'/'+outdir+'/'+MCname+'_CFmse_prior.txt','w') as f_out:
    f_out.write(tmp+',CF_tot\n')
    f_out.write(','.join([','.join([str(sigma[x]) for x in obsnames]),str(CFMSEprior)])+'\n')
with open(os.getcwd()+'/'+outdir+'/'+MCname+'_CFnse_prior.txt','w') as f_out:        
    f_out.write('Sample,'+tmp+',CF_tot\n')
    f_out.write(','.join([','.join([str(pr_nse[x]) for x in obsnames]),str(np.sqrt(CFNSEprior))])+'\n')


if switch == 0:

    print 'Get posterior misfits...'
    # Loop over jobs
    tmp = ','.join(['CF_'+x for x in obsnames])
    with open(os.getcwd()+'/'+outdir+'/'+MCname+'_CFmse.txt','w') as f_out:
        f_out.write('Sample,'+tmp+',CF_tot\n')
    with open(os.getcwd()+'/'+outdir+'/'+MCname+'_CFnse.txt','w') as f_out:        
        f_out.write('Sample,'+tmp+',CF_tot\n')
    with open(os.getcwd()+'/'+outdir+'/'+MCname+'_parameters.txt','w') as f_out:
        f_in = simdir+'1/Parameters.txt'
        tmp = list(np.genfromtxt(f_in,delimiter=',',dtype= '|S20')[0])
        npar = len(tmp)-1
        print npar
        pnames = ','.join([str(tmp[idx]) for idx in range(1,npar+1)])
        f_out.write('Job,Iteration,Sample,'+pnames+'\n')

    # initialize costfunction
    Jmse = {}
    Jnse = {}

    itot = 1
    itot2 = 1
    itot3 = 1
    for i in range(1,njob+1):
        
        print 'Job', i,'...',

        for iobs in range(nobs):
            f_sim = simdir+str(i)+'/'+simfiles[iobs]
            tmp = np.genfromtxt(f_sim,delimiter=',',skip_header=1)

            if iobs==0:
                # Get the number of runs that worked
                js = np.genfromtxt(f_sim,delimiter=',',skip_header=1,unpack=True,dtype=np.int)[0]
                js = [js[idx] for idx in range(len(js)) if js[idx]<=nit]
                #if len(tmp.shape)==1:
                #    nj = 1
                #else:
                nj = len(js)
                print 'nj =',nj,':',
                # initialize costfunction
                Jmsetot = [0]*nj
                Jnsetot = [0]*nj

                f_par = simdir+str(i)+'/Parameters.txt'
                tmp_par = np.genfromtxt(f_par,delimiter=',',skip_header=1)

            # initialize costfunction
            Jmse[obsnames[iobs]] = [0]*nj
            Jnse[obsnames[iobs]] = [0]*nj

            for j in range(1,nj+1):
                if iobs==0:
                    print js[j-1],
                # -- Cost functions
                #if nj==1:
                #    tmp2 = copy.copy(tmp)
                #else:
                tmp2 = tmp[j-1]

                # Crop between desired time frame (conversion factor as well)
                sim = [tmp2[idx+1]*simfct[iobs] for idx in range(lsim) if simt[idx]>=fitbeg[iobs] and simt[idx]<=fitend[iobs]]
                #lfit = len(sim)
                # Increment cost function              
                Jnse[obsnames[iobs]][j-1] = metrics.nash_sutcliff(sim,obs[obsnames[iobs]])
                Jmse[obsnames[iobs]][j-1] = metrics.mse(sim,obs[obsnames[iobs]])/sigma[obsnames[iobs]]

                Jmsetot[j-1] += Jmse[obsnames[iobs]][j-1]
                Jnsetot[j-1] += (1-Jnse[obsnames[iobs]][j-1])**2#/(1-pr_nse[obsnames[iobs]]
                
                if iobs ==0:
                    # -- Parameters 
                    tmp2 = [tmp_par[j-1][idx] for idx in range(1,npar+1)]
                    with open(os.getcwd()+'/'+outdir+'/'+MCname+'_parameters.txt','a') as f_out:
                        f_out.write(str(i)+','+str(js[j-1])+','+str(itot)+','+','.join([str(x) for x in tmp2])+'\n')
                    itot+=1

            print obsnames[iobs],
        print       

        # Write
        with open(os.getcwd()+'/'+outdir+'/'+MCname+'_CFmse.txt','a') as f_out:
            for j in range(nj):
                f_out.write(','.join([str(itot2),','.join([str(Jmse[x][j-1]) for x in obsnames]),str(Jmsetot[j-1])])+'\n')
                itot2+=1
        with open(os.getcwd()+'/'+outdir+'/'+MCname+'_CFnse.txt','a') as f_out:
            for j in range(nj):
                f_out.write(','.join([str(itot3),','.join([str(Jnse[x][j-1]) for x in obsnames]),str(np.sqrt(Jnsetot[j-1]))])+'\n')
                itot3+=1

    print 'Number of complete samples:',itot, itot2
    print

#############################################################################################
# Time series of the best runs
# ----------------------------

if switch == 1:

    # Metrics used
    metric = ['Streamflow','SatExt','SWC','T', # 1 constraints (type)
              'Stream&Sat','Stream&SWC','Stream&T','Sat&SWC','Sat&T','SWC&T', # 2 constraints
              'Stream&Sat&SWC','Stream&Sat&T','Stream&SWC&T','Sat&SWC&T', # 3 constraints
              'Total']#,'T']
    nmet = len(metric)

    # Locate indexes
    f_in = os.getcwd()+'/'+outdir+'/'+MCname+'_parameters.txt'
    iok = np.genfromtxt(f_in,delimiter=',',skip_header=1,unpack=True,dtype=int)[0]
    jok = np.genfromtxt(f_in,delimiter=',',skip_header=1,unpack=True,dtype=int)[1]
    #print iok
    #print jok
    f_in = os.getcwd()+'/'+outdir+'/'+MCname+'_CFmse.txt'
    Jtot = np.genfromtxt(f_in,delimiter=',',skip_header=1,unpack=True)[nobs+1]
    f_prior = os.getcwd()+'/'+outdir+'/'+MCname+'_CFmse_prior.txt'
    Jprior = np.genfromtxt(f_prior,delimiter=',',skip_header=1,unpack=True)[:nobs]

    # Using streamflow only 
    Jstream = np.genfromtxt(f_in,delimiter=',',skip_header=1,unpack=True)[1]
    # using saturation extent
    Jsatext = np.genfromtxt(f_in,delimiter=',',skip_header=1,unpack=True)[2]
    # using swc
    JSWCs = np.genfromtxt(f_in,delimiter=',',skip_header=1,unpack=True)[3:3+5]
    JSWC = JSWCs[0]+JSWCs[1]+JSWCs[2]+JSWCs[3]+JSWCs[4]
    # using transpiration
    JTNF = np.genfromtxt(f_in,delimiter=',',skip_header=1,unpack=True)[8]
    JTSF = np.genfromtxt(f_in,delimiter=',',skip_header=1,unpack=True)[9]
    JT = JTNF + JTSF
    #print Jmse
    #print Jns

    # Best runs selections -----------------------------------------
    print 'Fetch',nbest,' best runs'
    
    for im in range(nmet):
        # Using the MSE J
        if im==0:
            print '   -Using streamflow...'
            Jsort = Jstream
        if im==1:
            print '   -Using saturation extent...'
            Jsort = Jsatext
        if im==2:
            print '   -Using soil moisture...'
            Jsort = JSWC
        if im==3:
            print '   -Using transpiration...'
            Jsort = JT
        if im==4:
            print '   -Using streamflow and saturation...'
            Jsort = Jstream + Jsatext
        if im==5:
            print '   -Using streamflow and SWC...'
            Jsort = Jstream + JSWC/5
        if im==6:
            print '   -Using streamflow and transpiration...'
            Jsort = Jstream + JT/2
        if im==7:
            print '   -Using saturation and SWC...'
            Jsort = Jsatext + JSWC/5
        if im==8:
            print '   -Using saturation and transpiration...'
            Jsort = Jsatext + JT/2
        if im==9:
            print '   -Using SWC and transpiration...'
            Jsort = JSWC/5 + JT/2
        if im==10:
            print '   -Using streamflow, saturation and SWC...'
            Jsort = Jstream + Jsatext + JSWC/5
        if im==11:
            print '   -Using streamflow, saturation and transpiration...'
            Jsort = Jstream + Jsatext + JT/2
        if im==12:
            print '   -Using streamflow, SWC and transpiration...'
            Jsort = Jstream + JSWC/5 + JT/2
        if im==13:
            print '   -Using saturation, SWC and transpiration...'
            Jsort = Jsatext + JSWC/5 + JT/2
        if im==14:
            print '   -Using all...'
            Jsort = Jstream + Jsatext + JSWC/5 + JT/2


        idsort = np.argsort(np.asanyarray(Jsort)) # normal (lower = better)
        ibest = [iok[i] for i in idsort[range(nbest)]]
        jbest = [jok[i] for i in idsort[range(nbest)]]

        for i in range(nbest):
            print i,#ibest[i],jbest[i],
            for iobs in range(nobs2):
                #print obsnames2[iobs],
                f_Jout = os.getcwd()+'/'+outdir+'/'+MCname+'_CF'+metric[im]+'.'+str(nbest)+'best_'+obsnames2[iobs]+'_Ts.txt'
                if i==0:
                    with open(f_Jout,'w') as f_out:
                        f_out.write('J,'+','.join([str(it+1) for it in range(lsim)])+'\n')

                # Get simulation
                f_in = simdir+str(ibest[i])+'/'+simfiles2[iobs]
                if iobs==0:
                    js = np.genfromtxt(f_in,delimiter=',',skip_header=1,unpack=True,dtype=np.int)[0]
                    #print len(js),
                #k = [idx for idx in range(len(js)) if idx==jbest[i]-1][0] # quick fix, to be changed  with line below!!!
                    k = [idx for idx in range(len(js)) if js[idx]==jbest[i]][0]
                    #print k
                tmp = np.genfromtxt(f_in,delimiter=',',skip_header=1)[k]*simfct[iobs]
                sim = [tmp[x] for x in range(1,lsim+1)]
                # Add to list
                with open(f_Jout,'a') as f_out:
                    f_out.write(str(Jsort[idsort[i]])+','+','.join([str(val) for val in sim])+'\n')
        print

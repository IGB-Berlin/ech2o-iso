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
#nmet = len(metrics

#############################################################################################
# Simulations outputs : saving the metrics
# if switch == 0:

print 'Summarizing the metrics....'
print

obs={}
sigma = {}
#pr_nse = {}
nok = 0
#CFMSEprior = 0
#CFNSEprior = 0

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
    #pr_nse[obsnames[iobs]] = metrics.nash_sutcliff(prior,obs[obsnames[iobs]])
    #CFMSEprior += sigma[obsnames[iobs]]*len(prior)
    #CFNSEprior += (1-pr_nse[obsnames[iobs]])**2

    print obsnames[iobs], len(obs[obsnames[iobs]]), len(prior)
print sigma#[str(sigma[obs]) for obs in obsnames]
#print pr_nse

tmp = ','.join(['CF_'+x for x in obsnames])
with open(os.getcwd()+'/'+outdir+'/'+MCname+'_CFmse_prior.txt','w') as f_out:
    f_out.write(tmp+'\n')
    f_out.write(','.join([str(sigma[x]) for x in obsnames])+'\n')
#with open(os.getcwd()+'/'+outdir+'/'+MCname+'_CFnse_prior.txt','w') as f_out:        
#    f_out.write('Sample,'+tmp+'\n')
#    f_out.write(','.join([str(pr_nse[x]) for x in obsnames])+'\n')

# ------------------------------------------------------------------------------------

if switch == 0:

    print 'Get posterior misfits...'
    # Loop over jobs
    tmp = ','.join(['CF_'+x for x in obsnames])
    with open(os.getcwd()+'/'+outdir+'/'+MCname+'_CFmse.txt','w') as f_out:
        f_out.write('Sample,'+tmp+'\n')
    #with open(os.getcwd()+'/'+outdir+'/'+MCname+'_CFnse.txt','w') as f_out:        
    #    f_out.write('Sample,'+tmp+',CF_tot\n')
    with open(os.getcwd()+'/'+outdir+'/'+MCname+'_parameters.txt','w') as f_out:
        f_in = simdir+'1/Parameters.txt'
        tmp = list(np.genfromtxt(f_in,delimiter=',',dtype= '|S20')[0])
        npar = len(tmp)-1
        print npar
        pnames = ','.join([str(tmp[idx]) for idx in range(1,npar+1)])
        f_out.write('Job,Iteration,Sample,'+pnames+'\n')

    # initialize costfunction
    Jmse = {}
    #Jnse = {}

    itot = 1
    itot2 = 1
    #itot3 = 1
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
                print 'nj =',nj,':'
                # initialize costfunction
                #Jmsetot = [0]*nj
                #Jnsetot = [0]*nj

                f_par = simdir+str(i)+'/Parameters.txt'
                tmp_par = np.genfromtxt(f_par,delimiter=',',skip_header=1)

            # initialize costfunction
            Jmse[obsnames[iobs]] = [0]*nj
            #Jnse[obsnames[iobs]] = [0]*nj

            for j in range(1,nj+1):
                #if iobs==0:
                #    print js[j-1],
                # -- Cost functions
                #if nj==1:
                #    tmp2 = copy.copy(tmp)
                #else:(conversion factor as well)
                tmp2 = tmp[j-1]*simfct[iobs]

                # Crop between desired time frame 
                sim = [tmp2[idx+1] for idx in range(lsim) if simt[idx]>=fitbeg[iobs] and simt[idx]<=fitend[iobs]]
                #lfit = len(sim)
                # Increment cost function              
                #Jnse[obsnames[iobs]][j-1] = metrics.nash_sutcliff(sim,obs[obsnames[iobs]])
                Jmse[obsnames[iobs]][j-1] = metrics.mse(sim,obs[obsnames[iobs]])#/sigma[obsnames[iobs]]
                if j==1:
                    print np.mean(sim), np.mean(obs[obsnames[iobs]]), len(sim), len(obs[obsnames[iobs]]),
                #Jmsetot[j-1] += Jmse[obsnames[iobs]][j-1]
                #Jnsetot[j-1] += (1-Jnse[obsnames[iobs]][j-1])**2#/(1-pr_nse[obsnames[iobs]]
                
                if iobs ==0:
                    # -- Parameters 
                    tmp3 = [tmp_par[j-1][idx] for idx in range(1,npar+1)]
                    with open(os.getcwd()+'/'+outdir+'/'+MCname+'_parameters.txt','a') as f_out:
                        f_out.write(str(i)+','+str(j)+','+str(itot)+','+','.join([str(x) for x in tmp3])+'\n')
                    itot+=1

            print obsnames[iobs]
        print       

        # Write
        with open(os.getcwd()+'/'+outdir+'/'+MCname+'_CFmse.txt','a') as f_out:
            for j in range(1,nj+1):
                f_out.write(','.join([str(i),str(j),str(itot2),','.join([str(Jmse[x][j-1]) for x in obsnames])])+'\n')
                itot2+=1
        #with open(os.getcwd()+'/'+outdir+'/'+MCname+'_CFnse.txt','a') as f_out:
        #    for j in range(nj):
        #        f_out.write(','.join([str(itot3),','.join([str(Jnse[x][j-1]) for x in obsnames]),str(np.sqrt(Jnsetot[j-1]))])+'\n')
        #        itot3+=1

    print 'Number of complete samples:',itot, itot2
    print

#############################################################################################
# Time series of the best runs
# ----------------------------

if switch == 1:

    # Metrics used
    metric = ['Streamflow','SatExt',
              'SWC_Peat','SWC_Gley','SWC_Podzol','SWC_NF','SWC_SF','SWCall',
              'T_NF','T_SF','Tall',
              'TotalData','TotalDtype']
    #'Stream&Sat','Stream&SWC','Stream&T','Sat&SWC','Sat&T','SWC&T', # 2 constraints
    #'Stream&Sat&SWC','Stream&Sat&T','Stream&SWC&T','Sat&SWC&T', # 3 constraints
    #          'Total']#,'T']
    nmet = len(metric)

    # Locate indexes
    f_in = os.getcwd()+'/'+outdir+'/'+MCname+'_parameters.txt'
    iok = np.genfromtxt(f_in,delimiter=',',skip_header=1,unpack=True,dtype=int)[0]
    jok = np.genfromtxt(f_in,delimiter=',',skip_header=1,unpack=True,dtype=int)[1]
    sok = np.genfromtxt(f_in,delimiter=',',skip_header=1,unpack=True,dtype=int)[2]
    
    lok = len(iok)
    #print iok
    #print jok
    f_in = os.getcwd()+'/'+outdir+'/'+MCname+'_CFmse.txt'
    iok2 = np.genfromtxt(f_in,delimiter=',',skip_header=1,unpack=True,dtype=int)[0]
    jok2 = np.genfromtxt(f_in,delimiter=',',skip_header=1,unpack=True,dtype=int)[1]
    sok2 = np.genfromtxt(f_in,delimiter=',',skip_header=1,unpack=True,dtype=int)[2]
    #Jtot = np.genfromtxt(f_in,delimiter=',',skip_header=1,unpack=True)[nobs+1]
    f_prior = os.getcwd()+'/'+outdir+'/'+MCname+'_CFmse_prior.txt'
    Jprior = np.genfromtxt(f_prior,delimiter=',',skip_header=1,unpack=True)[:nobs]

    # Using streamflow only 
    Jstream = np.genfromtxt(f_in,delimiter=',',skip_header=1,unpack=True)[3]
    Jstream = Jstream/Jprior[0]
    # using saturation extent
    Jsatext = np.genfromtxt(f_in,delimiter=',',skip_header=1,unpack=True)[4]
    Jsatext = Jsatext/Jprior[1]
    # using swc
    JSWCs = np.genfromtxt(f_in,delimiter=',',skip_header=1,unpack=True)[5:5+5]
    JSWC = JSWCs[0]/Jprior[2]+JSWCs[1]/Jprior[3]+JSWCs[2]/Jprior[4]+JSWCs[3]/Jprior[5]+JSWCs[4]/Jprior[6]
    # using transpiration
    JTNF = np.genfromtxt(f_in,delimiter=',',skip_header=1,unpack=True)[10]
    JTSF = np.genfromtxt(f_in,delimiter=',',skip_header=1,unpack=True)[11]
    JT = JTNF/Jprior[7] + JTSF/Jprior[8]
    #print Jmse
    #print Jns

    # Best runs selections -----------------------------------------
    print 'Fetch',nbest,' best runs'
    
    for im in range(nmet):
        # Using the MSE J
        if im==0:
            print '   -Using streamflow...'
            Jsort = Jstream
            Jprior2 = Jprior[0]
            iobs2 = [0]
        if im==1:
            print '   -Using saturation extent...'
            Jsort = Jsatext
            Jprior2 = Jprior[1]
            iobs2 = [1]
        if im==2:
            print '   -Using soil moisture (peat)...'
            Jsort = JSWCs[0]
            Jprior2 = Jprior[2]
            iobs2 = [2,3,4,5,6]
        if im==3:
            print '   -Using soil moisture (gley)...'
            Jsort = JSWCs[1]
            Jprior2 = Jprior[3]
            iobs2 = [2,3,4,5,6]
        if im==4:
            print '   -Using soil moisture (podzol)...'
            Jsort = JSWCs[2]
            Jprior2 = Jprior[4]
            iobs2 = [2,3,4,5,6]
        if im==5:
            print '   -Using soil moisture (NF forest)...'
            Jsort = JSWCs[3]
            Jprior2 = Jprior[5]
            iobs2 = [2,3,4,5,6]
        if im==6:
            print '   -Using soil moisture (SF forest)...'
            Jsort = JSWCs[4]
            Jprior2 = Jprior[6]
            iobs2 = [2,3,4,5,6]
        if im==7:
            print '   -Using soil moisture (all)...'
            Jsort = JSWC
            Jprior2 = 5
            iobs2 = [2,3,4,5,6]
        if im==8:
            print '   -Using transpiration (NF forest)...'
            Jsort = JTNF
            Jprior2 = Jprior[7]
            iobs2 = [7,8]
        if im==9:
            print '   -Using transpiration (SF forest)...'
            Jsort = JTSF
            Jprior2 = Jprior[8]
            iobs2 = [7,8]
        if im==10:
            print '   -Using transpiration (all)...'
            Jsort = JT
            Jprior2 = 2
            iobs2 = [7,8]
        if im==11:
            print '   -Using all (equal weight of data streams)...'
            Jsort = Jstream + Jsatext + JSWC + JT
            Jprior2 = nobs
            iobs2 = range(nobs2)
        if im==12:
            print '   -Using all (equal weight of data types)...'
            Jsort = Jstream + Jsatext + JSWC/5 + JT/2
            Jprior2 = 4
            iobs2 = range(nobs2)


        Jsort_max = np.sort(np.asanyarray(Jsort))[nbest-1]    
        Jsort_min = np.sort(np.asanyarray(Jsort))[0]
        print Jsort_min,Jsort_max,Jprior2

        Jselected = [Jsort[ix] for ix in range(lok) if Jsort[ix]<=Jsort_max]
        iok3 = [iok[ix] for ix in range(lok) if (Jsort[ix]<=Jsort_max and iok[ix]==iok2[ix])]
        jok3 = [jok[ix] for ix in range(lok) if (Jsort[ix]<=Jsort_max and jok[ix]==jok2[ix])]
        sok3 = [sok[ix] for ix in range(lok) if (Jsort[ix]<=Jsort_max and sok[ix]==sok2[ix])]

        print iok3
        print jok3
        print sok3

        for iobs in iobs2:
            print obsnames2[iobs],
            f_Jout = os.getcwd()+'/'+outdir+'/'+MCname+'_CF'+metric[im]+'.'+str(nbest)+'best_'+obsnames2[iobs]+'_Ts.txt'
            with open(f_Jout,'w') as f_out:
                f_out.write('J,'+','.join([str(it+1) for it in range(lsim)])+'\n')

            for i in range(nbest):
                # Get all simulation of this job
                if i==0 or iok3[i]!=iok3[i-1]:
                    f_in = simdir+str(iok3[i])+'/'+simfiles2[iobs]
                    js = np.genfromtxt(f_in,delimiter=',',skip_header=1,unpack=True,dtype=np.int)[0]
                    tmp = np.genfromtxt(f_in,delimiter=',',skip_header=1)*simfct[iobs]

                sim = [tmp[jok3[i]-1][x] for x in range(1,lsim+1)]
                # Add to list
                with open(f_Jout,'a') as f_out:
                    f_out.write(str(Jselected[i])+','+','.join([str(val) for val in sim])+'\n')
        print

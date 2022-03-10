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
swpar = 1
swsim = 0
nbest = 30

# -- Output directory
outdir = 'Outputs'

# --- Definition ---
# Path to data
obsdir = '/users/s01ks6/Data/BB'
# Simulation subdirectories
MCname = 'MC8'
# Observations
obsnames = ['Streamflow','SaturationArea','SWC_Peat','SWC_Gley','SWC_Podzol',
            #'SWC_NFforest',
            'SWC_SFforest','T_NFforest','T_SFforest',
            'NetRad_WS1','NetRad_WS2','NetRad_WS3']#,'SWC_forest','T_forest']
nobs = len(obsnames)

# Corresponding simulations outputs
simdir = os.getcwd()+'/'+MCname+'_sampling.'
simfiles = ['Streamflow_hist.tab','SaturationArea_hist.tab',
            'SWC_Peat_hist.tab','SWC_Gley_hist.tab','SWC_Podzol_hist.tab',
            #'SWC_NFforest_hist.tab',
            'SWC_SFforest_hist.tab','T_NFforest_hist.tab','T_SFforest_hist.tab',
            'NetRad_WS1_hist.tab','NetRad_WS2_hist.tab','NetRad_WS3_hist.tab']
sim_order = [6,1,8,2,5,3,7,4,9]
# Length of the simulated outputs
lsim = 1941
simt = [datetime(2011,6,1)+timedelta(days=x) for x in range(lsim)]

#prdir = os.getcwd()+'/Results_prior'
#prfiles = ['Streamflow.tab','BasinSummary.txt',
#           'SoilMoistureAv.tab','SoilMoistureAv.tab',#'SoilMoistureAv.tab',
#           'SoilMoistureAv.tab','SoilMoistureAv.tab','Transpiration[0].tab','Transpiration[0].tab',
#           'NetRadToC.tab','NetRadToC.tab','NetRadToC.tab']
#prcol = [1,12,2,3,4,5,6,5,6,7,8,9]
# Number of locations for simulated time series
#nts = 9#len(prcol)
# Unit conversion sim --> obs
simfct = [1,1,1,1,1,1,#1,
          8.64e7,8.64e7,1,1,1] # sim T is in m/s --> mm/d
# Obs filenames
obsfiles = ['/Q/BB_discharge_daily_01062011-20092016.csv',
            '/SatExt/SatExt_daily_01062011-30092014.csv',
            '/SWC/VSM_Site1_daily_01062011-07092016.csv',
            '/SWC/VSM_Site2_daily_20042011-27092016.csv',
            '/SWC/VSM_Site3_daily_21042011-14102016.csv',
            '/SWC/VSM_SFForest_daily_24022015-05102016.csv',
            #'/SWC/VSM_Plantation_daily_10102012-14102016.csv',
            '/Transpiration/T_NFforest_08072015-28092015.csv',
            '/Transpiration/T_SFforest_01042016-22092016.csv',
            '/NetRad/station1_chickencage_daily_17072014-08082016.csv',
            '/NetRad/station2_bog_daily_17072014-03082016.csv',
            '/NetRad/station3_hilltop_daily_17042015-12072016.csv']
#obsbeg = [datetime(2011,6,1),datetime(2011,6,1),datetime(2011,4,20),
#          datetime(2011,4,21),datetime(2015,2,24),datetime(2015,7,8)]
#obsend = [datetime(2016,7,4),datetime(2016,9,9),datetime(2016,9,27),
#          datetime(2016,10,14),datetime(2016,10,5),datetime(2015,9,28)]
obscol = [1,1,1,4,4,4,#4,
          1,1,1,1,1]
# Number of samples
nit = 1000
# There is a bug in the outputs (for now), so we skip a few runs...
rngit = range(nit)#range(27,nit)
# Number of parallel runs
njob = 100

# Time frames used for the metrics
fitbeg = [datetime(2012,11,1),datetime(2012,11,1),datetime(2012,11,1),
          datetime(2012,11,1),datetime(2012,11,1),datetime(2015,2,24),
          #datetime(2012,11,1),
          datetime(2015,7,8),datetime(2016,4,1),
          datetime(2014,7,17),datetime(2014,7,17),datetime(2015,4,17)]
fitend = [datetime(2014,10,31),datetime(2013,10,31),datetime(2014,10,31),
          datetime(2014,10,31),datetime(2014,10,31),datetime(2016,2,23),
          #datetime(2014,10,31),
          datetime(2015,9,28),datetime(2016,9,22),
          datetime(2015,7,16),datetime(2015,7,16),datetime(2016,4,15)]
# length (in years) of the fit period (used for the total cost function)
lfit = [730,365,730,730,730,365,#730,
        83,175,365,365,365]
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
sigma2 = {}
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
    #if obsnames[iobs] == 'SaturationArea':
    #    hskip=1
    #    idx = prcol[iobs]-1
    #else:
    #    hskip= nts+3
    #    idx = np.argsort(np.array(sim_order))[prcol[iobs]-1]+1

    #tmp = (np.genfromtxt(prdir+'/'+prfiles[iobs],delimiter='\t',
    # skip_header=hskip,unpack=True)[idx])*simfct[iobs]

    #prior = [tmp[idx] for idx in range(lsim) if simt[idx]>=fitbeg[iobs] and simt[idx]<=fitend[iobs]]
    #sigma[obsnames[iobs]] = metrics.mse(prior,obs[obsnames[iobs]])

    sigma2[obsnames[iobs]] = np.nanmax(obs[obsnames[iobs]]) - np.nanmin(obs[obsnames[iobs]])

    print obsnames[iobs], len(obs[obsnames[iobs]])#, len(prior)
#print sigma#[str(sigma[obs]) for obs in obsnames]
#print pr_nse
print sigma2

#tmp = ','.join(['CF_'+x for x in obsnames])
#with open(os.getcwd()+'/'+outdir+'/'+MCname+'_CFmse_prior.txt','w') as f_out:
#    f_out.write(tmp+'\n')
#    f_out.write(','.join([str(sigma[x]) for x in obsnames])+'\n')

tmp = ','.join([x for x in obsnames])
with open(os.getcwd()+'/'+outdir+'/'+MCname+'_Obs_amp.txt','w') as f_out:
    f_out.write(tmp+'\n')
    f_out.write(','.join([str(sigma2[x]) for x in obsnames])+'\n')


#sys.exit()
# ------------------------------------------------------------------------------------

if switch == 0:

    print 'Get posterior misfits...'
    # Loop over jobs
    tmp = ','.join(['CF_'+x for x in obsnames])
    with open(os.getcwd()+'/'+outdir+'/'+MCname+'_CFmse.txt','w') as f_out:
        f_out.write('Sample,'+tmp+'\n')
    with open(os.getcwd()+'/'+outdir+'/'+MCname+'_parameters.txt','w') as f_out:
        f_in = simdir+'1/Parameters.txt'
        tmp = list(np.genfromtxt(f_in,delimiter=',',dtype= '|S20')[0])
        npar = len(tmp)-1
        print npar
        pnames = ','.join([str(tmp[idx]) for idx in range(1,npar+1)])
        f_out.write('Job,Iteration,Sample,'+pnames+'\n')

    # initialize costfunction
    Jmse = {}

    itot = 1
    itot2 = 1

    for i in range(1,njob+1):
        
        print 'Job', i,'...',

        for iobs in range(nobs):
            f_sim = simdir+str(i)+'/'+simfiles[iobs]
            tmp = np.genfromtxt(f_sim,delimiter=',',skip_header=1)

            if iobs==0:
                # Get the number of runs that worked
                js = np.genfromtxt(f_sim,delimiter=',',skip_header=1,unpack=True,dtype=np.int)[0]
                js = [js[idx] for idx in range(len(js)) if js[idx]<=nit]
                nj = len(js)
                print 'nj =',nj,':'

                f_par = simdir+str(i)+'/Parameters.txt'
                tmp_par = np.genfromtxt(f_par,delimiter=',',skip_header=1)

            # initialize costfunction
            Jmse[obsnames[iobs]] = [0]*nj

            for j in range(1,nj+1):

                # -- Cost functions
                tmp2 = tmp[j-1]*simfct[iobs]

                # Crop between desired time frame 
                sim = [tmp2[idx+1] for idx in range(lsim) if simt[idx]>=fitbeg[iobs] and simt[idx]<=fitend[iobs]]
                # Increment cost function              
                Jmse[obsnames[iobs]][j-1] = metrics.sse(sim,obs[obsnames[iobs]])#/sigma[obsnames[iobs]]
                if j==1:
                    print np.mean(sim), np.mean(np.ma.masked_array(obs[obsnames[iobs]],np.isnan(obs[obsnames[iobs]]))), len(sim), len(obs[obsnames[iobs]]), Jmse[obsnames[iobs]][0],
                
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

    print 'Number of complete samples:',itot, itot2
    print

#############################################################################################
# Time series of the best runs
# ----------------------------

if switch == 1:

    # Metrics used
    metric = ['Streamflow','SatExt',
              'SWC_Peat','SWC_Gley','SWC_Podzol',#'SWC_NF',
              'SWC_SF','SWCall',
              'T_NF','T_SF','Tall','NR_WS1','NR_WS2','NR_WS3','NRall',
              'Total']#,'TotalLgW']
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
    f_amp = os.getcwd()+'/'+outdir+'/'+MCname+'_Obs_amp.txt'
    Oamp = np.genfromtxt(f_amp,delimiter=',',skip_header=1,unpack=True)[:nobs]

    # Using streamflow only 
    Jstream = np.genfromtxt(f_in,delimiter=',',skip_header=1,unpack=True)[3]
    Jstream = Jstream/(Oamp[0]**2)
    # using saturation extent
    Jsatext = np.genfromtxt(f_in,delimiter=',',skip_header=1,unpack=True)[4]
    Jsatext = Jsatext/(Oamp[1]**2)
    # using swc
    JSWCs = np.genfromtxt(f_in,delimiter=',',skip_header=1,unpack=True)[5:5+4]
    JSWCs = np.transpose(np.transpose(JSWCs)/Oamp[[2,3,4,5]]**2)
    JSWC = np.sum(JSWCs,axis=0)#JSWCs[0]/(Oamp[2]**2)+JSWCs[1]/(Oamp[3]**2)+JSWCs[2]/(Oamp[4]**2)+JSWCs[3]/(Oamp[5]**2)+JSWCs[4]/(Oamp[6]**2)
    # using transpiration
    JTs = np.genfromtxt(f_in,delimiter=',',skip_header=1,unpack=True)[[9,10]]
    JTs = np.transpose(np.transpose(JTs)/Oamp[[6,7]]**2)
    #JTSF = np.genfromtxt(f_in,delimiter=',',skip_header=1,unpack=True)[11]/Oamp[8]**2
    JT = np.sum(JTs,axis=0)#JTNF/(Oamp[7]**2) + JTSF/(Oamp[8]**2)
    # using transpiration
    JNRs = np.genfromtxt(f_in,delimiter=',',skip_header=1,unpack=True)[[11,12,13]]
    JNRs = np.transpose(np.transpose(JNRs)/Oamp[[8,9,10]]**2)
    #JNR1 = np.genfromtxt(f_in,delimiter=',',skip_header=1,unpack=True)[10]
    #JNR2 = np.genfromtxt(f_in,delimiter=',',skip_header=1,unpack=True)[11]
    #JNR3 = np.genfromtxt(f_in,delimiter=',',skip_header=1,unpack=True)[12]
    JNR = np.sum(JNRs,axis=0)#JNR1/(Oamp[9]**2) + JNR2/(Oamp[10]**2) + JNR3/(Oamp[11]**2)
    #print Jmse
    #print Jns

    # Best runs selections -----------------------------------------
    print 'Fetch',nbest,' best runs'
    
    for im in range(0,1):#,nmet):
        # Using the MSE J
        if im==0:
            print '   -Using streamflow...'
            Jsort = Jstream
            iobs2 = [0]
        if im==1:
            print '   -Using saturation extent...'
            Jsort = Jsatext
            iobs2 = [1]
        if im==2:
            print '   -Using soil moisture (peat)...'
            Jsort = JSWCs[0]
            iobs2 = [2,3,4,5]
        if im==3:
            print '   -Using soil moisture (gley)...'
            Jsort = JSWCs[1]
            iobs2 = [2,3,4,5]
        if im==4:
            print '   -Using soil moisture (podzol)...'
            Jsort = JSWCs[2]
            iobs2 = [2,3,4,5]
        #if im==5:
        #    print '   -Using soil moisture (NF forest)...'
        #    Jsort = JSWCs[3]
        #    iobs2 = [2,3,4,5,6]
        if im==5:
            print '   -Using soil moisture (SF forest)...'
            Jsort = JSWCs[3]
            iobs2 = [2,3,4,5]
        if im==6:
            print '   -Using soil moisture (all)...'
            Jsort = JSWC
            iobs2 = [2,3,4,5]
        if im==7:
            print '   -Using transpiration (NF forest)...'
            Jsort = JTs[0]
            iobs2 = [6,7]
        if im==8:
            print '   -Using transpiration (SF forest)...'
            Jsort = JTs[1]
            iobs2 = [6,7]
        if im==9:
            print '   -Using transpiration (all)...'
            Jsort = JT
            iobs2 = [6,7]
        if im==10:
            print '   -Using net radiation (WS1)...'
            Jsort = JNRs[0]
            iobs2 = [8,9,10]
        if im==11:
            print '   -Using net radiation (WS2)...'
            Jsort = JNRs[1]
            iobs2 = [8,9,10]
        if im==12:
            print '   -Using net radiation (WS3)...'
            Jsort = JNRs[2]
            iobs2 = [8,9,10]
        if im==13:
            print '   -Using net radiation (all)...'
            Jsort = JNR
            iobs2 = [8,9,10]
        if im==14:
            print '   -Using all ...'#(equal weight of data streams)...'
            Jsort = Jstream + Jsatext + JSWC + JT + JNR
            iobs2 = range(nobs)
        #if im==16:
        #    print '   -Using all (weighted by information content = record length)...'
        #    Jsort = np.sum(np.transpose(np.concatenate([[Jstream],[Jsatext],
        #                                                [JSWCs[0]],[JSWCs[1]],[JSWCs[2]],[JSWCs[3]],[JSWCs[4]],
        #                                                [JTs[0]],[JTs[1]],[JNRs[0]],[JNRs[1]],[JNRs[2]]],axis=0))*lfit,axis=1)
        #    Jprior2 = 4
        #    iobs2 = range(nobs)


        Jsort_max = np.sort(np.asanyarray(Jsort))[nbest-1]    
        Jsort_min = np.sort(np.asanyarray(Jsort))[0]
        print Jsort_min,Jsort_max

        Jselected = [Jsort[ix] for ix in range(lok) if Jsort[ix]<=Jsort_max]
        iok3 = [iok[ix] for ix in range(lok) if (Jsort[ix]<=Jsort_max and iok[ix]==iok2[ix])]
        jok3 = [jok[ix] for ix in range(lok) if (Jsort[ix]<=Jsort_max and jok[ix]==jok2[ix])]
        sok3 = [sok[ix] for ix in range(lok) if (Jsort[ix]<=Jsort_max and sok[ix]==sok2[ix])]

        print iok3
        print jok3
        print sok3

        # Parameters
        if swpar == 1:
            f_in = simdir+'1/Parameters.txt'
            tmp = list(np.genfromtxt(f_in,delimiter=',',dtype= '|S20')[0])
            npar = len(tmp)-1
            pnames = ','.join([str(tmp[idx]) for idx in range(1,npar+1)])

            f_Jout = os.getcwd()+'/'+outdir+'/'+MCname+'_CF'+metric[im]+'.'+str(nbest)+'best_paramsets.txt'
            with open(f_Jout,'w') as f_out:
                f_out.write('J,'+pnames+'\n')

            for i in range(nbest):
                # Get all simulation of this job
                if i==0 or iok3[i]!=iok3[i-1]:
                    f_in = simdir+str(iok3[i])+'/Parameters.txt'
                    tmp = np.genfromtxt(f_in,delimiter=',',skip_header=1)

                par = [tmp[jok3[i]-1][x] for x in range(1,npar+1)]
                # Add to list
                with open(f_Jout,'a') as f_out:
                    f_out.write(str(Jselected[i])+','+','.join([str(val) for val in par])+'\n')
        
        if swsim == 1:
            # Time series
            for iobs in iobs2:
                print obsnames[iobs],
                f_Jout = os.getcwd()+'/'+outdir+'/'+MCname+'_CF'+metric[im]+'.'+str(nbest)+'best_'+obsnames[iobs]+'_Ts.txt'
                with open(f_Jout,'w') as f_out:
                    f_out.write('J,'+','.join([str(it+1) for it in range(lsim)])+'\n')

                for i in range(nbest):
                    # Get all simulation of this job
                    if i==0 or iok3[i]!=iok3[i-1]:
                        f_in = simdir+str(iok3[i])+'/'+simfiles[iobs]
                        js = np.genfromtxt(f_in,delimiter=',',skip_header=1,unpack=True,dtype=np.int)[0]
                        tmp = np.genfromtxt(f_in,delimiter=',',skip_header=1)*simfct[iobs]

                    sim = [tmp[jok3[i]-1][x] for x in range(1,lsim+1)]
                    # Add to list
                    with open(f_Jout,'a') as f_out:
                        f_out.write(str(Jselected[i])+','+','.join([str(val) for val in sim])+'\n')
        print

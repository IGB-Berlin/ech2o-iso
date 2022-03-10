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
from optparse import OptionParser

# --- Subroutine(s)
import eval_measures as metrics

parser = OptionParser()

# Switch : fetching obs metrics or corresponding parameters / best runs ?
# 0 = metrics + par, 1 = best runs
parser.add_option("--switch",dest="switch",metavar="SWITCH")
switch = int(options.switch)

if switch==1:
    parser.add_option("--nbest",dest="nbest",metavar="NBEST")
    parser.add_option("--swpar",dest="swpar",metavar="SWPAR")
    parser.add_option("--swsim",dest="swsim",metavar="SWSIM")

(options, args) = parser.parse_args()

if switch==1:
    nbest = int(options.nbest)
    swpar = int(options.swpar)
    swsim = int(options.swsim)

#switch = 1
#swpar = 1
#swsim = 1
#nbest = 30

# -- Output directory
outdir = 'Outputs'

# --- Definition ---
# Path to data
obsdir = '/users/s01ks6/Data/BB'
# Simulation subdirectories
MCname = 'MC9'
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
osig = {}
oave = {}
nok = 0

if switch == 0:

    print 'Get observations...'
    print

    for iobs in range(nobs):

        # -- Get the obs
        f_obs = obsdir+obsfiles[iobs]
        tmp = np.genfromtxt(f_obs,delimiter=',',skip_header=1,unpack=True, dtype= '|S10')[0]
        obst = [datetime.strptime(a, '%d/%m/%Y') for a in tmp]
        lobs = len(obst)
        tmp = np.genfromtxt(f_obs,delimiter=',',skip_header=1,unpack=True)[obscol[iobs]]
        # Crop obs between desired time frame
        obs[obsnames[iobs]] = [tmp[idx] for idx in range(lobs) if obst[idx]>=fitbeg[iobs] and obst[idx]<=fitend[iobs]]


    print 'Get sample misfits...'
    # -- Headers
    # Metrics
    with open(os.getcwd()+'/'+outdir+'/'+MCname+'_KGEmod.txt','w') as f_out:
        f_out.write('Job,Iteration,Sample,'+','.join([x for x in obsnames])+'\n')
    with open(os.getcwd()+'/'+outdir+'/'+MCname+'_KGEcorr.txt','w') as f_out:
        f_out.write('Job,Iteration,Sample,'+','.join([x for x in obsnames])+'\n')
    with open(os.getcwd()+'/'+outdir+'/'+MCname+'_KGEbeta.txt','w') as f_out:
        f_out.write('Job,Iteration,Sample,'+','.join([x for x in obsnames])+'\n')
    with open(os.getcwd()+'/'+outdir+'/'+MCname+'_KGEgamma.txt','w') as f_out:
        f_out.write('Job,Iteration,Sample,'+','.join([x for x in obsnames])+'\n')
    # Parameters
    with open(os.getcwd()+'/'+outdir+'/'+MCname+'_parameters.txt','w') as f_out:
        f_in = simdir+'1/Parameters.txt'
        tmp = list(np.genfromtxt(f_in,delimiter=',',dtype= '|S20')[0])
        npar = len(tmp)-1
        print npar
        pnames = ','.join([str(tmp[idx]) for idx in range(1,npar+1)])
        f_out.write('Job,Iteration,Sample,'+pnames+'\n')

    # initialize costfunction
    KGEmod = {}
    KGEcorr = {}
    KGEbeta = {}    
    KGEgamma = {}

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
            KGEmod[obsnames[iobs]] = [0]*nj
            KGEcorr[obsnames[iobs]] = [0]*nj
            KGEgamma[obsnames[iobs]] = [0]*nj
            KGEbeta[obsnames[iobs]] = [0]*nj

            for j in range(1,nj+1):

                # -- Cost functions
                tmp2 = tmp[j-1]*simfct[iobs]

                # Crop between desired time frame 
                sim = [tmp2[idx+1] for idx in range(lsim) if simt[idx]>=fitbeg[iobs] and simt[idx]<=fitend[iobs]]
                # Increment cost function              
                KGEcorr[obsnames[iobs]][j-1] = metrics.corr(sim,obs[obsnames[iobs]])
                KGEbeta[obsnames[iobs]][j-1] = metrics.rmu(sim,obs[obsnames[iobs]])
                KGEgamma[obsnames[iobs]][j-1] = metrics.rstd(sim,obs[obsnames[iobs]])/KGEbeta[obsnames[iobs]][j-1]
                KGEmod[obsnames[iobs]][j-1] = 1- np.sqrt((1 - KGEcorr[obsnames[iobs]][j-1])**2
                                                         + (1 - KGEgamma[obsnames[iobs]][j-1])**2
                                                         + (1 - KGEbeta[obsnames[iobs]][j-1])**2)

                # A few prints
                if j==1:
                    print obsnames[iobs]
                    print np.mean(sim), np.mean(np.ma.masked_array(obs[obsnames[iobs]],np.isnan(obs[obsnames[iobs]]))), len(sim), len(obs[obsnames[iobs]])
                    print KGEcorr[obsnames[iobs]][0], KGEgamma[obsnames[iobs]][0], KGEbeta[obsnames[iobs]][0], KGEmod[obsnames[iobs]][0]
                
                if iobs ==0:
                    # -- Parameters 
                    tmp3 = [tmp_par[j-1][idx] for idx in range(1,npar+1)]
                    with open(os.getcwd()+'/'+outdir+'/'+MCname+'_parameters.txt','a') as f_out:
                        f_out.write(str(i)+','+str(j)+','+str(itot)+','+','.join([str(x) for x in tmp3])+'\n')
                    itot+=1
        print       

        # Write
        with open(os.getcwd()+'/'+outdir+'/'+MCname+'_KGEmod.txt','a') as f_out:
            for j in range(1,nj+1):
                f_out.write(','.join([str(i),str(j),str(itot2),','.join([str(KGEmod[x][j-1]) for x in obsnames])])+'\n')
                itot2+=1

        with open(os.getcwd()+'/'+outdir+'/'+MCname+'_KGEcorr.txt','a') as f_out:
            for j in range(1,nj+1):
                f_out.write(','.join([str(i),str(j),str(itot2-nj+j),','.join([str(KGEcorr[x][j-1]) for x in obsnames])])+'\n')

        with open(os.getcwd()+'/'+outdir+'/'+MCname+'_KGEgamma.txt','a') as f_out:
            for j in range(1,nj+1):
                f_out.write(','.join([str(i),str(j),str(itot2-nj+j),','.join([str(KGEgamma[x][j-1]) for x in obsnames])])+'\n')

        with open(os.getcwd()+'/'+outdir+'/'+MCname+'_KGEbeta.txt','a') as f_out:
            for j in range(1,nj+1):
                f_out.write(','.join([str(i),str(j),str(itot2-nj+j),','.join([str(KGEbeta[x][j-1]) for x in obsnames])])+'\n')


    print 'Number of complete samples:',itot, itot2
    print

#############################################################################################
# Time series of the best runs
# ----------------------------

if switch == 1:

    # Metrics used


    # Locate indexes
    f_in = os.getcwd()+'/'+outdir+'/'+MCname+'_'+str(nbest)+'best.jobs.txt'
    iok = np.genfromtxt(f_in,delimiter=',',skip_header=1,unpack=True,dtype=int)
    metric = list(np.genfromtxt(f_in,delimiter=',',dtype= '|S20')[0])
    nmet = len(metric)

    f_in = os.getcwd()+'/'+outdir+'/'+MCname+'_'+str(nbest)+'best.iter.txt'
    jok = np.genfromtxt(f_in,delimiter=',',skip_header=1,unpack=True,dtype=int)
    f_in = os.getcwd()+'/'+outdir+'/'+MCname+'_'+str(nbest)+'best.sample.txt'
    sok = np.genfromtxt(f_in,delimiter=',',skip_header=1,unpack=True,dtype=int)
    
    print 'Fetch',nbest,' best runs'

    indm = [0,1,2,3,4,5,6,7,8,9,10,11,12,23]
    #indm = [23]
    
    for im in indm:
        print '   -Using '+metric[im]+'...'

        iok2 = iok[im]
        jok2 = jok[im]

        # Parameters
        if swpar == 1:
            print 'parameters...',
            f_in = simdir+'1/Parameters.txt'
            tmp = list(np.genfromtxt(f_in,delimiter=',',dtype= '|S20')[0])
            npar = len(tmp)-1
            pnames = ','.join([str(tmp[idx]) for idx in range(1,npar+1)])

            f_Jout = os.getcwd()+'/'+outdir+'/'+MCname+'_KGEortho-'+metric[im]+'.'+str(nbest)+'bestParams.txt'
            with open(f_Jout,'w') as f_out:
                f_out.write(pnames+'\n')

            if type(iok2)!=np.int64:
                for i in range(nbest):
                    # Get all simulation of this job
                    if i==0 or iok2[i]!=iok2[i-1]:
                        f_in = simdir+str(iok2[i])+'/Parameters.txt'
                        tmp = np.genfromtxt(f_in,delimiter=',',skip_header=1)

                    par = [tmp[jok2[i]-1][x] for x in range(1,npar+1)]
                    # Add to list
                    with open(f_Jout,'a') as f_out:
                        f_out.write(','.join([str(val) for val in par])+'\n')
            
            else:
                 f_in = simdir+str(iok2)+'/Parameters.txt'
                 tmp = np.genfromtxt(f_in,delimiter=',',skip_header=1)
                 par = [tmp[jok2-1][x] for x in range(1,npar+1)]
                 # Add to list
                 with open(f_Jout,'a') as f_out:
                    f_out.write(','.join([str(val) for val in par])+'\n')
        
        if swsim == 1:
            print 'time series...',
            # Time series
            f_Jout = os.getcwd()+'/'+outdir+'/'+MCname+'_KGEortho-'+metric[im]+'.'+str(nbest)+'bestTs.txt'
            with open(f_Jout,'w') as f_out:
                f_out.write('Ts,'+','.join([str(it+1) for it in range(lsim)])+'\n')
            for iobs in range(nobs):
                print obsnames[iobs],

                if type(iok2)!=np.int64:
                    for i in range(nbest):
                        # Get all simulation of this job
                        if i==0 or iok2[i]!=iok2[i-1]:
                            f_in = simdir+str(iok2[i])+'/'+simfiles[iobs]
                            tmp = np.genfromtxt(f_in,delimiter=',',skip_header=1)
                        sim = [tmp[jok2[i]-1][x] for x in range(1,lsim+1)]

                    #f_in = simdir+str(iok2[i])+'/'+simfiles[iobs]
                    #tmp = np.genfromtxt(f_in,delimiter=',',skip_header=1)[jok2[i]]*simfct[iobs]

                        # Add to list
                        with open(f_Jout,'a') as f_out:
                            f_out.write(obsnames[iobs]+','+','.join([str(val*simfct[iobs]) for val in sim])+'\n')
                
                else:
                    # Get all simulation of this job
                    f_in = simdir+str(iok2)+'/'+simfiles[iobs]
                    tmp = np.genfromtxt(f_in,delimiter=',',skip_header=1)
                    sim = [tmp[jok2-1][x] for x in range(1,lsim+1)]
                    # Add to list
                    with open(f_Jout,'a') as f_out:
                        f_out.write(obsnames[iobs]+','+','.join([str(val*simfct[iobs]) for val in sim])+'\n')
        print

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
parser.add_option("--nbest",dest="nbest",metavar="NBEST")
parser.add_option("--metric",dest="metric",metavar="metric")
parser.add_option("--swpar",dest="swpar",metavar="SWPAR")
parser.add_option("--swsim",dest="swsim",metavar="SWSIM")
parser.add_option("--ext",dest="ext",metavar="EXT")

(options, args) = parser.parse_args()

switch = int(options.switch)
MCname = copy.copy(options.ext)

if switch==1:
    nbest = int(options.nbest)
    namemet = copy.copy(options.metric)
    swpar = int(options.swpar)
    swsim = int(options.swsim)

#switch = 1
#swpar = 1
#swsim = 1
#nbest = 30

# -- Output directory
outdir = 'Outputs'

# -- Date set-up
tstart = datetime(2013,2,21) # start of 'useful' sim
leff = 1265 # useful sim length
spinup = 1095 # spinup length
lsim = leff + spinup # total sim length
simt = [tstart-timedelta(days=spinup)+timedelta(days=x) for x in range(lsim)]

# --- Definition ---
# Path to data
obsdir = '/users/s01ks6/Data/BB'
# Simulation subdirectories
#MCname = 'MC12'
# Observations
obsnames = ['Streamflow',
            'SWC_Peat','SWC_Gley','SWC_Podzol','SWC_SForest','SWC_NForest',
            'T_SForest','T_NForest','T_SHeather','T_NHeather','ET_SHeather','ET_NHeather',
            'NetRad_WS1','NetRad_WS2','NetRad_WS3',
            'dD_Stream','d18O_Stream',
            'dD_DW1','d18O_DW1','GWD_DW1',
            'dD_DW2','d18O_DW2','GWD_DW2',
            'dD_DW3','d18O_DW3','GWD_DW3',
            'dD_DW4','d18O_DW4','GWD_DW4']

#'SaturationArea']#,'SWC_forest','T_forest']
nobs = len(obsnames)

# Corresponding simulations outputs
simdir = os.getcwd()+'/'+MCname+'_sampling.'
simfiles = ['Streamflow_hist.tab',
            'SWC_Peat_hist.tab','SWC_Gley_hist.tab','SWC_Podzol_hist.tab','SWC_SForest_hist.tab','SWC_NForest_hist.tab',
            'T_SForest_hist.tab','T_NForest_hist.tab','T_SHeather_hist.tab','T_NHeather_hist.tab','ET_SHeather_hist.tab','ET_NHeather_hist.tab',
            'NetRad_WS1_hist.tab','NetRad_WS2_hist.tab','NetRad_WS3_hist.tab',
            'dD_Stream_hist.tab','d18O_Stream_hist.tab',
            'dD_DW1_hist.tab','d18O_DW1_hist.tab','GWD_DW1_hist.tab',
            'dD_DW2_hist.tab','d18O_DW2_hist.tab','GWD_DW2_hist.tab',
            'dD_DW3_hist.tab','d18O_DW3_hist.tab','GWD_DW3_hist.tab',
            'dD_DW3_hist.tab','d18O_DW4_hist.tab','GWD_DW4_hist.tab']
            #'SaturationArea_hist.tab']

# Unit conversion sim --> obs
simfct = [1,
          1,1,1,1,1,
          8.64e7,8.64e7,8.64e7,8.64e7,8.64e7,8.64e7,
          1,1,1,
          1,1,
          1,1,1,
          1,1,1,
          1,1,1,
          1,1,1]
#          1] # sim T is in m/s --> mm/d
# Obs filenames
obsfiles = ['/Q/BB_discharge_daily_01062011-01022017.csv',
            '/SWC/VSM_Peat_daily_02062011-06092016.csv',
            '/SWC/VSM_Gley_daily_21042011-26092016.csv',
            '/SWC/VSM_Podzol_daily_22042011-04012017.csv',
            '/SWC/VSM_SForest_daily_11022012-13102016.csv',
            '/SWC/VSM_NForest_daily_25022015-04102016.csv',
            '/ET/T_SForest_08072015-28092015.csv',
            '/ET/T_NForest_01042016-22092016.csv',
            '/ET/ET_T_SHeather_31072015-04082016.csv',
            '/ET/ET_T_NHeather_31072015-22092016.csv',
            '/ET/ET_T_SHeather_31072015-04082016.csv',
            '/ET/ET_T_NHeather_31072015-22092016.csv',
            '/NetRad/station1_chickencage_daily_17072014-08082016.csv',
            '/NetRad/station2_bog_daily_17072014-03082016.csv',
            '/NetRad/station3_hilltop_daily_17042015-12072016.csv',
            '/Isotopes/dD_Stream_01062011-19092016.csv',
            '/Isotopes/d18O_Stream_01062011-19092016.csv',
            '/Groundwater/dD_deeperwells.csv',
            '/Groundwater/d18O_deeperwells.csv',
            '/Groundwater/GWDmetre_DW_daily_09072015-19092016.csv',
            '/Groundwater/dD_deeperwells.csv',
            '/Groundwater/d18O_deeperwells.csv',
            '/Groundwater/GWDmetre_DW_daily_09072015-19092016.csv',
            '/Groundwater/dD_deeperwells.csv',
            '/Groundwater/d18O_deeperwells.csv',
            '/Groundwater/GWDmetre_DW_daily_09072015-19092016.csv',
            '/Groundwater/dD_deeperwells.csv',
            '/Groundwater/d18O_deeperwells.csv',
            '/Groundwater/GWDmetre_DW_daily_09072015-19092016.csv']
            #'/SatExt/SatExt_daily_01062011-30092014.csv']
#obsbeg = [datetime(2011,6,1),datetime(2011,6,1),datetime(2011,4,20),
#          datetime(2011,4,21),datetime(2015,2,24),datetime(2015,7,8)]
#obsend = [datetime(2016,7,4),datetime(2016,9,9),datetime(2016,9,27),
#          datetime(2016,10,14),datetime(2016,10,5),datetime(2015,9,28)]
obscol = [1,#1,
          1,4,4,4,4, #SWC
          1,1,1,1,2,2,# T/ET
          1,1,1, #NR
          1,1, #iso
          1,1,1,1,
          2,2,2,2,
          3,3,3,3,
          4,4,4,4] #GWD
#1]
# Number of samples
nit = 1000
# There is a bug in the outputs (for now), so we skip a few runs...
rngit = range(nit)#range(27,nit)
# Number of parallel runs
njob = 100

# Time frames used for the metrics
fitbeg = [datetime(2013,2,21),
          datetime(2013,2,21),datetime(2013,2,21),datetime(2013,2,21),datetime(2013,2,1),datetime(2015,2,25),
          datetime(2015,7,8),datetime(2016,4,1),datetime(2015,7,31),datetime(2015,7,31),datetime(2015,7,31),datetime(2015,7,31),
          datetime(2014,7,17),datetime(2014,7,17),datetime(2015,4,17),
          datetime(2013,2,21),datetime(2013,2,21),
          datetime(2015,6,9),datetime(2015,6,9),datetime(2015,7,9),
          datetime(2015,6,9),datetime(2015,6,9),datetime(2015,7,9),
          datetime(2015,6,9),datetime(2015,6,9),datetime(2015,7,9),
          datetime(2015,6,9),datetime(2015,6,9),datetime(2015,7,9)]
          
          #datetime(2012,11,1)]
fitend = [datetime(2015,2,20),
          datetime(2015,2,20),datetime(2015,2,20),datetime(2015,2,20),datetime(2015,2,20),datetime(2016,2,23),
          datetime(2015,9,28),datetime(2016,8,8),datetime(2015,10,31),datetime(2015,10,31),datetime(2015,10,31),datetime(2015,10,31),
          datetime(2015,7,16),datetime(2015,7,16),datetime(2016,4,15),
          datetime(2015,2,20),datetime(2015,2,20),
          datetime(2016,8,8),datetime(2016,8,8),datetime(2016,8,8),
          datetime(2016,8,8),datetime(2016,8,8),datetime(2016,8,8),
          datetime(2016,8,8),datetime(2016,8,8),datetime(2016,8,8),
          datetime(2016,8,8),datetime(2016,8,8),datetime(2016,8,8)]
          #datetime(2013,10,31)]

# Metrics
#metrics = ['CF_MSE','CF_NSE']
#nmet = len(metrics

#############################################################################################
# Simulations outputs : saving the metrics
# if switch == 0:

#print
#print 'Summarizing the metrics....'
print

obs={}
obst={}
nok = 0

if switch == 0:

    print 'Get observations...'
    print

    for iobs in range(nobs):
        
        print obsnames[iobs]
        # -- Get the obs
        f_obs = obsdir+obsfiles[iobs]
        tmp = np.genfromtxt(f_obs,delimiter=',',skip_header=1,unpack=True, dtype= '|S10')[0]
        tmpt = np.array([datetime.strptime(a, '%d/%m/%Y') for a in tmp])
        lobs = len(tmpt)
        tmp = np.genfromtxt(f_obs,delimiter=',',skip_header=1,unpack=True)[obscol[iobs]]
        # Crop obs between desired time frame
        obs[obsnames[iobs]] = [tmp[idx] for idx in range(lobs) if tmpt[idx]>=fitbeg[iobs] and tmpt[idx]<=fitend[iobs]]
        obst[obsnames[iobs]] = np.array([tmpt[idx] for idx in range(lobs) if tmpt[idx]>=fitbeg[iobs] and tmpt[idx]<=fitend[iobs]])

        #tmp24 = [obst[idx+1]-obst[idx] for idx in range(lobs) if obst[idx]>=fitbeg[iobs] and obst[idx]<=fitend[iobs]-timedelta(days=1)]
        #print tmp24
                


    print
    print 'Get sample misfits...'
    # -- Headers
    # Metrics
    with open(os.getcwd()+'/'+outdir+'/'+MCname+'_KGE2012.txt','w') as f_out:
        f_out.write('Job,Iteration,Sample,'+','.join([x for x in obsnames])+'\n')
    with open(os.getcwd()+'/'+outdir+'/'+MCname+'_MAE.txt','w') as f_out:
        f_out.write('Job,Iteration,Sample,'+','.join([x for x in obsnames])+'\n')
    with open(os.getcwd()+'/'+outdir+'/'+MCname+'_RMSE.txt','w') as f_out:
        f_out.write('Job,Iteration,Sample,'+','.join([x for x in obsnames])+'\n')
    with open(os.getcwd()+'/'+outdir+'/'+MCname+'_corr.txt','w') as f_out:
        f_out.write('Job,Iteration,Sample,'+','.join([x for x in obsnames])+'\n')

    # Parameters
    with open(os.getcwd()+'/'+outdir+'/'+MCname+'_parameters.txt','w') as f_out:
        f_in = simdir+'21/Parameters.txt'
        tmp = list(np.genfromtxt(f_in,delimiter=',',dtype= '|S20')[0])
        npar = len(tmp)-1
        print
        print 'Length of parameter vector: '+str(npar)
        print
        pnames = ','.join([str(tmp[idx]) for idx in range(1,npar+1)])
        f_out.write('Job,Iteration,Sample,'+pnames+'\n')

    # initialize costfunction
    KGE = {}
    MAE = {}
    RMSE = {}   
    corr = {}
    
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
            KGE[obsnames[iobs]] = [0]*nj
            MAE[obsnames[iobs]] = [0]*nj
            RMSE[obsnames[iobs]] = [0]*nj
            corr[obsnames[iobs]] = [0]*nj
            
            for j in range(1,nj+1):

                # -- Cost functions
                tmp2 = tmp[j-1]*simfct[iobs]

                # Crop between desired time frame, and account for potential gaps in the obs
                sim = [tmp2[idx+1] for idx in range(lsim) if any(obst[obsnames[iobs]]==simt[idx])==True]
                if i==21 and j==1:
                    tmp24 = [simt[idx] for idx in range(lsim) if simt[idx] in obst[obsnames[iobs]]]
                    print
                    print tmp24

                # Increment cost function              
                KGE[obsnames[iobs]][j-1] = metrics.kling_gupta(sim,obs[obsnames[iobs]],method='2012')
                MAE[obsnames[iobs]][j-1] = metrics.meanabs(sim,obs[obsnames[iobs]])
                RMSE[obsnames[iobs]][j-1] = metrics.rmse(sim,obs[obsnames[iobs]])
                corr[obsnames[iobs]][j-1] = metrics.corr(sim,obs[obsnames[iobs]])

                # A few prints
                #if j==1:
                    #print obsnames[iobs]
                    #print np.mean(sim), np.mean(np.ma.masked_array(obs[obsnames[iobs]],np.isnan(obs[obsnames[iobs]]))), len(sim), len(obs[obsnames[iobs]])
                    #print KGE[obsnames[iobs]][0], MAE[obsnames[iobs]][0], RMSE[obsnames[iobs]][0]
                
                if iobs ==0:
                    # -- Parameters 
                    tmp3 = [tmp_par[j-1][idx] for idx in range(1,npar+1)]
                    with open(os.getcwd()+'/'+outdir+'/'+MCname+'_parameters.txt','a') as f_out:
                        f_out.write(str(i)+','+str(j)+','+str(itot)+','+','.join([str(x) for x in tmp3])+'\n')
                    itot+=1
        print       

        # Write
        with open(os.getcwd()+'/'+outdir+'/'+MCname+'_KGE2012.txt','a') as f_out:
            for j in range(1,nj+1):
                f_out.write(','.join([str(i),str(j),str(itot2),','.join([str(KGE[x][j-1]) for x in obsnames])])+'\n')
                itot2+=1

        with open(os.getcwd()+'/'+outdir+'/'+MCname+'_MAE.txt','a') as f_out:
            for j in range(1,nj+1):
                f_out.write(','.join([str(i),str(j),str(itot2-nj+j),','.join([str(MAE[x][j-1]) for x in obsnames])])+'\n')

        with open(os.getcwd()+'/'+outdir+'/'+MCname+'_RMSE.txt','a') as f_out:
            for j in range(1,nj+1):
                f_out.write(','.join([str(i),str(j),str(itot2-nj+j),','.join([str(RMSE[x][j-1]) for x in obsnames])])+'\n')

        with open(os.getcwd()+'/'+outdir+'/'+MCname+'_corr.txt','a') as f_out:
            for j in range(1,nj+1):
                f_out.write(','.join([str(i),str(j),str(itot2-nj+j),','.join([str(corr[x][j-1]) for x in obsnames])])+'\n')

    print 'Number of complete samples:',itot, itot2
    print

#############################################################################################
# Time series of the best runs
# ----------------------------

if switch == 1:

    # Metrics used
    

    # Locate indexes
    f_in = os.getcwd()+'/'+outdir+'/'+MCname+'_'+namemet+'.'+str(nbest)+'best.jobs.txt'
    iok = np.genfromtxt(f_in,delimiter=',',skip_header=1,unpack=True,dtype=int)
    metric = list(np.genfromtxt(f_in,delimiter=',',dtype= '|S20')[0])
    nmet = len(metric)

    f_in = os.getcwd()+'/'+outdir+'/'+MCname+'_'+namemet+'.'+str(nbest)+'best.iter.txt'
    jok = np.genfromtxt(f_in,delimiter=',',skip_header=1,unpack=True,dtype=int)
    f_in = os.getcwd()+'/'+outdir+'/'+MCname+'_'+namemet+'.'+str(nbest)+'best.sample.txt'
    sok = np.genfromtxt(f_in,delimiter=',',skip_header=1,unpack=True,dtype=int)
    
    print 'Fetch',nbest,' best runs'

    indm = [0,1,2,3,4,5,6,7,8,9,10,11,12,23]
    #indm = [23]#10,11,12,23]
    
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

            f_Jout = os.getcwd()+'/'+outdir+'/'+MCname+'_'+namemet+'-'+metric[im]+'.'+str(nbest)+'bestParams.txt'
            with open(f_Jout,'w') as f_out:
                f_out.write(pnames+'\n')

            for i in range(nbest):
                # Get all simulation of this job
                if i==0 or iok2[i]!=iok2[i-1]:
                    f_in = simdir+str(iok2[i])+'/Parameters.txt'
                    tmp = np.genfromtxt(f_in,delimiter=',',skip_header=1)

                par = [tmp[jok2[i]-1][x] for x in range(1,npar+1)]
                # Add to list
                with open(f_Jout,'a') as f_out:
                    f_out.write(','.join([str(val) for val in par])+'\n')
        
        if swsim == 1:
            print 'time series...',
            # Time series
            f_Jout = os.getcwd()+'/'+outdir+'/'+MCname+'_'+namemet+'-'+metric[im]+'.'+str(nbest)+'bestTs.txt'
            with open(f_Jout,'w') as f_out:
                f_out.write('Ts,'+','.join([str(it+1) for it in range(lsim)])+'\n')
            for iobs in range(nobs):
                print obsnames[iobs],
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
        print

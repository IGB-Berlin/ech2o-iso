#!/usr/bin/env python
# -*- coding: utf-8 -*-
# *************************************************
#
# Monte Carlo calibration algorithm for ECH2O
#
# -------
# Routine: Dynamic PCRaster Model for ECH2O
# -------
# Author: S. Kuppel
# Created on 10/2016
# -------------------------------------------------

from pcraster import *
from pcraster.framework import *
from subprocess import *
import numpy as np
import os, time, sys
import csv

# Model syntax in the PCraster dynamic model
class IsoModel(DynamicModel, MonteCarloModel):
    
  def __init__(self, number_of_timesteps,location_of_input_folder, number_of_samples):

    DynamicModel.__init__(self)
    MonteCarloModel.__init__(self)
    
    #input variables to the actual model call
    self.number_of_timesteps = number_of_timesteps
    self.input_folder = location_of_input_folder
    self.number_of_samples=number_of_samples
    
    setclone(self.input_folder+"/maps/mask_new.map")

  def premcloop(self):
    ######################### MAP INPUTS #######################################
    self.ldd = ldd(readmap(self.input_folder+"/maps/ldd_new.map"))              # drainage direction map
    self.dem=readmap(self.input_folder+"/maps/dem_new.map")                      # dem
    self.GRAD=slope(self.dem)                                                   # slope
    self.order = streamorder(self.ldd)                                          # streamorder
    self.upTotal = upstream(self.ldd, scalar(1))                                # upstream cells
    self.twi=readmap(self.input_folder+"/maps/twi.map")                         # topographic wetness index to differentiate between soil types
    self.dsm=readmap(self.input_folder+"/maps/dsm_new.map")                     # surface elevation including vegetation
    self.veg = self.dsm - self.dem	                                        # difference surface - ground elevation, vegetation height
        # self.soil=readmap(self.input_folder+"/maps/Grundlager2.map")              # soil type map, for this application will be derived from twi
    self.RCdsm = self.veg                                                       # PANIC fix to run interception module
    
    # select cells to produce model output: most output for the outlet using pit()
        
    # creation of maps which define the cells with model output are done with command line executables map2col and col2map.
    # in between the output/input texfile is modified in a text editor to have the desired output cell have a value of one, others 0
    # self.outputloc = "Input/maps/sel_top_cell.map"                              # topmost cell domain where snow output plots are created                                                         
    
  def initial(self):
    
    # OUTPUT WRITING 
    # write allways output for discharge and concentration
    self.dischargeTss=TimeoutputTimeseries("aa_discharge", self, pit(self.ldd),noHeader=False)
    self.concQTss=TimeoutputTimeseries("aa_concQ", self, pit(self.ldd),noHeader=False)
    if self.doAges:
        self.ageQTss=TimeoutputTimeseries("ageQ", self, pit(self.ldd),noHeader=False)
    
    # MODEL INITIALISATION
    # all parameter values, initial conditions and model output are generated in script 'initialPara_nse.py'
    initialPara_BB.iniParameters(self)
    
    if doMC:
        
        # PARAMETER PERTUBATIONS
        ### soil
        self.FC1 = numpy.random.uniform(300,1000)                                     # field capacity of wetland 
        self.FC2 = numpy.random.uniform(100,300) 			                #field capacity for forest soil
        self.FC = ifthenelse(self.twi>9, scalar(self.FC1),scalar(self.FC2)) 			#fieldcapacity, total water holding capacity of the soil                                
        self.ks = numpy.random.uniform(0.1,0.5)                               # soil outflow coeff from forest soils
        self.BetaSeepage = numpy.random.uniform(0.1,10)	                        #exponent in soil runoff generation equation    
        self.LP = numpy.random.uniform(0.5,1)
        self.Cflux = numpy.random.uniform(0.1,3)
        ### groundwater
        self.kG = numpy.random.uniform(-4,-2.5)	  	# recession constant baseflow, logarithmic sampling!!!!	[0.0001...~0.003]
        self.kG = 10**self.kG                                              # transform to parameter values again
        self.Ksat = numpy.random.uniform(0.1,10)
        ### snow
        self.deplOffset = numpy.random.uniform(0,-3.5)
        self.Efrac  = numpy.random.uniform(0,15)
        self.albPow = numpy.random.uniform(1,2)
        ### passive storages
        self.GWpas = numpy.random.uniform(0,1000)
        self.SMpas1 = numpy.random.uniform(0,300)
        self.fracSMpas = numpy.random.uniform(0.1,1)
        self.SMpas   = ifthenelse(self.twi>9,scalar(self.SMpas1),scalar(self.fracSMpas)*scalar(self.SMpas1))
            
                    
        ## write parameter values for further processing
        file_name_for_parameters = str(self.currentSampleNumber())+"/parameter.txt"
        parameter_file = open(file_name_for_parameters, "w")
        parameter_values = [float(self.FC1),
                            float(self.FC2),
                            float(self.ks),
                            float(self.BetaSeepage),
                            float(self.LP),
                            float(self.Cflux),
                            float(self.kG),
                            float(self.Ksat), 
                            float(self.deplOffset),
                            float(self.Efrac),
                            float(self.albPow),                  
                            float(self.GWpas),
                            float(self.SMpas1),
                            float(self.fracSMpas),
        ] 
        parameter_file.write('[FC1,FC2,ks,BetaSeepage,LP,Cflux,kG,Ksat,deplOffset,Efrac,albPow,GWpas,SMpas1,fracSMpas]\n') # write the variables above to a file as a header
        parameter_file.write(str(parameter_values))
        parameter_file.close()
    else:
        print 'single run'

  def dynamic(self):
    print 'timestep', self.currentTimeStep(), 'in sample number',self.currentSampleNumber()
    self.setQuiet(quiet=True) # disable timesteps from one
    
    # RUNNING THE MODEL
    # input timeseries from 2011-06-01...2014-10-15
    # specify hydrological input timeseries as model input    
    T = self.readDeterministic(self.input_folder+"/Output_evap/Tcorr")          # temperature [C] corrected to altitude, generated by Marjolein et al
    P = self.readDeterministic(self.input_folder+"/Output_evap/Pcorr")          # precip [mm/d] corrected to altitude, generated by Marjolein et al
    self.ETa= self.readDeterministic(self.input_folder+"/Output_evap/ET02")                    # potential ET, generated by Marjolein et al
    # input data for snowmelt is poor for the period 2003-2004!! see data files
    self.radRs      = self.readDeterministic(self.input_folder+"/output_Rs/RsCorr")            # read mean spatially distributed global radiation (MJ/d*m2) generated by SpatialRs.py
    self.WSin       = timeinputscalar(self.input_folder+"/Output_evap/tsWind2.tss",1)          # read windspeed timeseries, generated by Marjolein et al
    self.RH         = timeinputscalar(self.input_folder+"/Output_evap/tsRHlong2.tss",1)             # read relative humidity timeseries, generated by Marjolein et al
        
    # specify timeseries of precipitation isotope concentration used as model input
    Pconc= timeinputscalar(self.input_folder+"/PconcO18.tss", 1)                   ### FOR NOW a regression from daily temperature OR timeseries generated from monthly means and stds
    
    ##########################
    ###### Snow Module    
    # select between the two: degree-day or energy balance
    snowEbal.update(self, P, T, Pconc)
    #snowModule.update(self, P, T, Pconc)
    
    ##########################
    ###### Interception part    
    interception.update(self)
    
    #########################
    ## soil storage module
    
    soilStorage.update(self, T)
    
    #########################
    #####groundwater storage    
    groundwater.update(self)
    
    #########################
    ## routing module
    routing.update(self)
    
    #################
    #############
    ##Lateral flow in groundwater
    
    latflowGW.update(self)
    
    #######Reporting#########
    reporting_BB.reportAll(self)
    
    self.dischargeTss.sample(self.discharge)    # m3/s
    self.concQTss.sample(self.Qtest2)
    if self.doAges:
        self.ageQTss.sample(self.ageQ2)

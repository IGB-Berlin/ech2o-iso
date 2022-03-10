#/usr/bin/env python
# -*- coding: utf-8 -*-
# *************************************************
#
# Monte Carlo calibration algorithm for ECH2O
#
# -------
# Routine: Main program
# -------
# Author: S. Kuppel
# Created on 10/2016
# -------------------------------------------------

import numpy as np
from optparse import OptionParser
import os, time, sys, glob, copy

# --- Subroutine(s)
import MC_Tools

# ---------
#  OPTIONS
# ---------
parser = OptionParser()

# --- Definition ---
# Configuration file
parser.add_option("--file",dest="file",metavar="FILE",
                  help="Name of the file that defines the calibration to perform")
# Output subdirectory
parser.add_option("--outdir",dest="outdir",metavar="outdir",help="Output directory")

(options, args) = parser.parse_args()

# set working directory
class Config:
    PATH_MAIN = os.getcwd()+'/'

print
print '*********************************************************************************'
print 'Reorganizing outputs - freeing file space: '
print #'- construction of the trajectories'

#############################################################################################
# Initialization of the calibration
# ---------------------------------

# -- General config file (where parameters info is taken from)
if options.file != None:
    [file_py,ext]=os.path.splitext(options.file)
    if len(glob.glob(options.file)) ==0 : sys.exit('# STOP. The file that define the assimilation characteristics does not exist : \n   ...'+file_py+'.py')
    file = options.file
    frep = os.path.dirname(file)
    if frep == '':  file = os.path.join(Config.PATH_MAIN,file)
    print 'The user provided definition file is: '+options.file

# -- Output path
if options.outdir == None:
    options.outdir = os.getcwd()+'/'+os.path.splitext(options.file)[0]

# --- Defining the various PATHs
Config.PATH_OUT  = os.path.abspath(os.path.join(Config.PATH_MAIN,options.outdir))
#Config.PATH_EXEC = os.path.abspath(os.path.join(Config.PATH_OUT,'tmp'))

#-- Import classes and setup from the def file
exec('from ' + file_py + ' import *')

#############################################################################################
# Cleaning loop
# -------------------
for it in range(Opti.nit):
    
    itout = '%03i' % int(it+1)
    Config.PATH_EXEC = Config.PATH_OUT+'/'+itout
    if len(glob.glob(Config.PATH_EXEC))==0:
        break

    print 'Iteration',itout, 'of',Opti.nit

    if len(glob.glob(Config.PATH_EXEC+'/Streamflow.tab'))==0:
        continue
    # Group sampling outputs
    MC_Tools.manage_outputs(Data, Opti, Config, it)
    os.system('rm -fr '+Config.PATH_EXEC)

###############################################################################################
###############################################################################################

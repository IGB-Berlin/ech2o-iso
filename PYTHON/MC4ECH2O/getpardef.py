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
from itertools import chain

# ---------
#  OPTIONS
# ---------
parser = OptionParser()

# --- Definition ---
# Configuration file
parser.add_option("--file",dest="file",metavar="FILE",
                  help="Name of the file that defines the calibration to perform")
# Output subdirectory
parser.add_option("--outext",dest="outext",metavar="outext",help="Output extension")

(options, args) = parser.parse_args()

# set working directory
class Config:
    PATH_MAIN = os.getcwd()+'/'

print
print '*********************************************************************************'
print 'Get parameters setup: min, max, log... '
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

# --- Defining the various PATHs
Config.PATH_OUT  = os.path.abspath(os.path.join(Config.PATH_MAIN,'Outputs'))
Config.EXT_OUT = options.outext

#-- Import classes and setup from the def file
exec('from ' + file_py + ' import *')

#############################################################################################
    
#-- Parameters: from definition to all values
Paras.names = Paras.ref.keys()
Paras.names.sort()
Paras.n = len(Paras.names)

print Paras.names

# Read dictionary to get all param setup
Opti.min = []
Opti.max = []
Opti.log = []
Opti.names = []
for par in Paras.names:
    # Dimensions (soil or veg or 1)
    nr = Paras.ref[par]['soil']*(Site.ns-1)+Paras.ref[par]['veg']*(Site.nv-1)+1
    # Build vectors used in the optimisation
    Opti.min.append(np.repeat(Paras.ref[par]['min'],nr))
    Opti.max.append(np.repeat(Paras.ref[par]['max'],nr))
    Opti.log.append(np.repeat(Paras.ref[par]['log'],nr))
    # For outputs
    if Paras.ref[par]['soil']==1:
        Opti.names.append([par + '_' + s for s in Site.soils])
    elif Paras.ref[par]['veg']==1:
        Opti.names.append([par + '_' + s for s in Site.vegs])
    else:
        Opti.names.append([par])
    
Opti.min = list(chain(*Opti.min))
Opti.max = list(chain(*Opti.max))
Opti.log = list(chain(*Opti.log))
Opti.names = list(chain(*Opti.names))

print Opti.names

# Total number of variables
Opti.nvar = len(Opti.min)

# Outputs the relevant info
f_out = Config.PATH_OUT+'/'+Config.EXT_OUT+'_parameters_char.txt'
with open(f_out,'w') as w_out:
    w_out.write('Names,'+','.join([x for x in Opti.names])+'\n')
    w_out.write('Min,'+','.join([str(x) for x in Opti.min])+'\n')
    w_out.write('Max,'+','.join([str(x) for x in Opti.max])+'\n')
    w_out.write('Log,'+','.join([str(x) for x in Opti.log])+'\n')

###############################################################################################
###############################################################################################

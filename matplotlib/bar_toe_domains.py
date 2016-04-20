#!/bin/env python
# -*- coding: utf-8 -*-
"""
Python matplotlib
Make time emergence PDF bar plots from
 1) Hist vs. histNat runs
 2) using piControl stdev (TODO)

(c) Eric Guilyardi April 2016

TODO: - add arguments for variable and output type

"""
import os, glob
from netCDF4 import Dataset as open_ncfile
import matplotlib.pyplot as plt
from densit_matplot_lib import defVar, averageDom
import numpy as np

# -------------------------------------------------------------------------------
#                                Define work
# -------------------------------------------------------------------------------

inDir = '/Users/ericg/Projets/Density_bining/'
workh = 'Prod_density_april15/toe_histNat'
inDirh = inDir + workh

# Define variable  TODO: read as argument
varname = defVar('salinity')
varname = defVar('temp')
#varname = defVar('depth')
#varname = defVar('volume')
#varname = defVar('persist')
#varname = defVar('heatcontent')

# -------------------------------------------------------------------------------
# file inits

os.chdir(inDirh)
listFiles = glob.glob('cmip5.*')
var = varname['var']

# find dimensions
fi = open_ncfile(inDirh+'/'+listFiles[0])
print inDirh+'/'+listFiles[0]
isond0  = fi.variables['isondepth'] ; # Create variable handle
latN = isond0.shape[3]
levN = isond0.shape[2]
basN = isond0.shape[1]
timN = isond0.shape[0]

# init cumulated number of members
nruns = 0
nrunmax = 200
# init arrays for ToE
toe1 = np.ma.ones([nrunmax,basN,levN,latN], dtype='float32')*1.
# loop over files to fill up array
for i,file in enumerate(listFiles):
    ft    = open_ncfile(inDirh+'/'+file)
    model = file.split('.')[1]
    # read TOE (runN,basN,latN,levN)
    toeread = ft.variables[var+'ToE1'][:]
    # find number of members
    nMembers  = toeread.shape[0]
    nruns1 = nruns + nMembers
    print ' -> Model '+model+' with '+str(nMembers)+' members'
    toe1[nruns:nruns1,...] = toeread
    nruns = nruns1

# crop array
print
print 'Total number of members :',nruns
toe1 = toe1[0:nruns,...]
print toe1.shape

# do bar plot

# select domain

DomToEA1 = {'domain': [-40., -30, 25.5, 26.7], 'name': 'Southern ST', 'color': 'blue'}


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
#varname = defVar('temp')
#varname = defVar('depth')
#varname = defVar('volume')
#varname = defVar('persist')
#varname = defVar('heatcontent')

iniyear = 1861
deltay = 10.

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

levr = fi.variables['lev'][:]
lats = fi.variables['latitude'][:]

# init cumulated number of members
nruns = 0
nrunmax = 200
# init arrays for ToE
toe1,toe2 = [np.ma.ones([nrunmax,basN,levN,latN], dtype='float32')*1. for _ in range(2)]
# loop over files to fill up array
for i,file in enumerate(listFiles):
    ft    = open_ncfile(inDirh+'/'+file)
    model = file.split('.')[1]
    # read TOE (runN,basN,latN,levN)
    toe1read = ft.variables[var+'ToE1'][:]
    toe2read = ft.variables[var+'ToE2'][:]
    # find number of members
    nMembers  = toe1read.shape[0]
    nruns1 = nruns + nMembers
    print ' -> Model '+model+' with '+str(nMembers)+' members'
    toe1[nruns:nruns1,...] = toe1read
    toe2[nruns:nruns1,...] = toe2read
    nruns = nruns1

# crop array
print
print 'Total number of members :',nruns
toe1 = toe1[0:nruns,...]
toe2 = toe2[0:nruns,...]

# do bar plot
# select domain

DomToEA1 = {'domain': [-40., -30, 25.5, 26.7], 'name': 'Southern ST', 'color': 'blue'}
DomToEP1 = {'domain': [-30, -10, 25., 26]   , 'name': 'Southern ST', 'color': 'blue'}
DomToEI1 = {'domain': [-40, -15, 25.6, 26.8] , 'name': 'Southern ST', 'color': 'blue'}

varToEA1 = np.around(averageDom(toe1[:,1,:,:], 3, DomToEA1['domain'], lats, levr) + iniyear)
print min(varToEA1),max(varToEA1), np.around(np.average(varToEA1)), np.median(varToEA1),np.std(varToEA1)
varToEP1 = np.around(averageDom(toe1[:,2,:,:], 3, DomToEP1['domain'], lats, levr) + iniyear)
print min(varToEP1),max(varToEP1), np.around(np.average(varToEP1)), np.median(varToEP1),np.std(varToEP1)
varToEI1 = np.around(averageDom(toe1[:,3,:,:], 3, DomToEI1['domain'], lats, levr) + iniyear)
print min(varToEI1),max(varToEI1), np.around(np.average(varToEI1)), np.median(varToEI1),np.std(varToEI1)

varToE2A1 = np.around(averageDom(toe2[:,1,:,:], 3, DomToEA1['domain'], lats, levr) + iniyear)
print min(varToE2A1),max(varToE2A1), np.around(np.average(varToE2A1)), np.median(varToE2A1),np.std(varToE2A1)


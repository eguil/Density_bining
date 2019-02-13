#!/bin/env python
# -*- coding: utf-8 -*-
"""
Python matplotlib
Make time emergence PDF bar plots from
 1) Hist vs. histNat runs
 2) using piControl stdev (TODO)

(c) Eric Guilyardi April 2016

TODO: - add arguments for variable and output type

TODO: - weigth PDF with 1/N members
TODO: - ToE as function of number of members (use GISS or large ensemble)
TODO: - use piControl for noise rather than histNat as signal is hist-histNat

"""
import os, glob
from netCDF4 import Dataset as open_ncfile
import matplotlib.pyplot as plt
from densit_matplot_lib import defVar, averageDom
import numpy as np

# -------------------------------------------------------------------------------
#                                Define work
# -------------------------------------------------------------------------------

inDir = '/data/ericglod/Density_binning/'
workh = 'Prod_density_april15/toe_histNat'
#workh = 'Prod_density_april15/mme_hist'
inDirh = inDir + workh

# Define variable  TODO: read as argument
varname = defVar('salinity')
#varname = defVar('temp')
#varname = defVar('depth')
#varname = defVar('volume')
#varname = defVar('persist')
#varname = defVar('heatcontent')

iniyear = 1860
finalyear = 2005
deltay = 10.

# -------------------------------------------------------------------------------
# file inits

os.chdir(inDirh)
listFiles = glob.glob('cmip5.*_zon2D.nc')
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
DomToEP1 = {'domain': [-15, -10, 26, 26.3]   , 'name': 'Southern ST', 'color': 'blue'}
DomToEI1 = {'domain': [-40, -15, 25.6, 26.8] , 'name': 'Southern ST', 'color': 'blue'}

varToEA1 = np.around(averageDom(toe2[:,1,:,:], 3, DomToEA1['domain'], lats, levr) + iniyear)
#print min(varToEA1),max(varToEA1), np.around(np.average(varToEA1)), np.median(varToEA1),np.std(varToEA1)
varToEP1 = np.around(averageDom(toe2[:,2,:,:], 3, DomToEP1['domain'], lats, levr) + iniyear)
#print min(varToEP1),max(varToEP1), np.around(np.average(varToEP1)), np.median(varToEP1),np.std(varToEP1)
varToEI1 = np.around(averageDom(toe2[:,3,:,:], 3, DomToEI1['domain'], lats, levr) + iniyear)
#print min(varToEI1),max(varToEI1), np.around(np.average(varToEI1)), np.median(varToEI1),np.std(varToEI1)

varToE2A1 = np.around(averageDom(toe2[:,1,:,:], 3, DomToEA1['domain'], lats, levr) + iniyear)
print min(varToE2A1),max(varToE2A1), np.around(np.average(varToE2A1)), np.median(varToE2A1),np.std(varToE2A1)

ndecades = int((finalyear - iniyear)/deltay)
yearbins = np.arange(ndecades+1)*10+5+iniyear

ToE2A1Means,bins = np.histogram(varToEA1, yearbins)
center = (bins[:-1] + bins[1:]) / 2
ToE2P1Means,bins = np.histogram(varToEP1, yearbins)
ToE2I1Means,bins = np.histogram(varToEI1, yearbins)

fig, ax = plt.subplots(nrows=3, ncols=1)

rects1 = ax[0].bar(center, ToE2A1Means, deltay, color='b')
ax[0].set_ylabel('Members')
ax[0].set_title('ToE1 in '+DomToEA1['name']+' Atl')

rects2 = ax[1].bar(center, ToE2P1Means, deltay, color='b')
ax[1].set_ylabel('Members')
ax[1].set_title('ToE1 in '+DomToEP1['name']+' Pac')

rects3 = ax[2].bar(center, ToE2I1Means, deltay, color='b')
ax[2].set_ylabel('Members')
ax[2].set_title('ToE1 in '+DomToEI1['name']+' Ind')

plt.subplots_adjust(hspace=.5, wspace=5, left=0.04, right=0.86)


plt.show()
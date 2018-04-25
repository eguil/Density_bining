#!/bin/env python
# -*- coding: utf-8 -*-

"""
Python matplotlib
Plot density/latitude Time of Emergence (median and distribution range)

"""

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset as open_ncfile
from maps_matplot_lib import defVarmme, zonal_2D
from modelsDef import defModels
import glob
import os
import colormaps as cmaps

# ----- Workspace ------

indir_toe_rcphn = '/home/ysilvy/Density_bining/Yona_analysis/data/toe_zonal/toe_rcp85_histNat/'

models = defModels()

# ----- Work ------

varname = defVarmme('salinity'); v = 'S'

multStd = 2. # detect ToE at multStd std dev of piControl

use_piC = False # Over projection period, signal = RCP-average(histNat), noise = std(histNat)
# use_piC = True # Over projection period, signal = RCP-average(PiControl), noise = std(PiControl)

iniyear = 1860
finalyear = 2100
deltat = 10.

# density domain
rhomin = 21
rhomid = 26
rhomax = 28
domrho = [rhomin, rhomid, rhomax]

# ----- Variables ------
var = varname['var_zonal_w/bowl']
legVar = varname['legVar']
unit = varname['unit']

# Read latitude and density from original file
fileh_2d = '/data/ericglod/Density_binning/Prod_density_april15/mme_hist/' \
       'cmip5.multimodel_Nat.historical.ensm.an.ocn.Omon.density_zon2D.nc'
fh2d = open_ncfile(fileh_2d, 'r')
lat = fh2d.variables['latitude'][:]; latN = len(lat)
density = fh2d.variables['lev'][:]; levN = len(density)


# ----- Read ToE for each model ------

# == Historical vs historicalNat + RCP8.5 vs. historicalNat or RCP8.5 vs. PiControl ==

nruns = 0 # Initialize total number of runs
nrunmax = 100
nMembers = np.ma.empty(len(models)) # Initialize array for keeping number of members per model

# -- Initialize varToE containing ToE of all runs
varToEA = np.ma.masked_all((nrunmax, levN, latN))
varToEP = np.ma.masked_all((nrunmax, levN, latN))
varToEI = np.ma.masked_all((nrunmax, levN, latN))

# Loop over models
listfiles = glob.glob(indir_toe_rcphn + '/*.nc')
nmodels = len(listfiles)

for i in range(nmodels):

    file_toe = listfiles[i]
    ftoe = open_ncfile(file_toe, 'r')
    name = os.path.basename(file_toe).split('.')[1]
    # Read ToE (members, basin density, latitude)
    toe2read = ftoe.variables[var + 'ToE2'][:]
    nMembers[i] = toe2read.shape[0]
    print('- Reading ToE of %s with %d members'%(name,nMembers[i]))
    nruns1 = nruns + nMembers[i]

    # Save ToE
    varToEA[nruns:nruns1,:,:] = toe2read[:,1,:,:]
    varToEP[nruns:nruns1,:,:] = toe2read[:,2,:,:]
    varToEI[nruns:nruns1,:,:] = toe2read[:,3,:,:]

    nruns = nruns1

print('Total number of runs:', nruns)
varToEA = varToEA[0:nruns,:,:]
varToEP = varToEP[0:nruns,:,:]
varToEI = varToEI[0:nruns,:,:]

nruns = int(nruns)

# ----- Build plot variables ------

# # -- Mask
# var_mask = np.ma.getmask(np.ma.average(fh2d.variables[var][:], axis=0))
# varToEA = np.ma.array(varToEA, mask=var_mask[1,:,:])
# varToEP = np.ma.array(varToEP, mask=var_mask[2,:,:])
# varToEI = np.ma.array(varToEI, mask=var_mask[3,:,:])


# -- Compute median and range
# Median
medianToEA = np.ma.around(np.ma.median(varToEA, axis=0)) + iniyear
medianToEP = np.ma.around(np.ma.median(varToEP, axis=0)) + iniyear
medianToEI = np.ma.around(np.ma.median(varToEI, axis=0)) + iniyear
# 16th percentile
percentile16ToEA = np.percentile(varToEA, 16, axis=0)
percentile16ToEP = np.percentile(varToEP, 16, axis=0)
percentile16ToEI = np.percentile(varToEI, 16, axis=0)
# 84th percentile
percentile84ToEA = np.percentile(varToEA, 84, axis=0)
percentile84ToEP = np.percentile(varToEP, 84, axis=0)
percentile84ToEI = np.percentile(varToEI, 84, axis=0)
# 16-84% range
rangeToEA = np.ma.around(percentile84ToEA - percentile16ToEA)
rangeToEP = np.ma.around(percentile84ToEP - percentile16ToEP)
rangeToEI = np.ma.around(percentile84ToEI - percentile16ToEI)

# -- Mask points where signal has not emerged
medianToEA[medianToEA == finalyear] = np.ma.masked
medianToEP[medianToEP == finalyear] = np.ma.masked
medianToEI[medianToEI == finalyear] = np.ma.masked


# -- Create variable bundles
labBowl = ['0','0']
varAtl = {'name': 'Atlantic', 'ToE': medianToEA, 'bowl2': None, 'bowl1': None, 'labBowl': labBowl}
varPac = {'name': 'Pacific', 'ToE': medianToEP, 'bowl2': None, 'bowl1': None, 'labBowl': labBowl}
varInd = {'name': 'Indian', 'ToE': medianToEI, 'bowl2': None, 'bowl1': None, 'labBowl': labBowl}


# ----- Plot ToE ------

fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(17,5))

minmax = [iniyear, finalyear +1, deltat]
unit = 'ToE'
cmap = cmaps.viridis
levels = np.arange(minmax[0], minmax[1], minmax[2])

cnplot = zonal_2D(plt, 'ToE', axes[0, 0], axes[1, 0], 'left', lat, density, varAtl, domrho, cmap, levels)

cnplot = zonal_2D(plt, 'ToE', axes[0, 1], axes[1, 1], 'mid', lat, density, varPac, domrho, cmap, levels)

cnplot = zonal_2D(plt, 'ToE', axes[0, 2], axes[1, 2], 'right', lat, density, varInd, domrho, cmap, levels)

plt.subplots_adjust(hspace=.0001, wspace=0.05, left=0.04, right=0.86)

cb = plt.colorbar(cnplot, ax=axes.ravel().tolist(), ticks = levels[::2], fraction=0.015, shrink=2.0, pad=0.05)
cb.set_label('%s' % (unit,), fontweight='bold')

plotTitle = 'Multimodel ensemble median ToE for ' + legVar + ', hist+RCP8.5 vs. histNat [> ' + str(multStd) + ' std]'
# ADD NUMBER OF RUNS

plt.suptitle(plotTitle, fontweight='bold', fontsize=14, verticalalignment='top')

plt.figtext(.5,.02,'Computed by : zonal_toe_median_range.py',fontsize=9,ha='center')

plt.show()
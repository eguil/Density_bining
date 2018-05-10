#!/bin/env python
# -*- coding: utf-8 -*-

import os
import glob
from netCDF4 import Dataset as open_ncfile
import matplotlib.pyplot as plt
from maps_matplot_lib import defVarmme, remapToZ, zon_2Dz
from modelsDef import defModels
import numpy as np
import colormaps as cmaps

# ----- Workspace ------

indir_toe_rcphn = '/home/ysilvy/Density_bining/Yona_analysis/data/toe_zonal/toe_rcp85_histNat/'
indir_mme_rcp85 = '/data/ericglod/Density_binning/Prod_density_april15/mme_rcp85/'
indir_mme_hn = '/data/ericglod/Density_binning/Prod_density_april15/mme_histNat/'

models = defModels()

# ----- Work ------

varname = defVarmme('salinity'); v = 'S'

multStd = 2. # detect ToE at multStd std dev of piControl

use_piC = False # Over projection period, signal = RCP-average(histNat), noise = std(histNat)
# use_piC = True # Over projection period, signal = RCP-average(PiControl), noise = std(PiControl)

iniyear = 1860
finalyear = 2100
deltat = 10.

# ----- Variables ------
var = varname['var_zonal_w/bowl']
legVar = varname['legVar']

# Read latitude and density from original file
fileh_2d = '/data/ericglod/Density_binning/Prod_density_april15/mme_hist/' \
       'cmip5.multimodel_Nat.historical.ensm.an.ocn.Omon.density_zon2D.nc'
fh2d = open_ncfile(fileh_2d, 'r')
lat = fh2d.variables['latitude'][:]; latN = len(lat)
density = fh2d.variables['lev'][:]; levN = len(density)

valmask = 1.e20
basinN = 4

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

# -- Mask points
rangeToEA[rangeToEA == 0] = np.ma.masked
rangeToEP[rangeToEP == 0] = np.ma.masked
rangeToEI[rangeToEI == 0] = np.ma.masked

# -- Regroup basins
medianToEr = np.ma.masked_all((basinN,levN,latN))
medianToEr[1,:,:] = medianToEA
medianToEr[2,:,:] = medianToEP
medianToEr[3,:,:] = medianToEI

rangeToEr = np.ma.masked_all((basinN,levN,latN))
rangeToEr[1,:,:] = rangeToEA
rangeToEr[2,:,:] = rangeToEP
rangeToEr[3,:,:] = rangeToEI

# ----- Read volume and depth of isopycnals ------

# Read files
file1d_rcp85 = 'cmip5.multimodel_Nat.rcp85.ensm.an.ocn.Omon.density_zon1D.nc'
f1drcp85 = open_ncfile(indir_mme_rcp85+file1d_rcp85,'r')
file2d_rcp85 = 'cmip5.multimodel_Nat.rcp85.ensm.an.ocn.Omon.density_zon2D.nc'
f2drcp85 = open_ncfile(indir_mme_rcp85+file2d_rcp85,'r')
file_hn = 'cmip5.multimodel_Nat.historicalNat.ensm.an.ocn.Omon.density_zon1D.nc'
fhn = open_ncfile(indir_mme_hn+file_hn,'r')

# Read variables
depthr = f2drcp85.variables['isondepth'][-5:,:,:,:]
volumr = f2drcp85.variables['isonvol'][-5:,:,:,:]
bowl2z = f1drcp85.variables['ptopdepth'][-5:,:,:]
bowl1z = fhn.variables['ptopdepth'][-5:,:,:]

# Average
depthr = np.ma.average(depthr, axis=0)
volumr = np.ma.average(volumr, axis=0)
bowl2z = np.ma.average(bowl2z, axis=0)
bowl1z = np.ma.average(bowl1z, axis=0)

labBowl = ['histNat', 'RCP8.5']


# ----- Bathymetry ------

# Read masks
fmask = open_ncfile('/home/ysilvy/Density_bining/Yona_analysis/data/170224_WOD13_masks.nc','r')
basinmask = fmask.variables['basinmask3'][:] # (latitude, longitude)
depthmask = fmask.variables['depthmask'][:] # (latitude, longitude)
# Create basin masks
mask_a = basinmask != 1
mask_p = basinmask != 2
mask_i = basinmask != 3
# Take max value along lon for each latitude of each basin --> gives us the bathymetry
depthmask_a = np.ma.array(depthmask, mask=mask_a)
bathy_a = np.ma.max(depthmask_a, axis=1)
depthmask_p = np.ma.array(depthmask, mask=mask_p)
bathy_p = np.ma.max(depthmask_p, axis=1)
depthmask_i = np.ma.array(depthmask, mask=mask_i)
bathy_i = np.ma.max(depthmask_i, axis=1)

bathy = np.ma.masked_all((basinN,len(lat)))
bathy[1,:] = bathy_a
bathy[2,:] = bathy_p
bathy[3,:] = bathy_i


# ----- Remapping -----

# WOA13 grid
targetz = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80,
   85, 90, 95, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 375,
   400, 425, 450, 475, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950,
   1000, 1050, 1100, 1150, 1200, 1250, 1300, 1350, 1400, 1450, 1500, 1550,
   1600, 1650, 1700, 1750, 1800, 1850, 1900, 1950, 2000, 2100, 2200, 2300,
   2400, 2500, 2600, 2700, 2800, 2900, 3000, 3100, 3200, 3300, 3400, 3500,
   3600, 3700, 3800, 3900, 4000, 4100, 4200, 4300, 4400, 4500, 4600, 4700,
   4800, 4900, 5000, 5100, 5200, 5300, 5400, 5500]

depthrbis = depthr

medianToEz = np.ma.around(remapToZ(medianToEr.data, depthrbis, volumr, targetz, bowl1z, 'ToE', bathy))
print('Median remapping done')
rangeToEz = np.ma.around(remapToZ(rangeToEr.data, depthr, volumr, targetz, bowl1z, 'ToE', bathy))
print('Range remapping done')

# Mask
# medianToEz[medianToEz == finalyear] = np.ma.masked
# medianToEz[medianToEz < iniyear ] = np.ma.masked
rangeToEz[medianToEz > finalyear-20] = np.ma.masked # Mask points where median hasn't emerged

# -- Make variable bundles for each basin
varAtlmedian = {'name': 'Atlantic', 'var_change': medianToEz[1,:,:], 'bowl1': bowl1z[1,:], 'bowl2': bowl2z[1,:], 'labBowl': labBowl}
varPacmedian = {'name': 'Pacific', 'var_change': medianToEz[2,:,:], 'bowl1': bowl1z[2,:], 'bowl2': bowl2z[2,:], 'labBowl': labBowl}
varIndmedian = {'name': 'Indian', 'var_change': medianToEz[3,:,:], 'bowl1': bowl1z[3,:], 'bowl2': bowl2z[3,:], 'labBowl': labBowl}

varAtlrange = {'name': 'Atlantic', 'var_change': rangeToEz[1,:,:], 'bowl1': bowl1z[1,:], 'bowl2': bowl2z[1,:], 'labBowl': labBowl}
varPacrange = {'name': 'Pacific', 'var_change': rangeToEz[2,:,:], 'bowl1': bowl1z[2,:], 'bowl2': bowl2z[2,:], 'labBowl': labBowl}
varIndrange = {'name': 'Indian', 'var_change': rangeToEz[3,:,:], 'bowl1': bowl1z[3,:], 'bowl2': bowl2z[3,:], 'labBowl': labBowl}


# ----- Plot -----

domzed = [0,500,4000]

# -- Median

# -- Create figure and axes instances
fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(17, 5))

# -- Color map
cmap = 'jet_r'
# -- Unit
unit = 'ToE'
# -- Levels
minmax = [1950, 2080.01, deltat]
levels = np.arange(minmax[0], minmax[1], minmax[2])

# -- Contourf
cnplot = zon_2Dz(plt, axes[0,0], axes[1,0], 'left', lat, targetz, varAtlmedian,
                 None, cmap, levels, domzed)
cnplot[0].cmap.set_over('white')
cnplot[1].cmap.set_over('white')

cnplot = zon_2Dz(plt, axes[0,1], axes[1,1], 'mid', lat, targetz, varPacmedian,
                 None, cmap, levels, domzed)
cnplot[0].cmap.set_over('white')
cnplot[1].cmap.set_over('white')

cnplot = zon_2Dz(plt, axes[0,2], axes[1,2], 'right', lat, targetz, varIndmedian,
                 None, cmap, levels, domzed)
cnplot[0].cmap.set_over('white')
cnplot[1].cmap.set_over('white')

plt.subplots_adjust(hspace=.0001, wspace=0.05, left=0.04, right=0.86)

# -- Add colorbar
cb = plt.colorbar(cnplot[1], ax=axes.ravel().tolist(), ticks=levels, fraction=0.015, shrink=2.0, pad=0.05)
cb.set_label('%s' % (unit,), fontweight='bold')

# -- Add title
plotTitle = 'Multimodel ensemble median ToE for ' + legVar + ', hist+RCP8.5 vs. histNat [> ' + str(multStd) + ' std]' \
    '\n %d models, %d runs '%(nmodels,nruns)

plt.suptitle(plotTitle, fontweight='bold', fontsize=14, verticalalignment='top')

plt.figtext(.5,.02,'Computed by : remap_to_z_toe_median_range.py',fontsize=9,ha='center')

# figureDir = 'models/zonal_remaptoz/'
# plotName = 'remapping_median_toe_rcp85vshistNat'
# plt.savefig('/home/ysilvy/Density_bining/Yona_analysis/figures/'+figureDir+plotName+'.png', bbox_inches='tight')

# -- 16-84% range

# -- Create figure and axes instances
fig2, axes = plt.subplots(nrows=2, ncols=3, figsize=(17, 5))

# -- Color map
cmap = 'jet_r'
# -- Unit
unit = 'Years'
# -- Levels
minmax = [0, 121, deltat]
levels = np.arange(minmax[0], minmax[1], minmax[2])

# -- Contourf
cnplot2 = zon_2Dz(fig2, axes[0,0], axes[1,0], 'left', lat, targetz, varAtlrange,
                 None, cmap, levels, domzed)
cnplot2 = zon_2Dz(fig2, axes[0,1], axes[1,1], 'mid', lat, targetz, varPacrange,
                 None, cmap, levels, domzed)
cnplot2 = zon_2Dz(fig2, axes[0,2], axes[1,2], 'right', lat, targetz, varIndrange,
                 None, cmap, levels, domzed)

fig2.subplots_adjust(hspace=.0001, wspace=0.05, left=0.04, right=0.86)

# -- Add colorbar
cb = fig2.colorbar(cnplot2[1], ax=axes.ravel().tolist(), ticks=levels, fraction=0.015, shrink=2.0, pad=0.05)
cb.set_label('%s' % (unit,), fontweight='bold')

# -- Add title
plotTitle = '16-84% multimodel ensemble range of the ToE for ' + legVar + ', hist+RCP8.5 vs. histNat [> ' + str(multStd) + ' std]' \
    '\n %d models, %d runs '%(nmodels,nruns)

plt.suptitle(plotTitle, fontweight='bold', fontsize=14, verticalalignment='top')

plt.figtext(.5,.02,'Computed by : remap_to_z_toe_median_range.py',fontsize=9,ha='center')

# # figureDir = 'models/zonal_remaptoz/'
# # plotName = 'remapping_16-84range_toe_rcp85vshistNat'
# # plt.savefig('/home/ysilvy/Density_bining/Yona_analysis/figures/'+figureDir+plotName+'.png', bbox_inches='tight')
#
#
# plt.show()

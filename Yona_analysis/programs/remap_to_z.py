#!/bin/env python
# -*- coding: utf-8 -*-

import os
from netCDF4 import Dataset as open_ncfile
import matplotlib.pyplot as plt
from maps_matplot_lib import defVarmme, remapToZ, zon_2Dz, custom_div_cmap
import numpy as np


# -------------------------------------------------------------------------------
#                                Define work
# -------------------------------------------------------------------------------

inDir = '/data/ericglod/Density_binning/'
workh = 'Prod_density_april15/mme_hist'
workhn = 'Prod_density_april15/mme_histNat'
inDirh = inDir + workh
inDirhn = inDir + workhn

# output format
outfmt = 'view'
#outfmt = 'save'

# Define variable
varname = defVarmme('salinity'); v = 'S'
varname = defVarmme('temp'); v = 'T'
#varname = defVarmme('volume'); v = 'V'

minmax = varname['minmax_zonal']
clevsm = varname['clevsm_zonal']
legVar = varname['legVar']
unit = varname['unit']

iniyear = 1860
finalyear = 2005

plotName = 'cmip5_remap_test' + varname['var_zonal']


os.chdir(inDir)
fileh2d = 'cmip5.multimodel_Nat.historical.ensm.an.ocn.Omon.density_zon2D.nc'
fileh1d = 'cmip5.multimodel_Nat.historical.ensm.an.ocn.Omon.density_zon1D.nc'
filehn2d = 'cmip5.multimodel_All.historicalNat.ensm.an.ocn.Omon.density_zon2D.nc'
filehn1d = 'cmip5.multimodel_All.historicalNat.ensm.an.ocn.Omon.density_zon1D.nc'

var = varname['var_zonal_w/bowl']

# find dimensions
fh2d = open_ncfile(inDirh+'/'+fileh2d)
fh1d = open_ncfile(inDirh+'/'+fileh1d)
fhn2d = open_ncfile(inDirhn+'/'+filehn2d)
fhn1d = open_ncfile(inDirhn+'/'+filehn1d)

# -- Read variables
fieldhr = fh2d.variables[var][-5:,:,:,:]
fieldhnr = fhn2d.variables[var][-5:,:,:]
depthr = fh2d.variables['isondepth'][-5:,:,:,:]
volumr = fh2d.variables['isonvol'][-5:,:,:,:]
lat = fh2d.variables['latitude'][:]
bowlhz = fh1d.variables['ptopdepth'][-5:,:,:]
bowlhnz = fhn1d.variables['ptopdepth'][-5:,:,:]

valmask = 1.e20
basinN = 4

# -------------------------------------------------------------------------------
#                                Build variables
# -------------------------------------------------------------------------------

# -- Bathymetry
# Read masks
fmask = open_ncfile('/home/ysilvy/Density_bining/Yona_analysis/data/170224_WOD13_masks.nc','r')
basinmask = fmask.variables['basinmask3'][:] #(latitude, longitude)
depthmask = fmask.variables['depthmask'][:] #(latitude, longitude)
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

# -- Compute signal hist - histNat
vardiffr = np.ma.average(fieldhr, axis=0) - np.ma.average(fieldhnr, axis=0)
# -- Average other variables
depthr = np.ma.average(depthr, axis=0)
volumr = np.ma.average(volumr, axis=0)
bowlhz = np.ma.average(bowlhz, axis=0)
bowlhnz = np.ma.average(bowlhnz, axis=0)


# -- Target grid for remapping

#targetz = [0.,5.,15.,25.,35.,45.,55.,65.,75.,85.,95.,105.,116.,128.,142.,158.,181.,216.,272.,364.,511.,732.,1033.,1405.,1830.,2289.,2768.,3257.,3752.,4350.,4749.,5250.]

# WOA09 grid
# targetz = [0, 10, 20, 30, 50, 75, 100, 125, 150, 200, 250, 300, 400, 500, 600,
# 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1750, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 5500]

# WOA13 grid
targetz = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80,
   85, 90, 95, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 375,
   400, 425, 450, 475, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950,
   1000, 1050, 1100, 1150, 1200, 1250, 1300, 1350, 1400, 1450, 1500, 1550,
   1600, 1650, 1700, 1750, 1800, 1850, 1900, 1950, 2000, 2100, 2200, 2300,
   2400, 2500, 2600, 2700, 2800, 2900, 3000, 3100, 3200, 3300, 3400, 3500,
   3600, 3700, 3800, 3900, 4000, 4100, 4200, 4300, 4400, 4500, 4600, 4700,
   4800, 4900, 5000, 5100, 5200, 5300, 5400, 5500]


# -- Remap
fieldz = remapToZ(vardiffr.data, depthr.data, volumr.data, targetz, bowlhz, v, bathy)


# -- Make variable bundles for each basin
varAtl = {'name': 'Atlantic', 'var_change': fieldz[1,:,:], 'bowl1': bowlhnz[1,:], 'bowl2': bowlhz[1,:]}
varPac = {'name': 'Pacific', 'var_change': fieldz[2,:,:], 'bowl1': bowlhnz[2,:], 'bowl2': bowlhz[2,:]}
varInd = {'name': 'Indian', 'var_change': fieldz[3,:,:], 'bowl1': bowlhnz[3,:], 'bowl2': bowlhz[3,:]}


# -------------------------------------------------------------------------------
#                                Plot
# -------------------------------------------------------------------------------

domzed = [0,500,2000]

# -- Create figure and axes instances
fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(17, 5))

# -- color map
cmap = custom_div_cmap()

# -- levels
levels = np.linspace(minmax[0],minmax[1],minmax[2]) # Change
#levels = np.arange(33.5,35.5,0.2) # Mean salinity
#levels = np.arange(-2,30,2) # Mean temperature
#levels = np.arange(0,701,50) # Mean volume

# -- Contourf of signal
cnplot = zon_2Dz(plt, axes[0,0], axes[1,0], 'left', lat, targetz, varAtl,
                 clevsm, cmap, levels, domzed)
cnplot = zon_2Dz(plt, axes[0,1], axes[1,1], 'mid', lat, targetz, varPac,
                 clevsm, cmap, levels, domzed)
cnplot = zon_2Dz(plt, axes[0,2], axes[1,2], 'right', lat, targetz, varInd,
                 clevsm, cmap, levels, domzed)


plt.subplots_adjust(hspace=.0001, wspace=0.05, left=0.04, right=0.86)

# -- Add colorbar
cb = plt.colorbar(cnplot, ax=axes.ravel().tolist(), ticks=levels[::3], fraction=0.015, shrink=2.0, pad=0.05)
cb.set_label('%s (%s)' % (legVar, unit), fontweight='bold')

# add Title text
fig.suptitle('Remapping %s changes hist-histNat (mme) from rho to z' %(legVar,), fontsize=14, fontweight='bold')

plt.figtext(.5,.01,'Computed by : remap_to_z.py',fontsize=9,ha='center')

# -- Output  # TODO read as argument

if outfmt == 'view':
    plt.show()



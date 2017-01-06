#!/bin/env python
# -*- coding: utf-8 -*-

"""
Python matplotlib
Make density projections of temperature and salinity change between hist and histNat, averaged between 2000 and 2005

"""

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset as open_ncfile
from maps_matplot_lib import defVar, proj_map


# ------ Read files -----------

# Data path
indir = '/data/ericglod/Density_binning/Prod_density_april15/Raw/'
file_hist = 'cmip5.CanESM2.historical.r2i1p1.an.ocn.Omon.density.ver-1.nc'
file_histNat = 'cmip5.CanESM2.historicalNat.r2i1p1.an.ocn.Omon.density.ver-1.nc'
name = 'CanESM2_r2i1p1'
data_hist = indir + 'historical/correct/' + file_hist
data_histNat = indir + 'historicalNat/' + file_histNat

fhist = open_ncfile(data_hist,'r') #1850-->2005
fhistNat = open_ncfile(data_histNat,'r') #1850-->2012

# Read variables
density = fhist.variables['lev'][:]
lat = fhist.variables['latitude'][:]
lon = fhist.variables['longitude'][:]


# ------ Define variables ----------

# Choose which variable to work on
#varname = defVar('temp')
varname = defVar('salinity')

# Define variable properties
var = varname['var']
minmax = [-2.5,2.5,20] #varname['minmax']
clevsm = varname['clevsm']
legVar = varname['legVar']
unit = varname['unit']


# Define isopycnal on which we want to project the temperature/salinity/variable
isopyc1 = 24.0
isopyc2 = 25.0
isopyc3 = 26.7
isopyc4 = 27.5
isopyc1_idx = np.argmin(np.abs(density - isopyc1))
isopyc4_idx = np.argmin(np.abs(density - isopyc4))


var_hist = fhist.variables[var][:, slice(isopyc1_idx,isopyc4_idx+1), :]
var_histNat = fhistNat.variables[var][:, slice(isopyc1_idx,isopyc4_idx+1), :]
sliced_density = density[isopyc1_idx:isopyc4_idx+1]


# ----- Plot changes and means ---------

fig = plt.figure()
ax1 = fig.add_axes([0.05, 0.74, 0.85, 0.2])
ax2 = fig.add_axes([0.05, 0.51, 0.85, 0.2])
ax3 = fig.add_axes([0.05, 0.28, 0.85, 0.2])
ax4 = fig.add_axes([0.05, 0.05, 0.85, 0.2])

cmap = plt.get_cmap('bwr')  # red/white/blue difference map

bmap = proj_map('model', plt, ax1, minmax, clevsm, lat, lon, cmap, isopyc1, sliced_density, var_hist, var2 = var_histNat)
bmap = proj_map('model', plt, ax2, minmax, clevsm, lat, lon, cmap, isopyc2, sliced_density, var_hist, var2 = var_histNat)
bmap = proj_map('model', plt, ax3, minmax, clevsm, lat, lon, cmap, isopyc3, sliced_density, var_hist, var2 = var_histNat)
bmap = proj_map('model', plt, ax4, minmax, clevsm, lat, lon, cmap, isopyc4, sliced_density, var_hist, var2 = var_histNat)


cb = plt.colorbar(bmap[0], ax = (ax1, ax2, ax3, ax4), ticks=bmap[1], orientation='vertical')
cb.set_label('%s (%s)' % (legVar, unit), fontweight='bold')
plt.suptitle('%s changes between hist and histNat averaged between 2000 and 2005 (%s)' %(legVar,name), fontweight='bold', fontsize=16)
plt.show()

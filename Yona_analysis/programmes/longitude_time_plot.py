#!/bin/env python
# -*- coding: utf-8 -*-

"""
Python matplotlib
Make longitude/density diagrams of temperature and salinity change between 2010 and 1945, from observations

"""

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset as open_ncfile
from maps_matplot_lib import defVar, hovmoller
from matplotlib.ticker import MaxNLocator


# ------ Read files -----------

# Data path
indir = '/data/ericglod/Density_binning/Prod_density_obs_april16/'
#file = 'obs.EN4.historical.r0i0p0.mo.ocn.Omon.density.ver-1.latestXCorr.nc'
#name = 'EN4'
file = 'obs.Ishii.historical.r0i0p0.an.ocn.Omon.density.ver-1.latestXCorr.nc'
name = 'Ishii'

data = indir + file

f = open_ncfile(data,'r')

# Read variables
density = f.variables['lev'][:]
lat = f.variables['latitude'][:]
lon = f.variables['longitude'][:]
time = f.variables['time'][:]

# ------ Define variables ----------

# Choose which variable to work on
#varname = defVar('temp')
varname = defVar('salinity')

# Choose latitude
lat1 = -60
lat1_idx = np.argmin(np.abs(lat - lat1))

# Create new time array
if name=='Ishii':
    t_start = 1945
if name=='EN4':
    t_start = 1900
date = np.zeros(np.shape(time))
date[0] = t_start
for i in range(1,len(date)):
    date[i] = date[i-1] +1


# Define variable properties
var = varname['var']
clevsm = varname['clevsm']
legVar = varname['legVar']
unit = varname['unit']

var = f.variables[var][:,:,lat1_idx,:].squeeze()

# Define isopycnal on which we want to project the temperature/salinity/variable
isopyc1 = 24.0
isopyc2 = 25.0
isopyc3 = 26.7
isopyc4 = 27.5
isopyc1_idx = np.argmin(np.abs(density - isopyc1))
isopyc4_idx = np.argmin(np.abs(density - isopyc4))

var = var[:, isopyc1_idx:isopyc4_idx+1, :]
sliced_density = density[isopyc1_idx:isopyc4_idx+1]


# ----------- Build plots ------------


fig, ax = plt.subplots(4,1)
# ax1 = plt.gca()
# ax1 = fig.add_axes([0.05, 0.74, 0.85, 0.2])
# ax2 = fig.add_axes([0.05, 0.51, 0.85, 0.2])
# ax3 = fig.add_axes([0.05, 0.28, 0.85, 0.2])
# ax4 = fig.add_axes([0.05, 0.05, 0.85, 0.2])

cmap = cmap=plt.get_cmap('RdYlBu_r')

plot = hovmoller(plt, ax[0], lon, date, clevsm, legVar, unit, cmap, isopyc1, sliced_density, var)
plot = hovmoller(plt, ax[1], lon, date, clevsm, legVar, unit, cmap, isopyc2, sliced_density, var)
plot = hovmoller(plt, ax[2], lon, date, clevsm, legVar, unit, cmap, isopyc3, sliced_density, var)
plot = hovmoller(plt, ax[3], lon, date, clevsm, legVar, unit, cmap, isopyc4, sliced_density, var)

ax[3].set_xlabel('Longitude')

plt.suptitle('Hovmoller plot of %s at different densities, at lat=%d (%s)' %(legVar, lat1, name), fontweight='bold', fontsize=16)

plt.show()


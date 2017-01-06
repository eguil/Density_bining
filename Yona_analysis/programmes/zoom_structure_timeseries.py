#!/bin/env python
# -*- coding: utf-8 -*-

"""
Python matplotlib
Zoom on a specific structure (from density projection maps), and plot evolution with time of average temperature and salinity

"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from netCDF4 import Dataset as open_ncfile
from maps_matplot_lib import defVar, proj_map, proj_map_zonal_changes


# ------ Read files -----------

# Data path
indir = '/data/ericglod/Density_binning/Prod_density_obs_april16/'

file = 'obs.Ishii.historical.r0i0p0.an.ocn.Omon.density.ver-1.latestXCorr.nc'
name = 'Ishii'

#file = 'obs.EN4.historical.r0i0p0.mo.ocn.Omon.density.ver-1.latestXCorr.nc'
#name = 'EN4'

data = indir + file

f = open_ncfile(data,'r')

# Read variables
density = f.variables['lev'][:]
lat = f.variables['latitude'][:]
lon = f.variables['longitude'][:]
time = f.variables['time'][:]


# ------ Define and build variables ----------
# Choose 2D or 3D average
#av = '2D'
av = '3D'

temp = defVar('temp')
salinity = defVar('salinity')

# Define variable properties
var_T = temp['var']
legVar_T = temp['legVar']
unit_T = temp['unit']
var_S = salinity['var']
legVar_S = salinity['legVar']
unit_S = salinity['unit']

# Re-create time array
if name=='Ishii':
    t_start = 1945
if name=='EN4':
    t_start = 1900

date = np.zeros(np.shape(time))
date[0] = t_start
for i in range(1,len(date)):
    date[i] = date[i-1] +1

if av == '2D':
    # Define isopycnal on which we want to project the temperature/salinity
    isopyc = 26.5 ; isopyc_idx = np.argmin(np.abs(density - isopyc))
else :
    isopycmin = 25 ; isopycmin_idx = np.argmin(np.abs(density - isopycmin))
    isopycmax = 27 ; isopycmax_idx = np.argmin(np.abs(density - isopycmax))

# Shift lon to make it start at 20E so that pacific ocean is not cut in half
lon = np.roll(lon,-200)

# Define lat/lon box
latmin = -50 ; latmin_idx = np.argmin(np.abs(lat-latmin))
latmax = -5 ; latmax_idx = np.argmin(np.abs(lat-latmax))
lonmin = 150 ; lonmin_idx = np.argmin(np.abs(lon-lonmin))
lonmax = -90 ; lonmax_idx = np.argmin(np.abs(lon-lonmax))

# Read variables at chosen lats/lons/density
if av == '2D':
    temp = f.variables[var_T][:, isopyc_idx, latmin_idx:latmax_idx+1,:].squeeze()
    salinity = f.variables[var_S][:, isopyc_idx, latmin_idx:latmax_idx + 1, :].squeeze()
    temp = np.roll(temp, -200, axis=2)
    temp = temp[:, :, lonmin_idx:lonmax_idx + 1]
    salinity = np.roll(salinity, -200, axis=2)
    salinity = salinity[:, :, lonmin_idx:lonmax_idx + 1]
    # Average spatially to obtain 1D time series
    temp_av = np.ma.average(temp, axis=1)
    temp_av = np.ma.average(temp_av, axis=1)
    salinity_av = np.ma.average(salinity, axis=1)
    salinity_av = np.ma.average(salinity_av, axis=1)

else :
    temp = f.variables[var_T][:, isopycmin_idx:isopycmax_idx, latmin_idx:latmax_idx+1,:]
    salinity = f.variables[var_S][:, isopycmin_idx:isopycmax_idx, latmin_idx:latmax_idx + 1, :]
    temp = np.roll(temp, -200, axis=3)
    temp = temp[:, :, :, lonmin_idx:lonmax_idx + 1]
    salinity = np.roll(salinity, -200, axis=3)
    salinity = salinity[:, :, :, lonmin_idx:lonmax_idx + 1]
    # Average spatially to obtain 1D time series
    temp_av = np.ma.average(temp, axis=2)
    temp_av = np.ma.average(temp_av, axis=2)
    temp_av = np.ma.average(temp_av, axis=1)
    salinity_av = np.ma.average(salinity, axis=2)
    salinity_av = np.ma.average(salinity_av, axis=2)
    salinity_av = np.ma.average(salinity_av, axis=1)


# ------ Plot -------

fig, ax = plt.subplots(2)
ax[0].plot(date, temp_av)
ax[0].set_ylabel('Temperature', fontweight='bold')
plt.suptitle('Temperature and Salinity in the South Pacific Ocean between' + r' $\gamma = [%.1f:%.1f]$ ' %(isopycmin,isopycmax) + ' (%s)' %(name,),
                va='top', fontweight='bold')
ax[0].set_xlim([t_start,date[-1]])
ax[1].plot(date, salinity_av)
ax[1].set_ylabel('Salinity', fontweight='bold')
ax[1].set_xlabel('Time', fontweight='bold')
ax[1].set_xlim([t_start,date[-1]])
y_formatter = ScalarFormatter(useOffset=False)
ax[1].yaxis.set_major_formatter(y_formatter)

plt.show()

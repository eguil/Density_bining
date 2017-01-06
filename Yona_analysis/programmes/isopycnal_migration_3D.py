#!/bin/env python
# -*- coding: utf-8 -*-

"""
Python matplotlib
Compute salinity changes due to isopycnal migration and plot maps projected on isopycnals

"""

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset as open_ncfile
from maps_matplot_lib import custom_div_cmap
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1 import make_axes_locatable

# ----- Workspace ------

indir = '/data/ericglod/Density_binning/Prod_density_obs_april16/'
file = 'DurackandWijffels_GlobalOceanChanges-NeutralDensity_1950-2000_120209_11_46_11_beta.nc'
name = 'Durack & Wijffels'

data = indir + file

f = open_ncfile(data,'r')

# ----- Variables ------

# Read variables
density = f.variables['density'][:]
lat = f.variables['latitude'][:]
lon = f.variables['longitude'][:]

temp = {'var_change':'thetao_change', 'var_mean':'thetao_mean', 'minmax': [-0.65, 0.65, 14],
        'clevsm': np.arange(-2, 30, 1), 'legVar': "Temperature", 'unit': "C", 'longN': 'temp'}
salinity = {'var_change':'salinity_change', 'var_mean':'salinity_mean', 'minmax': [-0.3, 0.3, 16],
            'clevsm': np.arange(30, 40, .2), 'legVar': "Salinity", 'unit': "PSU", 'longN': 'salinity'}

#varname = temp
varname = salinity

# Define variable properties
minmax = varname['minmax']
clevsm = varname['clevsm']
legVar = varname['legVar']
unit = varname['unit']
var_mean = varname['var_mean']

var_attributes = f.variables[var_mean]
var_mean = f.variables[var_mean][:].squeeze()


# Define isopycnal on which we want to project the temperature/salinity/variable
isopyc1 = 24.0
isopyc2 = 25.0
isopyc3 = 26.7
isopyc4 = 27.5
isopyc1_idx = np.argmin(np.abs(density - isopyc1))
isopyc3_idx = np.argmin(np.abs(density - isopyc3))
isopyc4_idx = np.argmin(np.abs(density - isopyc4))
sliced_density = density[isopyc1_idx:isopyc4_idx+1]

# ----- Build partial derivative of mean field for density
latN = len(lat)
lonN = len(lon)
valmask = var_attributes.missing_value

var_mean_rshp = np.reshape(var_mean, (len(density),latN*lonN))
#print('var_mean_rshp shape: ', np.shape(var_mean_rshp))
var_deriv = np.ma.ones(np.shape(var_mean_rshp))*valmask
#print('var_deriv shape: ', np.shape(var_deriv))

var_mean_rshp1 = np.roll(var_mean_rshp, -1, axis=0)
density1 = np.roll(density, -1)
delta_sig = density1-density
delta_sig = np.tile(delta_sig,(latN*lonN,1)); delta_sig = np.transpose(delta_sig)
#print(np.shape(delta_sig))

var_deriv = (var_mean_rshp1 - var_mean_rshp)/delta_sig
var_deriv[-1,:] = valmask
var_deriv = np.reshape(var_deriv, (len(density), latN, lonN))
#print(np.shape(var_deriv))


# ------ Plot --------

cmap = plt.get_cmap('bwr')
#cmap = custom_div_cmap()

fig = plt.figure()

var_deriv = np.squeeze(var_deriv[isopyc3_idx,:])

# Create meshgrid
lon2d, lat2d = np.meshgrid(lon, lat)

levels = np.linspace(-2, 2, 16)

# Basemap
map = Basemap(projection='cyl', llcrnrlon=20, llcrnrlat=-70., urcrnrlon=380, urcrnrlat=70)
map.drawmapboundary(fill_color='1')
map.drawparallels(np.arange(-60, 61, 20.), labels=[1, 0, 1, 0], linewidth=0.5)
map.drawmeridians(np.arange(-180, 180, 60), labels=[0, 0, 0, 1], linewidth=0.5)
map.fillcontinents(color='black')

# Draw filled contours of diff
cnplot = map.contourf(lon2d, lat2d, var_deriv, levels=levels, cmap=cmap, latlon=True, extend='both')

# Make colorbar same height as map
ax = plt.gca()
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="3%", pad=0.15)
cb = plt.colorbar(cnplot, ticks=levels, cax=cax)
#cb.set_label('%s (%s)' % (legVar, unit), fontweight='bold')

plt.show()

#!/bin/env python
# -*- coding: utf-8 -*-

"""
Python matplotlib
Make latitude/density diagrams of temperature and salinity changes

"""


import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset as open_ncfile
from maps_matplot_lib import defVarDurack, custom_div_cmap

# ------ Define work and Read files -----------

name = 'Durack&Wijffels'

# Data path
indir = '/data/ericglod/Density_binning/Obs_Prod_density_april16/'
file = 'DurackandWijffels_GlobalOceanChanges-NeutralDensity_1950-2000_120209_11_46_11_beta.nc'

data = indir + file

f = open_ncfile(data,'r')

# Read variables
lat = f.variables['latitude'][:]
lon = f.variables['longitude'][:]

# ------ Define variables ----------

# Choose which variable to work on
varname = defVarDurack('salinity')
#varname = defVarDurack('temp')

# Choose latitude (one specific lat or average a range of lats)
#lon_choice = 'one'
lon_choice = 'range'

if lon_choice == 'one':
    lon1 = 60
    lon1_idx = np.argmin(np.abs(lon - lon1))
else:
    lonmin = 175
    lonmax = 185
    lonmin_idx = np.argmin(np.abs(lon - lonmin))
    lonmax_idx = np.argmin(np.abs(lon - lonmax))

# Define variables
var_diff = varname['var_change']
if lon_choice == 'one':
    var_diff = f.variables[var_diff][:,:,:,lon1_idx].squeeze()
else:
    var_diff = f.variables[var_diff][:,:,:,lonmin_idx:lonmax_idx+1].squeeze()
    var_diff = np.ma.average(var_diff, axis=2)
density = f.variables['density'][:]
#Masked values in white
var_diff[var_diff > 4] = 0

# Define variable properties
minmax = varname['minmax']
clevsm = varname['clevsm']
legVar = varname['legVar']
unit = varname['unit']

# density domain
rhomin = 21. ; rhomin_idx = np.argmin(np.abs(density - rhomin))
rhomid = 26. ; rhomid_idx = np.argmin(np.abs(density - rhomid))
rhomax = 28. ; rhomax_idx = np.argmin(np.abs(density - rhomax))


# ------- Build plot variables --------

# Create meshgrid
lat2d, density2d = np.meshgrid(lat, density)

# Levels for shade plot
levels = np.linspace(minmax[0],minmax[1],minmax[2])

# Colormap
#cmap = plt.get_cmap('bwr')
cmap = custom_div_cmap()

# -------- Plot diagram ----------------

fig, ax = plt.subplots(2,1, figsize=(11,8))


# ==== Upper panel ====

cnplot1 = ax[0].contourf(lat2d, density2d, var_diff, cmap=cmap, levels=levels, extend='both')

ax[0].set_ylim([rhomin,rhomid])
ax[0].invert_yaxis()
ax[0].tick_params(
        axis='x',  # changes apply to the x axis
        which='both',  # both major and minor ticks are affected
        bottom='off',  # ticks along the bottom edge are off
        labelbottom='off',
        top='off')

ax[0].axvline(x=0, color='black', ls='--')

# === Lower panel ====

cnplot2 = ax[1].contourf(lat2d, density2d, var_diff, cmap=cmap, levels=levels, extend='both')

ax[1].set_ylim([rhomid,rhomax])
ax[1].invert_yaxis()
ax[1].tick_params(
        axis='x',  # changes apply to the x axis
        which='both',  # both major and minor ticks are affected
        top='off')  # ticks along the bottom edge are off

ax[1].axvline(x=0, color='black', ls='--')

# Re-label x-axis
xlabels=['','60S','40S','20S','0','20N','40N','60N']
ax[1].set_xticklabels(xlabels)


# Remove intersecting tick at rhomid
yticks = ax[1].yaxis.get_major_ticks()
yticks[0].label1.set_visible(False)

# ====

fig.text(0.45, 0.04, 'Latitude', ha='center', fontweight='bold')
fig.text(0.06, 0.5, 'Density', va='center', rotation='vertical', fontweight='bold')

plt.subplots_adjust(hspace=.0001)

cb = plt.colorbar(cnplot2, ax=(ax[0],ax[1]), ticks=levels)
cb.set_label('%s (%s)' % (legVar, unit), fontweight='bold')


if lon_choice == 'one':
    # Converting longitude for title clarity
    if lon1 > 180:
        lon1 = lon1 - 360
    plt.suptitle('%s changes at lon=%d (%s)' %(legVar, lon1, name),
          fontweight='bold', fontsize=14, verticalalignment='top')
else :
    # Converting longitudes for title
    if lonmin>180:
        lonmin = lonmin-360
    if lonmax>180:
        lonmax = lonmax-360
    plt.suptitle('%s changes at longitudes averaged between %d and %d (%s)' % (legVar, lonmin, lonmax, name),
                 fontweight='bold', fontsize=14, verticalalignment='top')

plt.show()


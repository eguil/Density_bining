#!/bin/env python
# -*- coding: utf-8 -*-

"""
Python matplotlib
Make longitude/density diagrams of temperature and salinity change between 2010 and 1945 (or 2000 and 1950),
from observations or CMIP5 mme

"""

# === COMMENTS ===

# There seems to be an issue with the mme_hist file, the figure appeares all white, haven't figured out why
# Works fine for Ishii and Durack&Wijffels data

# =========

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset as open_ncfile
from maps_matplot_lib import defVar, custom_div_cmap, defVarDurack, defVarmme


# ------ Define work and Read files -----------

#name = 'Ishii'
#name = 'EN4'
#name = 'Durack&Wijffels'
name = 'mme_hist'

# Data path
indir = '/data/ericglod/Density_binning/Obs_Prod_density_april16/'
if name=='Ishii':
    file = 'obs.Ishii.historical.r0i0p0.an.ocn.Omon.density.ver-1.latestXCorr.nc'
if name=='EN4':
    file = 'obs.EN4.historical.r0i0p0.mo.ocn.Omon.density.ver-1.latestXCorr.nc'
if name=='Durack&Wijffels':
    file = 'DurackandWijffels_GlobalOceanChanges-NeutralDensity_1950-2000_120209_11_46_11_beta.nc'
if name=='mme_hist':
    indir = '/data/ericglod/Density_binning/Prod_density_april15/Raw/mme_hist/mme/'
    file = 'cmip5.multimodel_All.historical.ensm.an.ocn.Omon.density_3D.nc'

data = indir + file

f = open_ncfile(data,'r')

# Read variables
lat = f.variables['latitude'][:]
lon = f.variables['longitude'][:]

# ------ Define variables ----------

# -- Choose which variable to work on
varname = 'salinity'
#varname = 'temp'


# -- Choose latitude (one specific lat or average a range of lats)
#lat_choice = 'one'
lat_choice = 'range'

if lat_choice == 'one':
    lat1 = -20
    lat1_idx = np.argmin(np.abs(lat - lat1))
else:
    latmin = -25
    latmax = -15
    latmin_idx = np.argmin(np.abs(lat - latmin))
    latmax_idx = np.argmin(np.abs(lat - latmax))
    weights = np.cos(lat*np.pi/180)
    weights = weights[latmin_idx:latmax_idx+1]
    #print(np.shape(weights))
    #print(weights)

# -- Read and build variables
if name == 'Durack&Wijffels':
    varname = defVarDurack(varname)
    var_diff = varname['var_change']
    if lat_choice == 'one':
        var_diff = f.variables[var_diff][:,:,lat1_idx,:].squeeze()
    else:
        var_diff = f.variables[var_diff][:,:,latmin_idx:latmax_idx+1,:].squeeze()
        var_diff = np.ma.average(var_diff, axis=1, weights=weights)
    density = f.variables['density'][:]
    #Masked values in white
    var_diff[var_diff > 4] = 0

elif name == 'mme_hist':
    varname = defVarmme(varname)
    var = varname['var_global']
    if lat_choice == 'one':
        var = f.variables[var][88:,:,lat1_idx,:].squeeze()
    else:
        var = f.variables[var][88:,:,latmin_idx:latmax_idx+1,:]
        var = np.ma.average(var, axis=2, weights=weights)
    # Difference
    var_diff = np.ma.average(var[-5:,:], axis=0) - np.ma.average(var[0:5,:], axis=0)
    density = f.variables['lev'][:]


else:
    varname = defVar(varname)
    var = varname['var']
    if lat_choice == 'one':
        var = f.variables[var][:,:,lat1_idx,:].squeeze()
    else:
        var = f.variables[var][:,:,latmin_idx:latmax_idx+1,:]
        var = np.ma.average(var, axis=2, weights=weights)
    # Difference
    var_diff = np.ma.average(var[-5:,:], axis=0) - np.ma.average(var[0:5,:], axis=0)
    density = f.variables['lev'][:]


# -- Define variable properties
minmax = varname['minmax']
clevsm = varname['clevsm']
legVar = varname['legVar']
unit = varname['unit']

# -- Density domain
rhomin = 21. ; rhomin_idx = np.argmin(np.abs(density - rhomin))
rhomid = 26. ; rhomid_idx = np.argmin(np.abs(density - rhomid))
rhomax = 28. ; rhomax_idx = np.argmin(np.abs(density - rhomax))


# ------- Build plot variables --------

# Move the elements in var_diff so that data starts at 20°E
if name == 'Durack&Wijffels':
    var_diff = np.roll(var_diff,-10,axis=1)
else:
    var_diff = np.roll(var_diff,-200,axis=1)

# Create meshgrid
lon2d, density2d = np.meshgrid(lon, density)

# Levels for shade plot
levels = np.linspace(minmax[0],minmax[1],minmax[2])

# Colormap
cmap = custom_div_cmap()

# -------- Plot diagram ----------------

fig, ax = plt.subplots(2,1, figsize=(12,8))


# ==== Upper panel ====

cnplot1 = ax[0].contourf(lon2d, density2d, var_diff, cmap=cmap, levels=levels, extend='both')

ax[0].set_ylim([rhomin,rhomid])
ax[0].invert_yaxis()
ax[0].tick_params(
        axis='x',  # changes apply to the x axis
        which='both',  # both major and minor ticks are affected
        bottom='off',  # ticks along the bottom edge are off
        labelbottom='off',
        top='off')


# === Lower panel ====

cnplot2 = ax[1].contourf(lon2d, density2d, var_diff, cmap=cmap, levels=levels, extend='both')

ax[1].set_ylim([rhomid,rhomax])
ax[1].invert_yaxis()
ax[1].tick_params(
        axis='x',  # changes apply to the x axis
        which='both',  # both major and minor ticks are affected
        top='off')  # ticks along the bottom edge are off

# Re-label x-axis to match longitudes
if name=='Durack&Wijffels':
    xlocs = [10,40,70,100,130,160,190,220,250,280,310,340]
else:
    xlocs=[-169.5,-139.5,-109.5,-79.5,-49.5,-19.5,10.5,40.5,70.5,100.5,130.5,160.5]
xlabels=['30E','60E','90E','120E','150E','180','150W','120W','90W','60W','30W','0']
ax[1].set_xticks(xlocs)
ax[1].set_xticklabels(xlabels)

# Remove intersecting tick at rhomid
yticks = ax[1].yaxis.get_major_ticks()
yticks[0].label1.set_visible(False)

# ==== Plot settings

fig.text(0.45, 0.04, 'Longitude', ha='center', fontweight='bold')
fig.text(0.06, 0.5, 'Density', va='center', rotation='vertical', fontweight='bold')

plt.subplots_adjust(hspace=.00001)

cb = plt.colorbar(cnplot2, ax=(ax[0],ax[1]), ticks=levels)
cb.set_label('%s (%s)' % (legVar, unit), fontweight='bold')

if lat_choice == 'one':
    plt.suptitle('Longitude/density diagram of %s changes at lat=%d (%s)' %(legVar, lat1, name),
          fontweight='bold', fontsize=14, verticalalignment='top')
else :
    plt.suptitle('Longitude/density diagram of %s changes at lats averaged between %d and %d (%s)' % (legVar, latmin, latmax, name),
                 fontweight='bold', fontsize=14, verticalalignment='top')

plt.show()




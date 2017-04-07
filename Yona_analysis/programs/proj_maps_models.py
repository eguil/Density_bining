#!/bin/env python
# -*- coding: utf-8 -*-

"""
Python matplotlib
Make density projections lonxlat of temperature and salinity change between 2005 and 1950 from mme

"""

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset as open_ncfile
from maps_matplot_lib import defVarmme, proj_map, proj_map_zonal_changes, custom_div_cmap
from modelsDef import defModels

# ------ Read files -----------

# -- Choose single model or mme
name = 'mme_hist'
#name = 'mme_hist_histNat'
#name = 'ens_mean_hist'

if name == 'ens_mean_hist':
    models = defModels()
    model = models[14] # Choose model
    name = model['name']
    nb_members = model['props'][0]

    # Data path
    indirh = '/data/ericglod/Density_binning/Prod_density_april15/Raw/mme_hist/mme/'
    fileh = 'cmip5.' + name + '.historical.ensm.an.ocn.Omon.density.ver-' + model['file_end'] + '.nc'

else:
    indirh = '/data/ericglod/Density_binning/Prod_density_april15/Raw/mme_hist/mme/'
    fileh = 'cmip5.multimodel_All.historical.ensm.an.ocn.Omon.density_3D.nc'
    # for hist - histNat, mme or not ready so use ensemble mean for now
    #fileh = 'cmip5.GFDL-ESM2M.historical.ensm.an.ocn.Omon.density.ver-v20130226.nc'

    if name == 'mme_hist_histNat':
        indirhn = '/data/ericglod/Density_binning/Prod_density_april15/Raw/mme_histNat/mme/'
        filehn = 'cmip5.GFDL-ESM2M.historicalNat.ensm.an.ocn.Omon.density.ver-v20110601.nc'
        datahn = indirhn + filehn
        fhn = open_ncfile(datahn,'r')

datah = indirh + fileh
fh = open_ncfile(datah,'r')

# Read variables
lat = fh.variables['latitude'][:]
lon = fh.variables['longitude'][:]


# ------ Define work and variables ---------

varname = defVarmme('salinity'); v = 'S'
#varname = defVarmme('temp'); v = 'T'

density = np.array([19.0, 19.2, 19.4, 19.6, 19.8, 20.0, 20.2, 20.4, 20.6, 20.8, 21.0, 21.2,
    21.4, 21.6, 21.8, 22.0, 22.2, 22.4, 22.6, 22.8, 23.0, 23.2, 23.4, 23.6,
    23.8, 24.0, 24.2, 24.4, 24.6, 24.8, 25.0, 25.2, 25.4, 25.6, 25.8, 26.0,
    26.1, 26.2, 26.3, 26.4, 26.5, 26.6, 26.7, 26.8, 26.9, 27.0, 27.1, 27.2,
    27.3, 27.4, 27.5, 27.6, 27.7, 27.8, 27.9, 28.0, 28.1, 28.2, 28.3, 28.4, 28.5])


# Define variable properties
var = varname['var_global']
var_std = varname['var_global_std']
minmax = varname['minmax']
clevsm = varname['clevsm']
clevsm_bold = varname['clevsm_bold']
legVar = varname['legVar']
unit = varname['unit']

# Define isopycnal on which we want to project the temperature/salinity/other variable
isopyc1 = 23.5
isopyc2 = 24.5
isopyc3 = 25
isopyc4 = 27.5
isopyc1_idx = np.argmin(np.abs(density - isopyc1))
isopyc4_idx = np.argmin(np.abs(density - isopyc4))

if name == 'mme_hist_histNat':
    varhn = fhn.variables[var][-5:,isopyc1_idx:isopyc4_idx+1,:,:]

var = fh.variables[var][88:,isopyc1_idx:isopyc4_idx+1,:,:] #index 88 = year 1950
#var_std = f.variables[var_std][88:,isopyc1_idx:isopyc4_idx+1,:,:]

sliced_density = density[isopyc1_idx:isopyc4_idx+1]
#print(np.ma.nonzero(var))


# ----- Plot changes and means ---------

cmap = custom_div_cmap()

fig = plt.figure(figsize=(12,18))
ax1 = fig.add_axes([0.05, 0.74, 0.95, 0.2])
ax2 = fig.add_axes([0.05, 0.51, 0.95, 0.2])
ax3 = fig.add_axes([0.05, 0.28, 0.95, 0.2])
ax4 = fig.add_axes([0.05, 0.05, 0.95, 0.2])

if name != 'mme_hist_histNat':
    bmap = proj_map('hist', plt, ax1, minmax, clevsm, clevsm_bold, lat, lon, cmap, isopyc1, sliced_density, var)
    bmap = proj_map('hist', plt, ax2, minmax, clevsm, clevsm_bold, lat, lon, cmap, isopyc2, sliced_density, var)
    bmap = proj_map('hist', plt, ax3, minmax, clevsm, clevsm_bold, lat, lon, cmap, isopyc3, sliced_density, var)
    bmap = proj_map('hist', plt, ax4, minmax, clevsm, clevsm_bold, lat, lon, cmap, isopyc4, sliced_density, var)
else:
    bmap = proj_map('hist-histNat', plt, ax1, minmax, clevsm, clevsm_bold, lat, lon, cmap, isopyc1, sliced_density,
                    var, varhn)
    bmap = proj_map('hist-histNat', plt, ax2, minmax, clevsm, clevsm_bold, lat, lon, cmap, isopyc2, sliced_density,
                    var, varhn)
    bmap = proj_map('hist-histNat', plt, ax3, minmax, clevsm, clevsm_bold, lat, lon, cmap, isopyc3, sliced_density,
                    var, varhn)
    bmap = proj_map('hist-histNat', plt, ax4, minmax, clevsm, clevsm_bold, lat, lon, cmap, isopyc4, sliced_density,
                    var, varhn)


cb = plt.colorbar(bmap[0], ax = (ax1, ax2, ax3, ax4), ticks=bmap[1], orientation='vertical')
cb.set_label('%s (%s)' % (legVar, unit), fontweight='bold')

if name == 'mme_hist':
    plt.suptitle('%s changes (2000-1950), %s ' %(legVar, name), fontweight='bold', fontsize=14, va='top')
    plotName = name + '_' + v + 'changes'
elif name == 'mme_hist_histNat':
    plt.suptitle('%s changes %s (last 5 years)' %(legVar, name), fontweight='bold', fontsize=14, va='top')
    plotName = name + '_' + v + 'changes'
else:
    plt.suptitle('%s changes (2000-1950), %s ensemble mean (%d members)' %(legVar, name, nb_members), fontweight='bold',
                 fontsize=14, va='top')
    plotName = name + '_hist_' + v + 'changes'

plt.show()
#plt.savefig('/home/ysilvy/figures/models/density_projection_maps/hist/'+plotName+'.png')#, bbox_inches='tight')



#!/bin/env python
# -*- coding: utf-8 -*-

"""
Python matplotlib
Make lonxlat maps of temperature and salinity change from the Durack&Wijfells data, projected on specified isopycnals

"""

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset as open_ncfile
from maps_matplot_lib import proj_map, proj_map_zonal_changes, custom_div_cmap, defVarDurack

# ----- Workspace ------

indir = '/data/ericglod/Density_binning/Obs_Prod_density_april16/'
file = 'DurackandWijffels_GlobalOceanChanges-NeutralDensity_1950-2000_120209_11_46_11_beta.nc'
name = 'Durack & Wijffels'

data = indir + file

f = open_ncfile(data,'r')

# ----- Variables ------

# Read variables
density = f.variables['density'][:]
lat = f.variables['latitude'][:]
lon = f.variables['longitude'][:]


varname = defVarDurack('temp'); v='T'
#varname = defVarDurack('salinity'); v='S'


# -- Look at zonal differences or not
zonal_change = 'No'
#zonal_change = 'global' # Zonal mean is global
#zonal_change = 'basin' # Differentiate the basins when building the zonal mean


# Define variable properties
minmax = varname['minmax']
clevsm = varname['clevsm']
clevsm_bold = varname['clevsm_bold']
legVar = varname['legVar']
unit = varname['unit']
var_change = varname['var_change']
var_mean = varname['var_mean']
var_change_er = varname['var_change_er']

# Define isopycnal on which we want to project the temperature/salinity/variable
isopyc1 = 24.0
isopyc2 = 25.0
isopyc3 = 26.7
isopyc4 = 27.5
isopyc1_idx = np.argmin(np.abs(density - isopyc1))
isopyc4_idx = np.argmin(np.abs(density - isopyc4))

var_change = f.variables[var_change][:,slice(isopyc1_idx,isopyc4_idx+1), :, :].squeeze()
var_mean = f.variables[var_mean][:,slice(isopyc1_idx,isopyc4_idx+1), :, :].squeeze()
var_change_er = f.variables[var_change_er][:,slice(isopyc1_idx,isopyc4_idx+1), :, :].squeeze()
sliced_density = density[isopyc1_idx:isopyc4_idx+1]


# ----- Plot changes and means ---------

cmap = custom_div_cmap()


# -- Choose whether to plot the "basic" change, or to remove the zonal mean as well to look at zonal asymmetry

if zonal_change == 'No':
    fig = plt.figure(figsize=(12,20))
    ax1 = fig.add_axes([0.05, 0.74, 0.95, 0.2])
    ax2 = fig.add_axes([0.05, 0.51, 0.95, 0.2])
    ax3 = fig.add_axes([0.05, 0.28, 0.95, 0.2])
    ax4 = fig.add_axes([0.05, 0.05, 0.95, 0.2])

    bmap = proj_map('Durack', plt, ax1, minmax, clevsm, clevsm_bold, lat, lon, cmap, isopyc1, sliced_density,
                    var_change, var_mean, var_change_er)
    bmap = proj_map('Durack', plt, ax2, minmax, clevsm, clevsm_bold, lat, lon, cmap, isopyc2, sliced_density,
                    var_change, var_mean, var_change_er)
    bmap = proj_map('Durack', plt, ax3, minmax, clevsm, clevsm_bold, lat, lon, cmap, isopyc3, sliced_density,
                    var_change, var_mean, var_change_er)
    bmap = proj_map('Durack', plt, ax4, minmax, clevsm, clevsm_bold, lat, lon, cmap, isopyc4, sliced_density,
                    var_change, var_mean, var_change_er)

    cb = plt.colorbar(bmap[0], ax=(ax1,ax2,ax3,ax4),
                          ticks=bmap[1], orientation='vertical')
    cb.set_label('%s (%s)' % (legVar, unit), fontweight='bold')

    plt.suptitle('%s changes (%s)' %(legVar, name), fontweight='bold', fontsize=14)

    #plotName = 'Durack_' + v + 'changes'
    #plt.savefig('/home/ysilvy/figures/obs/density_projection_maps/'+plotName+'.png')


# -- Remove global zonal mean
elif zonal_change == 'global' :

    fig = plt.figure(figsize=(14,14))

    ax11 = fig.add_axes([0.05, 0.74, 0.65, 0.2])
    ax12 = fig.add_axes([0.75, 0.74, 0.15, 0.2])
    ax21 = fig.add_axes([0.05, 0.51, 0.65, 0.2])
    ax22 = fig.add_axes([0.75, 0.51, 0.15, 0.2])
    ax31 = fig.add_axes([0.05, 0.28, 0.65, 0.2])
    ax32 = fig.add_axes([0.75, 0.28, 0.15, 0.2])
    ax41 = fig.add_axes([0.05, 0.05, 0.65, 0.2])
    ax42 = fig.add_axes([0.75, 0.05, 0.15, 0.2])

    bmap = proj_map_zonal_changes('Durack', zonal_change, plt, ax11, ax12, minmax, clevsm, lat, lon, cmap, isopyc1, sliced_density, var_change, var_mean)
    bmap = proj_map_zonal_changes('Durack', zonal_change, plt, ax21, ax22, minmax, clevsm, lat, lon, cmap, isopyc2, sliced_density, var_change, var_mean)
    bmap = proj_map_zonal_changes('Durack', zonal_change, plt, ax31, ax32, minmax, clevsm, lat, lon, cmap, isopyc3, sliced_density, var_change, var_mean)
    bmap = proj_map_zonal_changes('Durack', zonal_change, plt, ax41, ax42, minmax, clevsm, lat, lon, cmap, isopyc4, sliced_density, var_change, var_mean)

    cb = plt.colorbar(bmap[0], ax=(ax11,ax12,ax21,ax22,ax31,ax32,ax41,ax42),
                      ticks=bmap[1], orientation='vertical', pad=0.05)
    cb.set_label('%s (%s)' % (legVar, unit), fontweight='bold')


    plt.suptitle('%s zonal asymmetry (map) and zonal mean of %s changes' %(name, legVar),
                 fontweight='bold', fontsize=13)


# -- Remove zonal mean for each basin
else :

    fig = plt.figure(figsize=(14,14))

    ax11 = fig.add_axes([0.05, 0.74, 0.65, 0.2])
    ax12 = fig.add_axes([0.75, 0.74, 0.15, 0.2])
    ax21 = fig.add_axes([0.05, 0.51, 0.65, 0.2])
    ax22 = fig.add_axes([0.75, 0.51, 0.15, 0.2])
    ax31 = fig.add_axes([0.05, 0.28, 0.65, 0.2])
    ax32 = fig.add_axes([0.75, 0.28, 0.15, 0.2])
    ax41 = fig.add_axes([0.05, 0.05, 0.65, 0.2])
    ax42 = fig.add_axes([0.75, 0.05, 0.15, 0.2])

    bmap = proj_map_zonal_changes('Durack', zonal_change, plt, ax11, ax12, minmax, clevsm, lat, lon, cmap, isopyc1,
                                  sliced_density, var_change, var_mean)
    bmap = proj_map_zonal_changes('Durack', zonal_change, plt, ax21, ax22, minmax, clevsm, lat, lon, cmap, isopyc2,
                                  sliced_density, var_change, var_mean)
    bmap = proj_map_zonal_changes('Durack', zonal_change, plt, ax31, ax32, minmax, clevsm, lat, lon, cmap, isopyc3,
                                  sliced_density, var_change, var_mean)
    bmap = proj_map_zonal_changes('Durack', zonal_change, plt, ax41, ax42, minmax, clevsm, lat, lon, cmap, isopyc4,
                                  sliced_density, var_change, var_mean)

    cb = plt.colorbar(bmap[0], ax=(ax11, ax12, ax21, ax22, ax31, ax32, ax41, ax42),
                      ticks=bmap[1], orientation='vertical', pad=0.05)
    cb.set_label('%s (%s)' % (legVar, unit), fontweight='bold')

    plt.suptitle(
        '%s zonal basin asymmetry (map) and global zonal mean (plot) of %s changes' % (name, legVar),
        fontweight='bold', fontsize=13)

plt.show()





#!/bin/env python
# -*- coding: utf-8 -*-

"""
Python matplotlib

"""

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset as open_ncfile
from maps_matplot_lib import defVarDurack, zonal_2D, custom_div_cmap

# ----- Workspace ------

indir = '/data/ericglod/Density_binning/Obs_Prod_density_april16/'
file = 'DurackandWijffels_GlobalOceanChanges-NeutralDensity_1950-2000_120209_11_46_11_beta.nc'
name = 'Durack & Wijffels'

data = indir + file

f = open_ncfile(data,'r')

# ----- Variables ------

# -- Read variables
density = f.variables['density'][:]
lat = f.variables['latitude'][:]

varname = defVarDurack('salinity')
#varname = defVarDurack('temp')

# -- Define variable properties
minmax = varname['minmax_zonal']
clevsm = varname['clevsm_zonal']
clevsm_bold = varname['clevsm_bold']
legVar = varname['legVar']
unit = varname['unit']

var_mean = varname['var_mean_zonal']
var_change = varname['var_change_zonal']
var_change_er = varname['var_change_zonal_er']

var_attributes = f.variables[var_mean]
var_mean = f.variables[var_mean][:].squeeze()
var_change = f.variables[var_change][:].squeeze()
var_change_er = f.variables[var_change_er][:].squeeze()

valmask = var_attributes.missing_value
latN = len(lat)
levN = len(density)
basinN = 4

basin = ['Global', 'Pacific', 'Atlantic', 'Indian']

# -- Determine field in 1950 and 2000 from mean field
var_1950 = var_mean - var_change/2
var_2000 = var_mean + var_change/2

# -- Density domain
rhomin = 21.
rhomid = 26.
rhomax = 28.
domrho = [rhomin, rhomid, rhomax]

# -- Build plot variables for total change
var_change_p = var_change[:,:,1].squeeze()
var_change_a = var_change[:,:,2].squeeze()
var_change_i = var_change[:,:,3].squeeze()
var_change_er_p = var_change_er[:,:,1].squeeze()
var_change_er_a = var_change_er[:,:,2].squeeze()
var_change_er_i = var_change_er[:,:,3].squeeze()
var_mean_p = var_mean[:,:,1].squeeze()
var_mean_a = var_mean[:,:,2].squeeze()
var_mean_i = var_mean[:,:,3].squeeze()


# ------------------------------------------------------------------------------
#                               Get bowl position
# ------------------------------------------------------------------------------

bowl = np.ma.masked_all((latN,4))
bowl_1950 = np.ma.masked_all((latN,4))
bowl_2000 = np.ma.masked_all((latN,4))


for ibasin in range(1,4):
    for ilat in range(latN):
        nomask = np.ma.flatnotmasked_edges(var_mean[:,ilat,ibasin]) #returns indices of the 1st and last unmasked values
        nomask_1950 = np.ma.flatnotmasked_edges(var_1950[:, ilat, ibasin])
        nomask_2000 = np.ma.flatnotmasked_edges(var_2000[:, ilat, ibasin])

        if nomask != None :
            bowl[ilat,ibasin] = density[nomask[0]]
        else :
            bowl[ilat,ibasin] = np.ma.masked

        if nomask_1950 != None:
            bowl_1950[ilat, ibasin] = density[nomask_1950[0]]
        else:
            bowl_1950[ilat, ibasin] = np.ma.masked

        if nomask_2000 != None:
            bowl_2000[ilat, ibasin] = density[nomask_2000[0]]
        else:
            bowl_2000[ilat, ibasin] = np.ma.masked


bowl_p = bowl[:,1]
bowl_a = bowl[:,2]
bowl_i = bowl[:,3]



# ------------------------------------------------------------------------------
#     Partial derivative of mean field for latitude dS/dy --> dvar_dy
# ------------------------------------------------------------------------------

# -- Initialize dvar_dy
spiciness_mean = np.ma.masked_all(np.shape(var_mean))
delta_lat = 1 # degrees
# Roll latitude
var_mean_plus1lat = np.roll(var_mean, -1, axis=1)
var_1950_plus1lat = np.roll(var_1950, -1, axis=1)
var_2000_plus1lat = np.roll(var_2000, -1, axis=1)

# -- Determine derivative dvar_dy (mean)
spiciness_mean = (var_mean_plus1lat - var_mean)/delta_lat
spiciness_mean[:,-1,:] = np.ma.masked
# Now 1950
spiciness_1950 = (var_1950_plus1lat - var_1950)/delta_lat
spiciness_1950[:,-1,:] = np.ma.masked
# 2000
spiciness_2000 = (var_2000_plus1lat - var_2000)/delta_lat
spiciness_2000[:,-1,:] = np.ma.masked
# Spiciness change
spiciness_change = spiciness_2000 - spiciness_1950

# -- Separate for each basin
spiciness_mean_p = spiciness_mean[:,:,1]
spiciness_mean_a = spiciness_mean[:,:,2]
spiciness_mean_i = spiciness_mean[:,:,3]


# ------------------------------------------------------------------------------
#                               Now plot
# ------------------------------------------------------------------------------

# -- Create variable bundles
varPac = {'name': 'Pacific', 'var_change': var_change_p, 'var_mean': var_mean_p, 'var_error': var_change_er_p,
          'var_1950': var_1950[:,:,1],  'var_2000': var_2000[:,:,1],
          'dvar_dy': spiciness_mean_p, 'spiciness_change': spiciness_change[:,:,1], 'bowl': bowl_p,
          'dvar_dy_1950': spiciness_1950[:,:,1], 'dvar_dy_2000': spiciness_2000[:,:,1]}
varAtl = {'name': 'Atlantic', 'var_change': var_change_a, 'var_mean': var_mean_a, 'var_error': var_change_er_a,
            'var_1950': var_1950[:,:,2],  'var_2000': var_2000[:,:,2],
          'dvar_dy': spiciness_mean_a, 'spiciness_change': spiciness_change[:,:,2], 'bowl': bowl_a,
          'dvar_dy_1950': spiciness_1950[:,:,2], 'dvar_dy_2000': spiciness_2000[:,:,2]}
varInd = {'name': 'Indian', 'var_change': var_change_i, 'var_mean': var_mean_i, 'var_error': var_change_er_i,
            'var_1950': var_1950[:,:,3],  'var_2000': var_2000[:,:,3],
            'dvar_dy': spiciness_mean_i, 'spiciness_change': spiciness_change[:,:,3], 'bowl': bowl_i,
          'dvar_dy_1950': spiciness_1950[:,:,3], 'dvar_dy_2000': spiciness_2000[:,:,3]}




# ==== dS/dy ====

fig1, axes = plt.subplots(nrows=2, ncols=3, figsize=(17,5))

levels = np.linspace(-0.1,0.1,16)
cmap = plt.get_cmap('bwr')

cnplot1 = zonal_2D(plt, 'dvar_dy', axes[0,0], axes[1,0], 'left', lat, density, varAtl,
                   domrho, cmap, levels)

cnplot1 = zonal_2D(plt, 'dvar_dy', axes[0,1], axes[1,1], 'mid', lat, density, varPac,
                   domrho, cmap, levels)

cnplot1 = zonal_2D(plt, 'dvar_dy', axes[0,2], axes[1,2], 'right', lat, density, varInd,
                   domrho, cmap, levels)


plt.subplots_adjust(hspace=.0001, wspace=0.05, left=0.04, right=0.86)

cb = plt.colorbar(cnplot1, ax=axes.ravel().tolist(), ticks=levels[::3])

plt.suptitle('dS/dy (%s)' %(name,),
          fontweight='bold', fontsize=14, verticalalignment='top')

#plt.close()

# ==== dS/dy change ====

fig2, axes = plt.subplots(nrows=2, ncols=3, figsize=(17,5))

levels = np.linspace(-0.1,0.1,16)
cmap = custom_div_cmap()

cnplot2 = zonal_2D(plt, 'spiciness_change', axes[0,0], axes[1,0], 'left', lat, density, varAtl,
                   domrho, cmap, levels)

cnplot2 = zonal_2D(plt, 'spiciness_change', axes[0,1], axes[1,1], 'mid', lat, density, varPac,
                   domrho, cmap, levels)

cnplot2 = zonal_2D(plt, 'spiciness_change', axes[0,2], axes[1,2], 'right', lat, density, varInd,
                   domrho, cmap, levels)


plt.subplots_adjust(hspace=.0001, wspace=0.05, left=0.04, right=0.86)

cb = plt.colorbar(cnplot1, ax=axes.ravel().tolist(), ticks=levels[::3])

plt.suptitle('dS/dy change (%s)' %(name,),
          fontweight='bold', fontsize=14, verticalalignment='top')

#plt.close()

# ==== Mean dS/dy contours 1950/2000 ====

fig3, axes = plt.subplots(nrows=2, ncols=3, figsize=(17,5))

levels = np.arange(-0.1,0.101,0.05)
cmap = None

zonal_2D(plt, 'spiciness_fields', axes[0,0], axes[1,0], 'left', lat, density, varAtl,
         domrho, cmap, levels, clevsm, clevsm_bold)

zonal_2D(plt, 'spiciness_fields', axes[0,1], axes[1,1], 'mid', lat, density, varPac,
         domrho, cmap, levels, clevsm, clevsm_bold)

zonal_2D(plt, 'spiciness_fields', axes[0,2], axes[1,2], 'right', lat, density, varInd,
         domrho, cmap, levels, clevsm, clevsm_bold)


plt.subplots_adjust(hspace=.0001, wspace=0.05, left=0.04, right=0.86)

plt.suptitle('dS/dy contours in 1950 and in 2000',
          fontweight='bold', fontsize=14, verticalalignment='top')


plt.show()
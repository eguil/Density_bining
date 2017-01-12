#!/bin/env python
# -*- coding: utf-8 -*-

"""
Python matplotlib
Compute salinity changes due to isopycnal migration and plot density/latitude diagrams for each basin

"""

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset as open_ncfile
from maps_matplot_lib import defVarmme, zonal_2D
from scipy import interpolate

# ----- Workspace ------

indir = '/data/ericglod/Density_binning/Prod_density_april15/mme_hist/'
file_2D = 'cmip5.multimodel_Nat.historical.ensm.an.ocn.Omon.density_zon2D.nc'
file_1D = 'cmip5.multimodel_Nat.historical.ensm.an.ocn.Omon.density_zon1D.nc'
name = 'mme_hist'

data_2D = indir + file_2D
data_1D = indir + file_1D

f2D = open_ncfile(data_2D,'r')
f1D = open_ncfile(data_1D,'r')

# ----- Variables ------

# -- Read variables
lat = f2D.variables['latitude'][:]
density = np.array([19.0, 19.2, 19.4, 19.6, 19.8, 20.0, 20.2, 20.4, 20.6, 20.8, 21.0, 21.2,
    21.4, 21.6, 21.8, 22.0, 22.2, 22.4, 22.6, 22.8, 23.0, 23.2, 23.4, 23.6,
    23.8, 24.0, 24.2, 24.4, 24.6, 24.8, 25.0, 25.2, 25.4, 25.6, 25.8, 26.0,
    26.1, 26.2, 26.3, 26.4, 26.5, 26.6, 26.7, 26.8, 26.9, 27.0, 27.1, 27.2,
    27.3, 27.4, 27.5, 27.6, 27.7, 27.8, 27.9, 28.0, 28.1, 28.2, 28.3, 28.4, 28.5])

varname = defVarmme('salinity')
#varname = defVarmme('temp')

# -- Define variable properties
minmax = varname['minmax_zonal']
clevsm = varname['clevsm']
clevsm_bold = varname['clevsm_bold']
legVar = varname['legVar']
unit = varname['unit']

var = varname['var_zonal']

var_attributes = f2D.variables[var]
var = f2D.variables[var][89:,:,:,:]

# -- Determine mean field and change
var_mean = np.ma.average(var,axis=0)
var_change = np.ma.average(var[-5:,:,:,:],axis=0) -np.ma.average(var[0:5,:,:,:],axis=0)

# -- Determine field in 1950 and 2000
var_1950 = np.ma.average(var[0:5,:,:,:], axis=0)
var_2000 = np.ma.average(var[-5:,:,:,:],axis=0)

valmask = var_attributes.missing_value
latN = len(lat)
levN = len(density)
basinN = 4

basin = ['Global', 'Atlantic', 'Pacific', 'Indian']

delta_t = 50 # years

# -- Density domain
rhomin = 21.
rhomid = 26.
rhomax = 28.
domrho = [rhomin, rhomid, rhomax]

# -- Build plot variables for total change
var_change_p = var_change[1,:,:].squeeze()
var_change_a = var_change[2,:,:].squeeze()
var_change_i = var_change[3,:,:].squeeze()
var_mean_p = var_mean[1,:,:].squeeze()
var_mean_a = var_mean[2,:,:].squeeze()
var_mean_i = var_mean[3,:,:].squeeze()

# -- Bowl position
bowl = f1D.variables['ptopsigma'][89:,:,:]
bowl_1950 = np.ma.average(bowl[0:5,:,:], axis=0)
bowl_2000 = np.ma.average(bowl[-5:,:,:], axis=0)
bowl = np.ma.average(bowl, axis=0)

bowl_p = bowl[2,:]
bowl_a = bowl[1,:]
bowl_i = bowl[3,:]


# ------------------------------------------------------------------------------
#       Determine the different terms of the migration-driven change
# ------------------------------------------------------------------------------


# ==== Partial derivative of mean field for density dS/dsigma --> dvar_dsig ===

# -- Initialize dvar_dsig
dvar_dsig = np.ma.masked_all(np.shape(var_mean))

# -- Roll arrays
var_mean_plus1sig = np.roll(var_mean, -1, axis=1)
density1 = np.roll(density, -1)
# -- Substract density arrays to get delta sigma and transform into 2d array
delta_sig_1d = density1-density
delta_sig = np.tile(delta_sig_1d,(latN,1)); delta_sig = np.transpose(delta_sig)
delta_sig = np.tile(delta_sig,(basinN,1,1))

# -- Determine derivative dvar/disgma
dvar_dsig = (var_mean_plus1sig - var_mean)/delta_sig
dvar_dsig[:,-1,:] = np.ma.masked

# -- Separate for each basin
dvar_dsig_p = dvar_dsig[2,:,:]
dvar_dsig_a = dvar_dsig[1,:,:]
dvar_dsig_i = dvar_dsig[3,:,:]


# ==== Partial derivative of mean field for latitude dS/dy --> dvar_dy ===

# -- Initialize dvar_dy
dvar_dy = np.ma.masked_all(np.shape(var_mean))
delta_lat = 1 # degrees
# Roll latitude
var_mean_plus1lat = np.roll(var_mean, -1, axis=2)

# -- Determine derivative dvar_dy
dvar_dy = (var_mean_plus1lat - var_mean)/delta_lat
dvar_dy[:,:,-1] = np.ma.masked

# -- Separate for each basin
dvar_dy_p = dvar_dy[2,:,:]
dvar_dy_a = dvar_dy[1,:,:]
dvar_dy_i = dvar_dy[3,:,:]



# ==== Partial derivative of latitude for time dy/dt --> dy_dt ===

# -- Define max distance (in lat degrees) between same var_1950 and var_2000
max_lat_shift = 5

# -- Initialize dy/dt to masked value
dy_dt = np.ma.masked_all((basinN,levN,latN))

# -- HR latitude
step = 1
lat_hr = np.arange(-89.5,89.501,step)

# -- Start calculation
for ibasin in range(1,4):

    for isig in range(levN):

        if np.any(var_1950[ibasin,isig,:] != valmask):

            # -- Interpolate on a finer latitude grid
            z = var_2000[ibasin,isig, :].squeeze()
            # Find contiguous unmasked data of var_2000 at isig to interpolate
            slice_z = np.ma.flatnotmasked_contiguous(z)
            if len(slice_z) != 0:
                # print(ibasin, isig, len(slice_z))
                var_2000_hr = np.ma.masked_all(len(lat_hr))
                for j in range(len(slice_z)):
                    #print(basin[ibasin], density[isig])
                    lat1 = lat[slice_z[j]]
                    z1 = z[slice_z[j]]
                    if len(lat1) > 2:
                        #print(lat1[0],lat1[-1],z1[0],z1[-1])
                        lat1_hr = np.arange(lat1[0], lat1[-1] + step/2, step)
                        #print(len(lat1_hr))
                        interp = interpolate.interp1d(lat1, z1, bounds_error=False)
                        var_2000_hr1 = interp(lat1_hr)
                        i1 = np.searchsorted(lat_hr, lat1_hr[0])
                        i2 = np.searchsorted(lat_hr, lat1_hr[-1])
                        # # Clean up exceptions
                        # if (69.999 <= lat1_hr[-1] <= 70.0001) :
                        #     var_2000_hr1 = var_2000_hr1[:-1]
                        #     print(np.shape(var_2000_hr1))
                        var_2000_hr[i1:i2 + 1] = var_2000_hr1
                        # print(lat1_hr[0], lat1_hr[-1], var_2000_hr[i1], var_2000_hr[i2], len(var_2000_hr))


                for ilat in range(latN):

                    if (var_1950[ibasin,isig,ilat] != valmask):

                        # Look at the slope of S(y) : - if dS/dy is close to zero, only look at the vertical shift
                        if np.ma.abs(dvar_dy[ibasin,isig,ilat]) <= 0.01 :
                            dy_dt[ibasin,isig,ilat] = 0
                        else:

                            # Find latitude indices where abs(lat-lat[ilat])< threshold
                            ilat_interval = np.flatnonzero(np.ma.abs(lat_hr[:] - lat[ilat]) < max_lat_shift)

                            if len(ilat_interval) != 0:
                                lat_interval = lat_hr[ilat_interval]
                                # print('Hello')
                                # print(ibasin,density[isig],lat[ilat],var_1950[ibasin,isig,ilat])
                                # print(lat_interval)
                                # print(var_2000[ibasin,isig,ilat_interval])

                                ilat2 = ilat_interval[np.ma.argmin(np.ma.abs(var_2000_hr[ilat_interval] - var_1950[ibasin,isig,ilat]))]
                                lat2 = lat_hr[ilat2]
                                # Save value to plot
                                # delta_var_min[isig,ilat,ibasin] = np.ma.abs(var_2000_hr[ilat2]-var_1950[ibasin,isig,ilat])
                                # print(delta_var_min[isig,ilat,ibasin])
                                # print(lat2, lat2 - lat[ilat])
                                dy_dt[ibasin, isig, ilat] = -(lat2 - lat[ilat])


# -- Separate for each basin
dy_dt_p = dy_dt[2,:,:]
dy_dt_a = dy_dt[1,:,:]
dy_dt_i = dy_dt[3,:,:]



# ==== Partial derivative of sigma for time dsigma/dt ===

# -- Define max distance (in density units) between same var_1950 and var_2000
max_density_shift1 = 0.25
max_density_shift2 = 0.25

# -- Initialize dsigma/dt to masked value
dsig_dt = np.ma.masked_all((basinN,levN,latN))

# -- HR density
step = 0.005
sig_hr1 = np.arange(19,26,step)
sig_hr2 = np.arange(26,28.501,step/2)
density_hr = np.hstack((sig_hr1,sig_hr2))

# -- Start calculation

for ibasin in range(1,4):

    for ilat in range(latN):

        if np.any(var_1950[ibasin,:,ilat] != valmask):

            # -- Interpolate on a finer density grid
            z = var_2000[ibasin,:, ilat].squeeze()
            # Find contiguous unmasked data of var_2000 at ilat to interpolate
            slice_z = np.ma.flatnotmasked_contiguous(z)
            if len(slice_z) != 0:
                # print(ibasin, ilat, len(slice_z))
                var_2000_hr = np.ma.masked_all(len(density_hr))
                density1 = density[slice_z[0]]
                z1 = z[slice_z[0]]
                if len(density1) > 2:
                    # print(density1[0],density1[-1],z1[0],z1[-1])
                    i1 = np.searchsorted(density_hr, density1[0])
                    i2 = np.searchsorted(density_hr, density1[-1])
                    density1_hr = density_hr[i1:i2 + 1]
                    # print(len(density1_hr))
                    interp = interpolate.interp1d(density1, z1, bounds_error=False)
                    var_2000_hr1 = interp(density1_hr)
                    var_2000_hr[i1:i2 + 1] = var_2000_hr1
                    # print(density1_hr[0], density1_hr[-1], var_2000_hr[i1], var_2000_hr[i2], len(var_2000_hr))

                for isig in range(levN):

                    if density[isig] < 26.5 :
                        max_density_shift = max_density_shift1
                    else :
                        max_density_shift = max_density_shift2

                    if (var_1950[ibasin,isig,ilat] != valmask) :

                        # Look at the slope of dS/dsigma - if it is close to zero, only look at the horizontal shift
                        if np.ma.abs(dvar_dsig[ibasin,isig,ilat]) <= 0.25:
                            dsig_dt[ibasin,isig,ilat] = 0
                        else:

                            # Find density indices where abs(density-density[isig])<0.5 kg.m-3
                            isig_interval = np.flatnonzero(np.ma.abs(density_hr[:] - density[isig]) < max_density_shift)

                            if len(isig_interval) != 0:
                                sig_interval = density_hr[isig_interval]
                                # print(ibasin,density[isig],lat[ilat],var_1950[ibasin,isig,ilat])
                                # print(sig_interval)
                                # print(var_2000[ibasin,isig_interval,ilat])

                                isig2 = isig_interval[np.ma.argmin(np.ma.abs(var_2000_hr[isig_interval] - var_1950[ibasin,isig,ilat]))]
                                sig2 = density_hr[isig2]
                                # Save value to plot
                                # delta_var_min[ibasin,isig,ilat] = np.ma.abs(var_2000[ibasin,isig2,ilat]-var_1950[ibasin,isig,ilat])
                                # print(delta_var_min[ibasin,isig,ilat])
                                dsig_dt[ibasin, isig, ilat] = -(sig2 - density[isig])

# -- Separate for each basin
dsig_dt_p = dsig_dt[2,:,:]
dsig_dt_a = dsig_dt[1,:,:]
dsig_dt_i = dsig_dt[3,:,:]



# ------------------------------------------------------------------------------
#       Determine residual change : total change - migration-driven change
# ------------------------------------------------------------------------------

var_change_isopmig_p = dvar_dy_p*dy_dt_p + dvar_dsig_p*dsig_dt_p
var_change_res_p = var_change_p - var_change_isopmig_p
var_change_isopmig_a = dvar_dy_a*dy_dt_a + dvar_dsig_a*dsig_dt_a
var_change_res_a = var_change_a - var_change_isopmig_a
var_change_isopmig_i = dvar_dy_i*dy_dt_i + dvar_dsig_i*dsig_dt_i
var_change_res_i = var_change_i - var_change_isopmig_i


# ------------------------------------------------------------------------------
#                               Now plot
# ------------------------------------------------------------------------------

# -- Create variable bundles
varPac = {'name': 'Pacific', 'var_change': var_change_p, 'var_mean': var_mean_p,
          'var_1950': var_1950[2,:,:],  'var_2000': var_2000[2,:,:],
          'var_change_res': var_change_res_p, 'dvar_dsig': dvar_dsig_p, 'dvar_dy': dvar_dy_p,
          'dy_dt': dy_dt_p, 'dsig_dt': dsig_dt_p, 'bowl': bowl_p,}
varAtl = {'name': 'Atlantic', 'var_change': var_change_a, 'var_mean': var_mean_a,
            'var_1950': var_1950[1,:,:],  'var_2000': var_2000[1,:,:],
          'var_change_res': var_change_res_a, 'dvar_dsig': dvar_dsig_a, 'dvar_dy': dvar_dy_a,
          'dy_dt': dy_dt_a, 'dsig_dt': dsig_dt_a, 'bowl': bowl_a,}
varInd = {'name': 'Indian', 'var_change': var_change_i, 'var_mean': var_mean_i,
            'var_1950': var_1950[3,:,:],  'var_2000': var_2000[3,:,:],
          'var_change_res': var_change_res_i, 'dvar_dsig': dvar_dsig_i, 'dvar_dy': dvar_dy_i,
          'dy_dt': dy_dt_i, 'dsig_dt': dsig_dt_i, 'bowl': bowl_i,}



# ==== dS/dsigma ====

# fig6, axes = plt.subplots(nrows=2, ncols=3, figsize=(17,5))
#
# cnplot6 = zonal_2D(plt, 'dvar_dsig', axes[0,0], axes[1,0], 'left', lat, density, varAtl, minmax, domrho)
#
# cnplot6 = zonal_2D(plt, 'dvar_dsig', axes[0,1], axes[1,1], 'mid', lat, density, varPac, minmax, domrho)
#
# cnplot6 = zonal_2D(plt, 'dvar_dsig', axes[0,2], axes[1,2], 'right', lat, density, varInd, minmax, domrho)
#
#
# plt.subplots_adjust(hspace=.0001, wspace=0.05, left=0.04, right=0.86)
#
# cb = plt.colorbar(cnplot6[0], ax=axes.ravel().tolist())
#
# plt.suptitle('dS/dsigma (%s)' %(name,),
#           fontweight='bold', fontsize=14, verticalalignment='top')


# ==== dS/dy ====

# fig9, axes = plt.subplots(nrows=2, ncols=3, figsize=(17,5))
#
# cnplot9 = zonal_2D(plt, 'dvar_dy', axes[0,0], axes[1,0], 'left', lat, density, varAtl, minmax, domrho)
#
# cnplot9 = zonal_2D(plt, 'dvar_dy', axes[0,1], axes[1,1], 'mid', lat, density, varPac, minmax, domrho)
#
# cnplot9 = zonal_2D(plt, 'dvar_dy', axes[0,2], axes[1,2], 'right', lat, density, varInd, minmax, domrho)
#
#
# plt.subplots_adjust(hspace=.0001, wspace=0.05, left=0.04, right=0.86)
#
# cb = plt.colorbar(cnplot9[0], ax=axes.ravel().tolist())
#
# plt.suptitle('dS/dy (%s)' %(name,),
#           fontweight='bold', fontsize=14, verticalalignment='top')


# ==== dsigma/dt ====

fig7, axes = plt.subplots(nrows=2, ncols=3, figsize=(17,5))

cnplot7 = zonal_2D(plt, 'dsig_dt', axes[0,0], axes[1,0], 'left', lat, density, varAtl, minmax, domrho)

cnplot7 = zonal_2D(plt, 'dsig_dt', axes[0,1], axes[1,1], 'mid', lat, density, varPac, minmax, domrho)

cnplot7 = zonal_2D(plt, 'dsig_dt', axes[0,2], axes[1,2], 'right', lat, density, varInd, minmax, domrho)


plt.subplots_adjust(hspace=.0001, wspace=0.05, left=0.04, right=0.86)

cb = plt.colorbar(cnplot7[0], ax=axes.ravel().tolist(), ticks=cnplot7[1][::3])

plt.suptitle('dsigma/dt (%s)' %(name,),
          fontweight='bold', fontsize=14, verticalalignment='top')


# ==== dy/dt ====

fig8, axes = plt.subplots(nrows=2, ncols=3, figsize=(17,5))

cnplot8 = zonal_2D(plt, 'dy_dt', axes[0,0], axes[1,0], 'left', lat, density, varAtl, minmax, domrho)

cnplot8 = zonal_2D(plt, 'dy_dt', axes[0,1], axes[1,1], 'mid', lat, density, varPac, minmax, domrho)

cnplot8 = zonal_2D(plt, 'dy_dt', axes[0,2], axes[1,2], 'right', lat, density, varInd, minmax, domrho)


plt.subplots_adjust(hspace=.0001, wspace=0.05, left=0.04, right=0.86)

cb = plt.colorbar(cnplot8[0], ax=axes.ravel().tolist(), ticks=cnplot8[1][::3])

plt.suptitle('dy/dt (%s)' %(name,),
          fontweight='bold', fontsize=14, verticalalignment='top')


# ==== Density driven term ====

fig5, axes = plt.subplots(nrows=2, ncols=3, figsize=(17,5))

cnplot5 = zonal_2D(plt, 'isopyc_mig_sig', axes[0,0], axes[1,0], 'left', lat, density, varAtl, minmax, domrho)

cnplot5 = zonal_2D(plt, 'isopyc_mig_sig', axes[0,1], axes[1,1], 'mid', lat, density, varPac, minmax, domrho)

cnplot5 = zonal_2D(plt, 'isopyc_mig_sig', axes[0,2], axes[1,2], 'right', lat, density, varInd, minmax, domrho)


plt.subplots_adjust(hspace=.0001, wspace=0.05, left=0.04, right=0.86)

cb = plt.colorbar(cnplot5[0], ax=axes.ravel().tolist(), ticks=cnplot5[1][::3])
cb.set_label('%s (%s)' % (legVar, unit), fontweight='bold')

plt.suptitle('%s changes due to density-driven isopycnal migration (%s)' %(legVar, name),
         fontweight='bold', fontsize=14, verticalalignment='top')


# ==== Latitude driven term ====

fig4, axes = plt.subplots(nrows=2, ncols=3, figsize=(17,5))

cnplot4 = zonal_2D(plt, 'isopyc_mig_lat', axes[0,0], axes[1,0], 'left', lat, density, varAtl, minmax, domrho)

cnplot4 = zonal_2D(plt, 'isopyc_mig_lat', axes[0,1], axes[1,1], 'mid', lat, density, varPac, minmax, domrho)

cnplot4 = zonal_2D(plt, 'isopyc_mig_lat', axes[0,2], axes[1,2], 'right', lat, density, varInd, minmax, domrho)


plt.subplots_adjust(hspace=.0001, wspace=0.05, left=0.04, right=0.86)

cb = plt.colorbar(cnplot4[0], ax=axes.ravel().tolist(), ticks=cnplot4[1][::3])
cb.set_label('%s (%s)' % (legVar, unit), fontweight='bold')

plt.suptitle('%s changes due to latitude-driven isopycnal migration (%s)' %(legVar, name),
          fontweight='bold', fontsize=14, verticalalignment='top')


# ==== Mean contours 1950/2000 ====

fig14, axes = plt.subplots(nrows=2, ncols=3, figsize=(17,5))

zonal_2D(plt, 'mean_fields', axes[0,0], axes[1,0], 'left', lat, density, varAtl, minmax, domrho, clevsm)

zonal_2D(plt, 'mean_fields', axes[0,1], axes[1,1], 'mid', lat, density, varPac, minmax, domrho, clevsm)

zonal_2D(plt, 'mean_fields', axes[0,2], axes[1,2], 'right', lat, density, varInd, minmax, domrho, clevsm)


plt.subplots_adjust(hspace=.0001, wspace=0.05, left=0.04, right=0.86)

plt.suptitle('Mean %s field in 1950 and in 2000' %(legVar,),
          fontweight='bold', fontsize=14, verticalalignment='top')


plt.show()
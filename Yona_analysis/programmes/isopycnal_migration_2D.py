#!/bin/env python
# -*- coding: utf-8 -*-

"""
Python matplotlib
Compute salinity changes due to isopycnal migration and plot density/latitude diagrams for each basin

"""

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset as open_ncfile
from maps_matplot_lib import defVarDurack, zonal_2D
from scipy import interpolate

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

delta_t = 50 # years

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
            #print(ibasin, ilat, nomask[0], density[nomask[0]])
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
#       Determine the different terms of the migration-driven change
# ------------------------------------------------------------------------------


# ==== Partial derivative of mean field for density dS/dsigma --> dvar_dsig ===

# -- Initialize dvar_dsig
dvar_dsig = np.ma.ones(np.shape(var_mean))*valmask
#print('dvar_dsig shape: ', np.shape(dvar_dsig))
# -- Roll arrays
var_mean_plus1sig = np.roll(var_mean, -1, axis=0)
density1 = np.roll(density, -1)
# -- Substract density arrays to get delta sigma and transform into 2d array
delta_sig_1d = density1-density
delta_sig = np.tile(delta_sig_1d,(latN,1)); delta_sig = np.transpose(delta_sig)
delta_sig = np.tile(delta_sig,(basinN,1,1)); delta_sig = np.transpose(delta_sig, (1,2,0))
#print(np.shape(delta_sig))
#print(delta_sig[15,23,0],delta_sig[15,23,1],delta_sig[15,23,2],delta_sig[15,23,3])

# -- Determine derivative dvar/disgma
dvar_dsig = (var_mean_plus1sig - var_mean)/delta_sig
dvar_dsig[-1,:] = valmask
#print(np.shape(dvar_dsig))

# -- Mask
dvar_dsig[np.nonzero(dvar_dsig == valmask)] = np.ma.masked

# -- Separate for each basin
dvar_dsig_p = dvar_dsig[:,:,1]
dvar_dsig_a = dvar_dsig[:,:,2]
dvar_dsig_i = dvar_dsig[:,:,3]



# ==== Partial derivative of mean field for latitude dS/dy --> dvar_dy ===

# -- Initialize dvar_dy
dvar_dy = np.ma.ones(np.shape(var_mean))*valmask
delta_lat = 1 # degrees
# Roll latitude
var_mean_plus1lat = np.roll(var_mean, -1, axis=1)

# -- Determine derivative dvar_dy
dvar_dy = (var_mean_plus1lat - var_mean)/delta_lat
dvar_dy[:,-1,:] = valmask
#print(np.shape(dvar_dy))

# -- Mask
dvar_dy[np.nonzero(dvar_dy == valmask)] = np.ma.masked

# -- Separate for each basin
dvar_dy_p = dvar_dy[:,:,1]
dvar_dy_a = dvar_dy[:,:,2]
dvar_dy_i = dvar_dy[:,:,3]



# ==== Partial derivative of latitude for time dy/dt --> dy_dt and density for time dsigma/dt --> dsig_dt ===

### NEW METHOD ###

# -- Initialize dy/dt to masked value
dy_dt = np.ma.masked_all((levN,latN,basinN))
# -- Initialize dsigma/dt to masked value
dsig_dt = np.ma.masked_all((levN,latN,basinN))

# -- Define window for search
window_lat = 5 #degrees C
window_density = 0.5 #kg.m-3

# -- Define threshold for search
if varname['longN'] == 'salinity':
    epsilon = 0.01
else:
    epsilon = 0.5

# -- HR latitude
step = 0.01
lat_hr = np.arange(-70, 70 + step/2, step)

# -- Initialize var_2000_hr
var_2000_hr = np.ma.masked_all((levN,len(lat_hr),basinN))

## TESTER AVEC INTERPOLATION + FINE POUR COMPARER
# Start loop
for ibasin in range(1,4):

    # 2D interpolation
    interp  = interpolate.interp2d(lat,density,var_2000[:,:,ibasin])
    var_2000_hr[:,:,ibasin] = interp(lat_hr,density)

    # Start loops over density and latitude
    for isig in range(levN):
        for ilat in range(latN):

            if var_1950[isig,ilat,ibasin] != valmask:

                # Find latitude and density indices within the window [lat-window_lat;lat+window_lat],[sig-window_density;sig+window_density]
                ilat_window = np.flatnonzero(np.ma.abs(lat_hr[:]-lat[ilat]) < window_lat)
                isig_window = np.flatnonzero(np.ma.abs(density[:] - density[isig]) < window_density)
                if len(ilat_window) != 0 and len(isig_window) != 0:
                    lat_window = lat_hr[ilat_window]
                    sig_window = density[isig_window]
                    #print(basin[ibasin], density[isig], lat[ilat], var_1950[isig,ilat,ibasin])
                    # Find close salinity in 2000 to salinity in 1950 at isig,ilat within the window
                    window_iijj = np.ma.nonzero(np.ma.abs(var_2000_hr[isig_window[0]:isig_window[-1]+1,ilat_window[0]:ilat_window[-1]+1,ibasin]
                                                          - var_1950[isig,ilat,ibasin]) < epsilon)
                    # window_iijj gives us the indices of the closest salinities in isig_window, ilat_window
                    if len(window_iijj[0]) != 0 :
                        #print(window_iijj)
                        ii = isig_window[window_iijj[0]] # accessing the indices of the full density vector
                        jj = ilat_window[window_iijj[1]] # accessing the indices of the full hr latitude vector

                        delta_sig = density[ii] - density[isig]
                        delta_lat = lat_hr[jj] - lat[ilat]
                        dist = np.ma.sqrt(np.square(delta_sig) + np.square(delta_lat)) # 1D array of distances from the point in 1950

                        dsig_dt[isig,ilat,ibasin] = -delta_sig[np.ma.argmin(dist)]
                        dy_dt[isig,ilat,ibasin] = -delta_lat[np.ma.argmin(dist)]
                        #print(dsig_dt[isig,ilat,ibasin], dy_dt[isig,ilat,ibasin])



var_2000_hr[np.nonzero(var_2000_hr == valmask)] = np.ma.masked


# ==== Partial derivative of latitude for time dy/dt --> dy_dt ===

# # -- Determine minimum values for each sigma to obtain a threshold
# dvar_dy_min_p = np.ma.min(np.ma.abs(dvar_dy_p), axis=1)
# dvar_dy_min_a = np.ma.min(np.ma.abs(dvar_dy_a), axis=1)
# dvar_dy_min_i = np.ma.min(np.ma.abs(dvar_dy_i), axis=1)

# # -- Plot min
# plt.figure()
# plt.plot(density,dvar_dy_min_a,'r')
# plt.plot(density,dvar_dy_min_p,'b')
# plt.plot(density,dvar_dy_min_i,'g')
# plt.xlabel('density')
# plt.legend(('Atlantic','Pacific','Indian'))
# plt.title('dvar_dy min')

# # -- Define max distance (in lat degrees) between same var_1950 and var_2000
# max_lat_shift = 5
#
# # -- Initialize dy/dt to masked value
# dy_dt = np.ma.masked_all((levN,latN,basinN))
#
# delta_var_min = np.ma.masked_all((levN,latN,4))
#
# # -- Initialize variable to store and plot the interpolated data
# #var_2000_lat_hr = np.ma.masked_all((levN,len(lat_hr),basinN))
#
#
# ### OLD METHOD ###
#
# for ibasin in range(1,4):
#
#     for isig in range(levN):
#
#         if np.any(var_1950[isig,:,ibasin] != valmask):
#
#             # # -- Interpolate on a finer latitude grid
#             # z = var_2000[isig, :, ibasin].squeeze()
#             # # Find contiguous unmasked data of var_2000 at isig to interpolate
#             # slice_z = np.ma.flatnotmasked_contiguous(z)
#             # if len(slice_z) != 0:
#             #     # print(ibasin, isig, len(slice_z))
#             #     var_2000_hr = np.ma.masked_all(len(lat_hr))
#             #     for j in range(len(slice_z)):
#             #         #print(basin[ibasin], density[isig])
#             #         lat1 = lat[slice_z[j]]
#             #         z1 = z[slice_z[j]]
#             #         if len(lat1) > 2:
#             #             #print(lat1[0],lat1[-1],z1[0],z1[-1])
#             #             lat1_hr = np.arange(lat1[0], lat1[-1] + step/2, step)
#             #             #print(len(lat1_hr))
#             #             interp = interpolate.interp1d(lat1, z1, bounds_error=False)
#             #             var_2000_hr1 = interp(lat1_hr)
#             #             i1 = np.searchsorted(lat_hr, lat1_hr[0])
#             #             i2 = np.searchsorted(lat_hr, lat1_hr[-1])
#             #             # # Clean up exceptions
#             #             # if (69.999 <= lat1_hr[-1] <= 70.0001) :
#             #             #     var_2000_hr1 = var_2000_hr1[:-1]
#             #             #     print(np.shape(var_2000_hr1))
#             #             var_2000_hr[i1:i2 + 1] = var_2000_hr1
#             #             # print(lat1_hr[0], lat1_hr[-1], var_2000_hr[i1], var_2000_hr[i2], len(var_2000_hr))
#             #
#             #     # -- Save interpolated data to plot and check if it makes sense
#             #     var_2000_lat_hr[isig, :, ibasin] = var_2000_hr
#
#                 # if (ibasin == 2 and density[isig] == 27) or (ibasin==1 and density[isig]==26):
#                 #     # -- Plot S(y)
#                 #     plt.figure()
#                 #     plt.plot(lat, var_1950[isig, :, ibasin])
#                 #     plt.plot(lat_hr, var_2000_hr, ls='--')
#                 #     plt.xlabel('Latitude')
#                 #     plt.ylabel('Salinity')
#                 #     plt.legend(('1950','2000'))
#                 #     plt.title('S(y), sigma=%.2f, basin=%s' %(density[isig],basin[ibasin]))
#
#
#             for ilat in range(latN):
#
#                 if (var_1950[isig,ilat,ibasin] != valmask):
#
#                     # Look at the slope of dS/dy - if it is close to zero, only look at the vertical shift
#                     if np.ma.abs(dvar_dy[isig,ilat,ibasin]) <= 0.01 :
#                         dy_dt[isig,ilat,ibasin] =0
#                     else:
#
#                         # Find latitude indices where abs(lat_hr-lat[ilat])< 10 degrees
#                         ilat_interval = np.flatnonzero(np.ma.abs(lat[:] - lat[ilat]) < max_lat_shift)
#
#                         if len(ilat_interval) != 0:
#                             lat_interval = lat[ilat_interval]
#                             # print('Hello')
#                             # print(ibasin,density[isig],lat[ilat],var_1950[isig,ilat,ibasin])
#                             # print(lat_interval)
#                             # print(var_2000_hr[ilat_interval])
#
#                             ilat2 = ilat_interval[np.ma.argmin(np.ma.abs(var_2000[isig,ilat_interval,ibasin] - var_1950[isig,ilat,ibasin]))]
#                             lat2 = lat[ilat2]
#                             # Save value to plot
#                             delta_var_min[isig,ilat,ibasin] = np.ma.abs(var_2000[isig,ilat2,ibasin]-var_1950[isig,ilat,ibasin])
#                             # print(delta_var_min[isig,ilat,ibasin])
#                             # print(lat2, lat2 - lat[ilat])
#                             dy_dt[isig, ilat, ibasin] = -(lat2 - lat[ilat])
#                             # if ibasin ==1 and density[isig] == 26 and lat[ilat] == 35:
#                             #     print(basin[ibasin], density[isig], lat[ilat])
#                             #     print(dy_dt[isig,ilat,ibasin], var_1950[isig,ilat,ibasin], var_2000_hr[ilat2])






### OLD OLD METHOD ###

# for ibasin in range(1,4) :
#
#     if ibasin == 1 :
#         dvar_dy_min = dvar_dy_min_p
#     if ibasin == 2 :
#         dvar_dy_min = dvar_dy_min_a
#     if ibasin == 3 :
#         dvar_dy_min = dvar_dy_min_i
#
#     for isig in range(levN):
#         min = dvar_dy_min[isig]
#         # -- Test if min is masked (no data or just one value for dS/dy at isig)
#         if min != valmask:
#             if min < 0.01:
#                 min = 0.01
#             #print(ibasin,isig,min)
#             # -- Interpolate the salinity/temp on a finer latitude grid
#             z = var_2000[isig, :, ibasin].squeeze()
#             # Find contiguous unmasked data of var_2000 at isig to interpolate
#             slice_z = np.ma.flatnotmasked_contiguous(z)
#             if len(slice_z) != 0:
#                 #print(ibasin, isig, len(slice_z))
#                 var_2000_hr = np.ma.masked_all(len(lat_hr))
#                 for j in range(len(slice_z)):
#                     #print(ibasin, isig)
#                     lat1 = lat[slice_z[j]]
#                     z1 = z[slice_z[j]]
#                     if len(lat1) > 2 :
#                         #print(lat1[0],lat1[-1],z1[0],z1[-1])
#                         lat1_hr = np.arange(lat1[0], lat1[-1] + 0.005, 0.01)
#                         #print(len(lat1_hr))
#                         interp = interpolate.interp1d(lat1, z1, bounds_error=False)
#                         var_2000_hr1 = interp(lat1_hr)
#                         i1 = np.searchsorted(lat_hr, lat1_hr[0])
#                         i2 = np.searchsorted(lat_hr, lat1_hr[-1])
#                         var_2000_hr[i1:i2+1] = var_2000_hr1
#                         #print(lat1_hr[0], lat1_hr[-1], var_2000_hr[i1], var_2000_hr[i2], len(var_2000_hr))
#
#                 # -- Save interpolated data to plot and check if it makes sense
#                 var_2000_lat_hr[isig, :, ibasin] = var_2000_hr
#
#                 for ilat in range(latN):
#
#                     # -- Test if value is masked
#                     if var_1950[isig,ilat,ibasin] != valmask :
#
#                         # -- Determine latitude indices where salinity/temp in 2000 at isig is close to
#                         # salinity in 1950 at isig, ilat
#                         ind_lat2 = np.flatnonzero(np.ma.abs(var_2000_hr[:] - var_1950[isig,ilat,ibasin]) < min)
#                         if len(ind_lat2) != 0 :
#                             #print('Hello')
#                             #print(ibasin,density[isig],lat[ilat])
#                             #print("Treshold: %.4f" %(min,))
#                             #print(var_1950[isig,ilat,ibasin])
#                             #print(ind_lat2)
#                             # -- Keep closest point in ind_lat2 to ilat
#                             ilat2 = ind_lat2[np.argmin(np.abs(lat_hr[ind_lat2] - lat[ilat]))]
#                             lat2 = lat_hr[ilat2]
#                             #print(ilat2,lat2, var_2000_hr[ilat2])
#                             #print(lat[ilat],lat2)
#                             # -- Save np.min(np.abs(lat_hr[ind_lat2] - lat[ilat]) in array
#                             min_distance_lat[isig,ilat,ibasin] = lat2 - lat[ilat]
#                             #print(min_distance_lat[isig,ilat,ibasin])
#                             if np.ma.abs(min_distance_lat[isig,ilat,ibasin]) <= max_isopyc_shift:
#                                 #print('Good')
#                                 dy_dt[isig,ilat,ibasin] = (lat2 - lat[ilat])/delta_t
#                                 #print(min_distance_lat[isig,ilat,ibasin], dy_dt[isig,ilat,ibasin])
#                             # else:
#                             #     print('Bad')



# -- Separate for each basin
dy_dt_p = dy_dt[:,:,1]
dy_dt_a = dy_dt[:,:,2]
dy_dt_i = dy_dt[:,:,3]


##del var_2000_hr


# ==== Partial derivative of sigma for time dsigma/dt ===

# # -- Determine minimum values for each latitude to obtain a threshold
# dvar_dsig_min_p = np.ma.min(np.ma.abs(dvar_dsig_p), axis=0)
# dvar_dsig_min_a = np.ma.min(np.ma.abs(dvar_dsig_a), axis=0)
# dvar_dsig_min_i = np.ma.min(np.ma.abs(dvar_dsig_i), axis=0)

# # -- Plot min
# plt.figure()
# plt.plot(lat,dvar_dsig_min_a,'r')
# plt.plot(lat,dvar_dsig_min_p,'b')
# plt.plot(lat,dvar_dsig_min_i,'g')
# plt.xlabel('latitude')
# plt.legend(('Atlantic','Pacific','Indian'))
# plt.title('dvar_dsig min')

# # -- Define max distance (in density units) between same var_1950 and var_2000
# max_density_shift1 = 0.3
# max_density_shift2 = 0.3
#
# # -- Initialize dsigma/dt to masked value
# dsig_dt = np.ma.masked_all((levN,latN,basinN))
#
# delta_var_min = np.ma.masked_all((levN,latN,basinN))

# -- Initialize variable to store and plot the interpolated data
#var_2000_sig_hr = np.ma.masked_all((len(density_hr), latN, basinN))

### OLD METHOD ###

# for ibasin in range(1,4):
#
#     for ilat in range(latN):
#
#         if np.any(var_1950[:,ilat,ibasin] != valmask):
#
#             # # -- Interpolate on a finer density grid
#             # z = var_2000[:, ilat, ibasin].squeeze()
#             # # Find contiguous unmasked data of var_2000 at ilat to interpolate
#             # slice_z = np.ma.flatnotmasked_contiguous(z)
#             # if len(slice_z) != 0:
#             #     # print(ibasin, ilat, len(slice_z))
#             #     var_2000_hr = np.ma.masked_all(len(density_hr))
#             #     density1 = density[slice_z[0]]
#             #     z1 = z[slice_z[0]]
#             #     if len(density1) > 2:
#             #         # print(density1[0],density1[-1],z1[0],z1[-1])
#             #         i1 = np.searchsorted(density_hr, density1[0])
#             #         i2 = np.searchsorted(density_hr, density1[-1])
#             #         density1_hr = density_hr[i1:i2 + 1]
#             #         # print(len(density1_hr))
#             #         interp = interpolate.interp1d(density1, z1, bounds_error=False)
#             #         var_2000_hr1 = interp(density1_hr)
#             #         var_2000_hr[i1:i2 + 1] = var_2000_hr1
#             #         # print(density1_hr[0], density1_hr[-1], var_2000_hr[i1], var_2000_hr[i2], len(var_2000_hr))
#             #
#             #     # -- Save interpolated data to plot and check if it makes sense
#             #     var_2000_sig_hr[:, ilat, ibasin] = var_2000_hr
#
#                 # if ibasin == 2 and lat[ilat] == 0:
#                     # # -- Plot S(sigma)
#                     # plt.figure()
#                     # plt.plot(var_1950[:, ilat, ibasin], density)
#                     # plt.plot(var_2000_hr, density, ls='--')
#                     # plt.gca().invert_yaxis()
#                     # plt.xlabel('Salinity')
#                     # plt.ylabel('Density')
#                     # plt.legend(('1950','2000'))
#                     # plt.title('S(sigma)')
#
#             for isig in range(levN):
#
#                 if density[isig] < 26.5 :
#                     max_density_shift = max_density_shift1
#                 else :
#                     max_density_shift = max_density_shift2
#
#                 if (var_1950[isig,ilat,ibasin] != valmask) :
#
#                     # Look at the slope of dS/dsigma - if it is close to zero, only look at the horizontal shift
#                     if np.ma.abs(dvar_dsig[isig, ilat, ibasin]) <= 0.25:
#                         dsig_dt[isig, ilat, ibasin] = 0
#                     else:
#
#                         # Find density indices where abs(density_hr-density[isig])<0.5 kg.m-3
#                         isig_interval = np.flatnonzero(np.ma.abs(density[:] - density[isig]) < max_density_shift)
#
#                         if len(isig_interval) != 0:
#                             sig_interval = density[isig_interval]
#                             # print(ibasin,density[isig],lat[ilat],var_1950[isig,ilat,ibasin])
#                             # print(sig_interval)
#                             # print(var_2000_hr[isig_interval])
#
#                             isig2 = isig_interval[np.ma.argmin(np.ma.abs(var_2000[isig_interval,ilat,ibasin] - var_1950[isig,ilat,ibasin]))]
#                             sig2 = density[isig2]
#                             # Save value to plot
#                             delta_var_min[isig,ilat,ibasin] = np.ma.abs(var_2000[isig2,ilat,ibasin]-var_1950[isig,ilat,ibasin])
#                             # print(delta_var_min[isig,ilat,ibasin])
#                             dsig_dt[isig, ilat, ibasin] = -(sig2 - density[isig])
#
#                                 # if ibasin == 1 and lat[ilat] == 0 and 26.5 <= density[isig] <= 27 :
#                                 # if ibasin == 2 and density[isig] == 27.55 and lat[ilat] == -20:
#                                 #     print('Hello')
#                                 #     print(ibasin, lat[ilat], density[isig], var_1950[isig,ilat,ibasin])
#                                 #     print(sig2, sig2-density[isig], dsig_dt[isig,ilat,ibasin],delta_var_min[isig,ilat,ibasin])
#                                 #
#                                 # if ibasin == 1 and density[isig] == 25 and lat[ilat] == 0 :
#                                 #     print(ibasin, lat[ilat], density[isig], var_1950[isig,ilat,ibasin])
#                                 #     print(sig2, sig2 - density[isig], dsig_dt[isig, ilat, ibasin],
#                                 #       delta_var_min[isig, ilat, ibasin])





# ### OLD OLD METHOD ###
#
# for ibasin in range(1,4):
#
#     if ibasin == 1:
#         dvar_dsig_min = dvar_dsig_min_p
#     if ibasin == 2:
#         dvar_dsig_min = dvar_dsig_min_a
#     if ibasin == 3:
#         dvar_dsig_min = dvar_dsig_min_i
#
#
#     for ilat in range(latN):
#         min = dvar_dsig_min[ilat]
#         #print(ilat,min)
#         # Test if min is masked (no data or just one value for dS/dsigma at ilat)
#         if min != valmask:
#             # if min < 0.05:
#             #     min = 0.05
#             #print(ibasin, lat[ilat], min)
#             # -- Interpolate on a finer density grid
#             z = var_2000[:, ilat, ibasin].squeeze()
#             # Find contiguous unmasked data of var_2000 at ilat to interpolate
#             slice_z = np.ma.flatnotmasked_contiguous(z)
#             if len(slice_z) != 0:
#                 #print(ibasin, ilat, len(slice_z))
#                 var_2000_hr = np.ma.masked_all(len(density_hr))
#                 density1 = density[slice_z[0]]
#                 z1 = z[slice_z[0]]
#                 if len(density1) > 2 :
#                     #print(density1[0],density1[-1],z1[0],z1[-1])
#                     i1 = np.searchsorted(density_hr, density1[0])
#                     i2 = np.searchsorted(density_hr, density1[-1])
#                     density1_hr = density_hr[i1:i2+1]
#                     #print(len(density1_hr))
#                     interp = interpolate.interp1d(density1, z1, bounds_error=False)
#                     var_2000_hr1 = interp(density1_hr)
#                     var_2000_hr[i1:i2+1] = var_2000_hr1
#                     #print(density1_hr[0], density1_hr[-1], var_2000_hr[i1], var_2000_hr[i2], len(var_2000_hr))
#
#                 # -- Save interpolated data to plot and check if it makes sense
#                 var_2000_sig_hr[:, ilat, ibasin] = var_2000_hr

                # for isig in range(levN) :
                #
                #     # Test if value is masked
                #     if var_1950[isig,ilat,ibasin] != valmask :
                #         # Determine density indices where salinity/temp in 2000 at ilat is close to
                #         # salinity in 1950 at ilat, isig
                #
                #         var_min = np.min(np.abs(var_2000_hr[:]-var_1950[isig,ilat,ibasin]))
                #         if var_min != valmask :
                #             print(var_min)
                #             ind_sig2 = np.flatnonzero(np.ma.abs(var_2000_hr[:] - var_1950[isig,ilat,ibasin]) < 1.5*var_min)
                #             if len(ind_sig2) != 0 :
                #                 # print('Hello')
                #                 # print(ibasin,density[isig],lat[ilat])
                #                 # print("Treshold: %.4f" %(min,))
                #                 # print(var_1950[isig,ilat,ibasin])
                #                 # print(ind_sig2)
                #                 # Keep closest point in ind_sig2 to isig
                #                 isig2 = ind_sig2[np.argmin(np.abs(density_hr[ind_sig2] - density[isig]))]
                #                 sig2 = density_hr[isig2]
                #                 # print(sig2, var_2000_hr[isig2])
                #                 # Save ind_sig2[np.argmin(np.abs(ind_sig2 - isig))] in array
                #                 min_distance_sig[isig, ilat, ibasin] = sig2 - density[isig]
                #                 # print(min_distance_sig[isig, ilat, ibasin])
                #                 if np.ma.abs(min_distance_sig[isig,ilat,ibasin]) <= max_isopyc_shift:
                #                     # print('Good')
                #                     dsig_dt[isig,ilat,ibasin] = (sig2 - density[isig])/delta_t
                #                     #print(min_distance_sig[isig,ilat,ibasin], dsig_dt[isig,ilat,ibasin])
                #                 # else:
                #                 #     print('Bad')


# -- Separate for each basin
dsig_dt_p = dsig_dt[:,:,1]
dsig_dt_a = dsig_dt[:,:,2]
dsig_dt_i = dsig_dt[:,:,3]



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
varPac = {'name': 'Pacific', 'var_change': var_change_p, 'var_mean': var_mean_p, 'var_error': var_change_er_p,
          'var_1950': var_1950[:,:,1],  'var_2000': var_2000[:,:,1],
          'var_change_res': var_change_res_p, 'dvar_dsig': dvar_dsig_p, 'dvar_dy': dvar_dy_p,
          'dy_dt': dy_dt_p, 'dsig_dt': dsig_dt_p, 'bowl': bowl_p,
          #'var_2000_lat_hr': var_2000_lat_hr[:,:,1], 'var_2000_sig_hr': var_2000_sig_hr[:,:,1],
          #'delta_var_min': delta_var_min[:,:,1],
           'var_2000_hr': var_2000_hr[:,:,1]}
varAtl = {'name': 'Atlantic', 'var_change': var_change_a, 'var_mean': var_mean_a, 'var_error': var_change_er_a,
            'var_1950': var_1950[:,:,2],  'var_2000': var_2000[:,:,2],
          'var_change_res': var_change_res_a, 'dvar_dsig': dvar_dsig_a, 'dvar_dy': dvar_dy_a,
          'dy_dt': dy_dt_a, 'dsig_dt': dsig_dt_a, 'bowl': bowl_a,
          #'var_2000_lat_hr': var_2000_lat_hr[:, :, 2], 'var_2000_sig_hr': var_2000_sig_hr[:, :, 2],
          #'delta_var_min': delta_var_min[:,:, 2]
          'var_2000_hr': var_2000_hr[:,:,2]}
varInd = {'name': 'Indian', 'var_change': var_change_i, 'var_mean': var_mean_i, 'var_error': var_change_er_i,
            'var_1950': var_1950[:,:,3],  'var_2000': var_2000[:,:,3],
          'var_change_res': var_change_res_i, 'dvar_dsig': dvar_dsig_i, 'dvar_dy': dvar_dy_i,
          'dy_dt': dy_dt_i, 'dsig_dt': dsig_dt_i, 'bowl': bowl_i,
          #'var_2000_lat_hr': var_2000_lat_hr[:, :, 3], 'var_2000_sig_hr': var_2000_sig_hr[:, :, 3],
          #'delta_var_min': delta_var_min[:,:, 3]
          'var_2000_hr': var_2000_hr[:,:,3]}


# # ==== Total change ====
# fig1, axes = plt.subplots(nrows=2, ncols=3, figsize=(17,5))
#
# cnplot1 = zonal_2D(plt, 'total', axes[0,0], axes[1,0], 'left', lat, density, varAtl, minmax, domrho, clevsm, clevsm_bold)
#
# cnplot1 = zonal_2D(plt, 'total', axes[0,1], axes[1,1], 'mid', lat, density, varPac, minmax, domrho, clevsm, clevsm_bold)
#
# cnplot1 = zonal_2D(plt, 'total', axes[0,2], axes[1,2], 'right', lat, density, varInd, minmax, domrho, clevsm, clevsm_bold)
#
#
# plt.subplots_adjust(hspace=.0001, wspace=0.05, left=0.04, right=0.86)
#
# cb = plt.colorbar(cnplot1[0], ax=axes.ravel().tolist(), ticks=cnplot1[1][::3])
# cb.set_label('%s (%s)' % (legVar, unit), fontweight='bold')
#
# plt.suptitle('Total %s changes (%s)' %(legVar, name),
#           fontweight='bold', fontsize=14, verticalalignment='top')
#
#
# ==== Isopycnal migration term ====
fig2, axes = plt.subplots(nrows=2, ncols=3, figsize=(17,5))

cnplot2 = zonal_2D(plt, 'isopyc_mig', axes[0,0], axes[1,0], 'left', lat, density, varAtl, minmax, domrho)

cnplot2 = zonal_2D(plt, 'isopyc_mig', axes[0,1], axes[1,1], 'mid', lat, density, varPac, minmax, domrho)

cnplot2 = zonal_2D(plt, 'isopyc_mig', axes[0,2], axes[1,2], 'right', lat, density, varInd, minmax, domrho)


plt.subplots_adjust(hspace=.0001, wspace=0.05, left=0.04, right=0.86)

cb = plt.colorbar(cnplot2[0], ax=axes.ravel().tolist(), ticks=cnplot2[1][::3])
cb.set_label('%s (%s)' % (legVar, unit), fontweight='bold')

plt.suptitle('%s changes due to isopycnal migration (%s)' %(legVar, name),
          fontweight='bold', fontsize=14, verticalalignment='top')

plt.savefig('/home/ysilvy/Density_bining/Yona_analysis/figures/obs/zonal_ys/isopycnal_migration/'
            +'Durack_isopyc_migration_term.png', bbox_inches='tight')


# ==== Residual term ====
fig3, axes = plt.subplots(nrows=2, ncols=3, figsize=(17,5))

cnplot3 = zonal_2D(plt, 'residual', axes[0,0], axes[1,0], 'left', lat, density, varAtl, minmax, domrho)

cnplot3 = zonal_2D(plt, 'residual', axes[0,1], axes[1,1], 'mid', lat, density, varPac, minmax, domrho)

cnplot3 = zonal_2D(plt, 'residual', axes[0,2], axes[1,2], 'right', lat, density, varInd, minmax, domrho)


plt.subplots_adjust(hspace=.0001, wspace=0.05, left=0.04, right=0.86)

cb = plt.colorbar(cnplot3[0], ax=axes.ravel().tolist(), ticks=cnplot3[1][::3])
cb.set_label('%s (%s)' % (legVar, unit), fontweight='bold')

plt.suptitle('Residual %s changes (%s)' %(legVar, name),
          fontweight='bold', fontsize=14, verticalalignment='top')

plt.savefig('/home/ysilvy/Density_bining/Yona_analysis/figures/obs/zonal_ys/isopycnal_migration/'
            +'Durack_residual_term.png', bbox_inches='tight')


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

plt.savefig('/home/ysilvy/Density_bining/Yona_analysis/figures/obs/zonal_ys/isopycnal_migration/'
            +'Durack_lat_term.png', bbox_inches='tight')


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

plt.savefig('/home/ysilvy/Density_bining/Yona_analysis/figures/obs/zonal_ys/isopycnal_migration/'
            +'Durack_density_term.png', bbox_inches='tight')


# ==== dS/dsigma ====

fig6, axes = plt.subplots(nrows=2, ncols=3, figsize=(17,5))

cnplot6 = zonal_2D(plt, 'dvar_dsig', axes[0,0], axes[1,0], 'left', lat, density, varAtl, minmax, domrho)

cnplot6 = zonal_2D(plt, 'dvar_dsig', axes[0,1], axes[1,1], 'mid', lat, density, varPac, minmax, domrho)

cnplot6 = zonal_2D(plt, 'dvar_dsig', axes[0,2], axes[1,2], 'right', lat, density, varInd, minmax, domrho)


plt.subplots_adjust(hspace=.0001, wspace=0.05, left=0.04, right=0.86)

cb = plt.colorbar(cnplot6[0], ax=axes.ravel().tolist())

plt.suptitle('dS/dsigma (%s)' %(name,),
          fontweight='bold', fontsize=14, verticalalignment='top')

plt.savefig('/home/ysilvy/Density_bining/Yona_analysis/figures/obs/zonal_ys/isopycnal_migration/'
            +'Durack_dS_dsigma.png', bbox_inches='tight')


# # ==== dsigma/dt ====

fig7, axes = plt.subplots(nrows=2, ncols=3, figsize=(17,5))

cnplot7 = zonal_2D(plt, 'dsig_dt', axes[0,0], axes[1,0], 'left', lat, density, varAtl, minmax, domrho)

cnplot7 = zonal_2D(plt, 'dsig_dt', axes[0,1], axes[1,1], 'mid', lat, density, varPac, minmax, domrho)

cnplot7 = zonal_2D(plt, 'dsig_dt', axes[0,2], axes[1,2], 'right', lat, density, varInd, minmax, domrho)


plt.subplots_adjust(hspace=.0001, wspace=0.05, left=0.04, right=0.86)

cb = plt.colorbar(cnplot7[0], ax=axes.ravel().tolist(), ticks=cnplot7[1][::3])

plt.suptitle('dsigma/dt (%s)' %(name,),
          fontweight='bold', fontsize=14, verticalalignment='top')

plt.savefig('/home/ysilvy/Density_bining/Yona_analysis/figures/obs/zonal_ys/isopycnal_migration/'
            +'Durack_dsig_dt.png', bbox_inches='tight')

# ==== dy/dt ====

fig8, axes = plt.subplots(nrows=2, ncols=3, figsize=(17,5))

cnplot8 = zonal_2D(plt, 'dy_dt', axes[0,0], axes[1,0], 'left', lat, density, varAtl, minmax, domrho)

cnplot8 = zonal_2D(plt, 'dy_dt', axes[0,1], axes[1,1], 'mid', lat, density, varPac, minmax, domrho)

cnplot8 = zonal_2D(plt, 'dy_dt', axes[0,2], axes[1,2], 'right', lat, density, varInd, minmax, domrho)


plt.subplots_adjust(hspace=.0001, wspace=0.05, left=0.04, right=0.86)

cb = plt.colorbar(cnplot8[0], ax=axes.ravel().tolist(), ticks=cnplot8[1][::3])

plt.suptitle('dy/dt (%s)' %(name,),
          fontweight='bold', fontsize=14, verticalalignment='top')

plt.savefig('/home/ysilvy/Density_bining/Yona_analysis/figures/obs/zonal_ys/isopycnal_migration/'
            +'Durack_dy_dt.png', bbox_inches='tight')


#
# ==== dS/dy ====

fig9, axes = plt.subplots(nrows=2, ncols=3, figsize=(17,5))

cnplot9 = zonal_2D(plt, 'dvar_dy', axes[0,0], axes[1,0], 'left', lat, density, varAtl, minmax, domrho)

cnplot9 = zonal_2D(plt, 'dvar_dy', axes[0,1], axes[1,1], 'mid', lat, density, varPac, minmax, domrho)

cnplot9 = zonal_2D(plt, 'dvar_dy', axes[0,2], axes[1,2], 'right', lat, density, varInd, minmax, domrho)


plt.subplots_adjust(hspace=.0001, wspace=0.05, left=0.04, right=0.86)

cb = plt.colorbar(cnplot9[0], ax=axes.ravel().tolist())

plt.suptitle('dS/dy (%s)' %(name,),
          fontweight='bold', fontsize=14, verticalalignment='top')

plt.savefig('/home/ysilvy/Density_bining/Yona_analysis/figures/obs/zonal_ys/isopycnal_migration/'
            +'Durack_dS_dy.png', bbox_inches='tight')


# ==== var_2000_sig_hr ====

# fig10, axes = plt.subplots(nrows=2, ncols=3, figsize=(17,5))
#
# cnplot10 = zonal_2D(plt, 'var_2000_sig_hr', axes[0,0], axes[1,0], 'left', lat, density_hr, varAtl, minmax, domrho)
#
# cnplot10 = zonal_2D(plt, 'var_2000_sig_hr', axes[0,1], axes[1,1], 'mid', lat, density_hr, varPac, minmax, domrho)
#
# cnplot10 = zonal_2D(plt, 'var_2000_sig_hr', axes[0,2], axes[1,2], 'right', lat, density_hr, varInd, minmax, domrho)
#
#
# plt.subplots_adjust(hspace=.0001, wspace=0.05, left=0.04, right=0.86)
#
# cb = plt.colorbar(cnplot10[0], ax=axes.ravel().tolist())
#
# plt.suptitle('2000 mean field interpolated on a finer density grid (%s)' %(name,),
#           fontweight='bold', fontsize=14, verticalalignment='top')
#
# ==== var_2000_hr ====

# fig11, axes = plt.subplots(nrows=2, ncols=3, figsize=(17,5))
#
# cnplot11 = zonal_2D(plt, 'var_2000_hr', axes[0,0], axes[1,0], 'left', lat_hr, density, varAtl, minmax, domrho)
#
# cnplot11 = zonal_2D(plt, 'var_2000_hr', axes[0,1], axes[1,1], 'mid', lat_hr, density, varPac, minmax, domrho)
#
# cnplot11 = zonal_2D(plt, 'var_2000_hr', axes[0,2], axes[1,2], 'right', lat_hr, density, varInd, minmax, domrho)
#
#
# plt.subplots_adjust(hspace=.0001, wspace=0.05, left=0.04, right=0.86)
#
# cb = plt.colorbar(cnplot11[0], ax=axes.ravel().tolist(), ticks = cnplot11[1])
#
# plt.suptitle('2000 mean field interpolated on a finer lat grid (%s)' %(name,),
#           fontweight='bold', fontsize=14, verticalalignment='top')


# ==== min_distance_lat ====

# fig11, axes = plt.subplots(nrows=2, ncols=3, figsize=(17,5))
#
# cnplot11 = zonal_2D(plt, 'min_dist_lat', axes[0,0], axes[1,0], 'left', lat, density, varAtl, minmax, domrho)
#
# cnplot11 = zonal_2D(plt, 'min_dist_lat', axes[0,1], axes[1,1], 'mid', lat, density, varPac, minmax, domrho)
#
# cnplot11 = zonal_2D(plt, 'min_dist_lat', axes[0,2], axes[1,2], 'right', lat, density, varInd, minmax, domrho)
#
#
# plt.subplots_adjust(hspace=.0001, wspace=0.05, left=0.04, right=0.86)
#
# cb = plt.colorbar(cnplot11[0], ax=axes.ravel().tolist())
#
# plt.suptitle('Min distance latitude',
#           fontweight='bold', fontsize=14, verticalalignment='top')

# ==== min_distance_sig ====

# fig12, axes = plt.subplots(nrows=2, ncols=3, figsize=(17,5))
#
# cnplot12 = zonal_2D(plt, 'min_dist_sig', axes[0,0], axes[1,0], 'left', lat, density, varAtl, minmax, domrho)
#
# cnplot12 = zonal_2D(plt, 'min_dist_sig', axes[0,1], axes[1,1], 'mid', lat, density, varPac, minmax, domrho)
#
# cnplot12 = zonal_2D(plt, 'min_dist_sig', axes[0,2], axes[1,2], 'right', lat, density, varInd, minmax, domrho)
#
#
# plt.subplots_adjust(hspace=.0001, wspace=0.05, left=0.04, right=0.86)
#
# cb = plt.colorbar(cnplot12[0], ax=axes.ravel().tolist())
#
# plt.suptitle('Min distance density',
#           fontweight='bold', fontsize=14, verticalalignment='top')


# ==== delta_var_min ====

# fig13, axes = plt.subplots(nrows=2, ncols=3, figsize=(17,5))
#
# cnplot13 = zonal_2D(plt, 'delta_var_min', axes[0,0], axes[1,0], 'left', lat, density, varAtl, minmax, domrho)
#
# cnplot13 = zonal_2D(plt, 'delta_var_min', axes[0,1], axes[1,1], 'mid', lat, density, varPac, minmax, domrho)
#
# cnplot13 = zonal_2D(plt, 'delta_var_min', axes[0,2], axes[1,2], 'right', lat, density, varInd, minmax, domrho)
#
#
# plt.subplots_adjust(hspace=.0001, wspace=0.05, left=0.04, right=0.86)
#
# cb = plt.colorbar(cnplot13[0], ax=axes.ravel().tolist(), ticks=cnplot13[1])
#
# plt.suptitle('Minimum abs(var_2000 - var1950)',
#           fontweight='bold', fontsize=14, verticalalignment='top')



# ==== Mean contours 1950/2000 ====

# fig14, axes = plt.subplots(nrows=2, ncols=3, figsize=(17,5))
#
# zonal_2D(plt, 'mean_fields', axes[0,0], axes[1,0], 'left', lat, density, varAtl, minmax, domrho, clevsm, clevsm_bold)
#
# zonal_2D(plt, 'mean_fields', axes[0,1], axes[1,1], 'mid', lat, density, varPac, minmax, domrho, clevsm, clevsm_bold)
#
# zonal_2D(plt, 'mean_fields', axes[0,2], axes[1,2], 'right', lat, density, varInd, minmax, domrho, clevsm, clevsm_bold)
#
#
# plt.subplots_adjust(hspace=.0001, wspace=0.05, left=0.04, right=0.86)
#
# plt.suptitle('Mean %s field in 1950 and in 2000' %(legVar,),
#           fontweight='bold', fontsize=14, verticalalignment='top')


#plt.show()
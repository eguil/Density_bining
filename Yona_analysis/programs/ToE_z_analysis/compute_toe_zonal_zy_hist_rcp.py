#!/bin/env python
# -*- coding: utf-8 -*-

"""
Compute ToE hist + RCP85 vs. histNat (or PiControl) in lat/depth domain for one model
Save ToE in output files
"""
import sys
sys.path.append("/home/ysilvy/Density_bining/")
sys.path.append("/home/ysilvy/Density_bining/Yona_analysis/programs/")
import numpy as np
from netCDF4 import Dataset as open_ncfile
from maps_matplot_lib import defVarmme
#from modelsDef import defModels 
from libToE import findToE
import glob
from functions_z_analysis import maptogamma
from binDensity import rhonGrid, maskVal

# ----- Workspace ------

indir_histrcp85 = '/data/ysilvy/CMIP5_annual/'
indir_histNat = '/data/ysilvy/CMIP5_annual/'

# ----- Work ------

#models = defModels() # Error: opening the one in Density_bining/ instead of Density_bining/Yona_analysis/programs/

varname = defVarmme('salinity'); v = 'S'
#varname = defVarmme('temperature'); v = 'T'

multStd = 2. # detect ToE at multStd std dev of histNat or PiControl

use_piC = False # Signal = (hist-histNat) + RCP8.5-average(histNat), noise = std(histNat)
# use_piC = True # Signal = hist + RCP8.5 - PiControl, noise = std(PiControl)

iniyear = 1860
finalyear = 2100
deltay = 10.

# Define variable properties
legVar = varname['legVar']
unit = varname['unit']

if v=='S':
    var = 'salinity'
elif v=='T':
    var = 'temperature'

model = {'name':'IPSL-CM5A-LR'  ,'props':[6,3,11,156], 'picontrol':[1000],'correctFile':[0,0,0],
          'file_end_hist':'v20111119', 'file_end_histNat':'v20120430',
          'hist-rcp85':['r2i1p1','r3i1p1','r4i1p1']}

valmask = 1.e20

# ----- Compute zonal ToE ------

print('Computing ToE for ',model['name'])

# Index of common time interval
tstart = model['props'][2]
tend = model['props'][3]

# Read histNat
filehn='so_thetao_gamma_Oan_IPSL-CM5A-LR_historicalNat_r2i1p1_185001-201012.nc'
fhn = open_ncfile(indir_histNat+filehn,'r')

# Read historical + rcp8.5
filehrcp = 'so_thetao_gamma_Oan_IPSL-CM5A-LR_historical-rcp85_r2i1p1_185001-210012.nc'
fhrcp = open_ncfile(indir_histrcp85+filehrcp,'r')

lat = fhn.variables['latitude'][:]
depthi = fhn.variables['lev'][:]
latN = len(lat)
depthN = len(depthi)
basinN=4
timN = 240

# Read var histNat
varhn_a = fhn.variables[var][tstart:tend,1,:,:].squeeze()
varhn_p = fhn.variables[var][tstart:tend,2,:,:].squeeze()
varhn_i = fhn.variables[var][tstart:tend,3,:,:].squeeze()

# Compute std of histNat
stdvarhn_a = np.ma.std(varhn_a, axis=0)
stdvarhn_p = np.ma.std(varhn_p, axis=0)
stdvarhn_i = np.ma.std(varhn_i, axis=0)

# Compute time average of the whole histNat series (signal over projection = RCP - mean(histNat))
meanvarhn_a = np.ma.mean(varhn_a, axis=0)
meanvarhn_p = np.ma.mean(varhn_p, axis=0)
meanvarhn_i = np.ma.mean(varhn_i, axis=0)

# Reorganise i,j dims in single dimension data (speeds up loops)
varnoise_a = np.reshape(stdvarhn_a, (depthN*latN))
varnoise_p = np.reshape(stdvarhn_p, (depthN*latN))
varnoise_i = np.reshape(stdvarhn_i, (depthN*latN))

# Read var hist + RCP8.5
varhrcp_a = fhrcp.variables[var][tstart-1:tend+95,1,:,:].squeeze()
varhrcp_p = fhrcp.variables[var][tstart-1:tend+95,2,:,:].squeeze()
varhrcp_i = fhrcp.variables[var][tstart-1:tend+95,3,:,:].squeeze()

# Initialize and fill var_signal for each basin (timN, depth, latitude)
varsignal_a = np.ma.masked_all((timN,depthN,latN))
varsignal_p = np.ma.masked_all((timN,depthN,latN))
varsignal_i = np.ma.masked_all((timN,depthN,latN))

varsignal_a[0:145,:,:] = varhrcp_a[0:145,:,:]-varhn_a
varsignal_p[0:145,:,:] = varhrcp_p[0:145,:,:]-varhn_p
varsignal_i[0:145,:,:] = varhrcp_i[0:145,:,:]-varhn_i
varsignal_a[145:,:,:] = varhrcp_a[145:,:,:]-meanvarhn_a
varsignal_p[145:,:,:] = varhrcp_p[145:,:,:]-meanvarhn_p
varsignal_i[145:,:,:] = varhrcp_i[145:,:,:]-meanvarhn_i

# Initialize variable to save signal
varsignal_end = np.ma.masked_all((basinN,depthN,latN))
# Save signal
varsignal_end[1,:,:] = np.ma.average(varsignal_a[-5:,:,:],axis=0)
varsignal_end[2,:,:] = np.ma.average(varsignal_p[-5:,:,:],axis=0)
varsignal_end[3,:,:] = np.ma.average(varsignal_i[-5:,:,:],axis=0)

# Reorganise i,j dims in single dimension data (speeds up loops)
varsignal_a = np.reshape(varsignal_a, (timN, depthN*latN))
varsignal_p = np.reshape(varsignal_p, (timN, depthN*latN))
varsignal_i = np.reshape(varsignal_i, (timN, depthN*latN))
    
# Initialize ToE for each basin (depth, lat)
toe1_a = np.ma.masked_all((depthN,latN))
toe1_p = np.ma.masked_all((depthN,latN))
toe1_i = np.ma.masked_all((depthN,latN))
toe2_a = np.ma.masked_all((depthN,latN))
toe2_p = np.ma.masked_all((depthN,latN))
toe2_i = np.ma.masked_all((depthN,latN))
# Initialize output variable
varToE1 = np.ma.masked_all((basinN,depthN,latN)) # (>1std) (basin,depth,latitude)
varToE2 = np.ma.masked_all((basinN,depthN,latN)) # (>2std)

# Compute ToE as last date when diff hist+RCP - histNat is larger than mult * stddev
toe2_a = np.ma.reshape(findToE(varsignal_a, varnoise_a, multStd),(depthN,latN))
toe2_p = np.ma.reshape(findToE(varsignal_p, varnoise_p, multStd),(depthN,latN))
toe2_i = np.ma.reshape(findToE(varsignal_i, varnoise_i, multStd),(depthN,latN))
toe1_a = np.ma.reshape(findToE(varsignal_a, varnoise_a, 1),(depthN,latN))
toe1_p = np.ma.reshape(findToE(varsignal_p, varnoise_p, 1),(depthN,latN))
toe1_i = np.ma.reshape(findToE(varsignal_i, varnoise_i, 1),(depthN,latN))

# Save in output variable
varToE1[1,:,:] = toe1_a
varToE1[2,:,:] = toe1_p
varToE1[3,:,:] = toe1_i
varToE2[1,:,:] = toe2_a
varToE2[2,:,:] = toe2_p
varToE2[3,:,:] = toe2_i

varToE1.fill_value = valmask
varToE2.fill_value = valmask 

# Mask points because when calculating ToE, masked points (e.g. bathy) are set to 240 (=no emergence)
idx=np.argwhere(varsignal_end.mask==True)
varToE2[idx[:,0],idx[:,1],idx[:,2]] = np.ma.masked
varToE1[idx[:,0],idx[:,1],idx[:,2]] = np.ma.masked


# ---- Map ToE to gamma ----

# Density grid
targetrho, s_sax, del_s, N_s = rhonGrid(19, 26, 28.501, 0.2, 0.1)
levN = len(targetrho)

# Define Gamma/depth relationship for mapping
gammaz = np.ma.average(fhrcp.variables['density'][tstart:tend+95,:,:,:],axis=0)

# Map to gamma
varToE1_gamma = maptogamma(varToE1,gammaz,targetrho)
varToE2_gamma = maptogamma(varToE2,gammaz,targetrho)

# ---- Save in output file ----

fileName = 'cmip5.'+model['name']+'.'+'r2i1p1.'+legVar+'_toe_zonal_rcp_histNat.nc' 

dir = '/home/ysilvy/Density_bining/Yona_analysis/data/toe_zonal_z/toe_rcp85_histNat/'
description = 'Time of Emergence hist+rcp8.5 vs. histNat. \n' \
              'The signal in each point is hist-histNat over the historical period and ' \
              'hist-timeaverage(histNat) over the projection period. ' \
              'The noise in each point is the max standard deviation of the historicalNat run \n' \
              'The ToE is computed by using once or twice the noise as the threshold. \n' \
              'ToE Computation is done on the original depth grid, then mapped to neutral density using a gamma to pressure relationship from the model.' 

fout = open_ncfile(dir+fileName,'w', format='NETCDF4')
fout.description = description

# dimensions
fout.createDimension('basin', 4)
fout.createDimension('depth', depthN)
fout.createDimension('density', levN)
fout.createDimension('latitude', latN)

# variables
basin = fout.createVariable('basin', 'i4', ('basin',))
depth = fout.createVariable('depth', 'f4', ('depth',))
density = fout.createVariable('density', 'f4', ('density',))
latitude = fout.createVariable('latitude', 'f4', ('latitude',))
ToE1_z = fout.createVariable(var+'ToE1_z', 'f4', ('basin','depth','latitude',))
ToE2_z = fout.createVariable(var+'ToE2_z', 'f4', ('basin','depth','latitude',))
ToE1_gamma = fout.createVariable(var+'ToE1_gamma', 'f4', ('basin','density','latitude',))
ToE2_gamma = fout.createVariable(var+'ToE2_gamma', 'f4', ('basin','density','latitude',))
varchange = fout.createVariable(var+'_change', 'f4', ('basin','depth','latitude',))

# data
basin[:] =  np.arange(0,basinN)
depth[:] = depthi
density[:] = targetrho
latitude[:] = lat
ToE1_z[:,:,:] = varToE1
ToE2_z[:,:,:] = varToE2
ToE1_gamma[:,:,:] = varToE1_gamma
ToE2_gamma[:,:,:] = varToE2_gamma
varchange[:,:,:] = varsignal_end

# units
basin.units = 'basin index'
depth.units = 'm'
density.units = 'kg.m-3'
latitude.units = ''
ToE2_z.units = 'Year'
ToE2_z.long_name = 'ToE (>2std) for ' + legVar + ' on original depth grid'
ToE1_z.units = 'Year'
ToE1_z.long_name = 'ToE (>1std) for ' + legVar + ' on original depth grid'
ToE2_gamma.units = 'Year'
ToE2_gamma.long_name = 'ToE (>2std) for ' + legVar + ' mapped to gamma'
ToE1_gamma.units = 'Year'
ToE1_gamma.long_name = 'ToE (>1std) for ' + legVar + ' mapped to gamma'
varchange.units = unit
varchange.long_name = legVar + ' signal averaged in the last five years'

fout.close()


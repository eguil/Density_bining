#!/bin/env python
# -*- coding: utf-8 -*-

"""
Compute ToE 1%CO2 vs. PiControl in lat/rho domain for all available models
Save ToE in output files
"""

import numpy as np
from netCDF4 import Dataset as open_ncfile
from maps_matplot_lib import defVarmme
from modelsDef import defModelsCO2piC
from libToE import findToE
import glob

# ----- Workspace ------

indir_CO2 = '/data/ericglod/Density_binning/Prod_density_april15/mme_1pctCO2/'
indir_piC = '/data/ericglod/Density_binning/Prod_density_april15/mme_piControl/'

models = defModelsCO2piC()

# ----- Work ------

varname = defVarmme('salinity'); v = 'S'

multStd = 2. # detect ToE at multStd std dev of histNat

# -- Choose random file to read only the basic variables and properties common to all files
file = 'cmip5.' + models[0]['name'] + '.1pctCO2.ensm.an.ocn.Omon.density.ver-' + models[0]['file_end_CO2'] + '_zon2D.nc'
f = open_ncfile(indir_CO2 + file,'r')

lat = f.variables['latitude'][:]; latN = lat.size
density = f.variables['lev'][:]; levN = density.size
time = f.variables['time'][:]; timN = time.size
var = varname['var_zonal_w/bowl']
basinN = 4

# Define variable properties
legVar = varname['legVar']
unit = varname['unit']

# ----- Compute zonal ToE for each model ------

for i, model in enumerate(models): # Loop on models

    print('Working on', model['name'])

    file_CO2 = glob.glob(indir_CO2+ 'cmip5.' + model['name'] + '.' + '*zon2D.nc')[0]
    file_piC = glob.glob(indir_piC+ 'cmip5.' + model['name'] + '.' + '*zon2D.nc')[0]
    fCO2 = open_ncfile(file_CO2,'r')
    fpiC = open_ncfile(file_piC,'r')

    # -- Read var 1pctCO2
    varCO2 = fCO2.variables[var][:]
    # -- Read var PiControl
    varpiC = fpiC.variables[var][-140:,:,:,:]
    # Read std of PiControl
    stdvarpiC = np.ma.std(varpiC, axis=0)

    # Save signal sign at the end
    varsignal_end = np.ma.masked_all((basinN,levN,latN))
    varsignal_end[1,:,:] = np.ma.average(varCO2[-5:,1,:,:],axis=0) - np.ma.average(varpiC[-5:,1,:,:],axis=0)
    varsignal_end[2,:,:] = np.ma.average(varCO2[-5:,2,:,:],axis=0) - np.ma.average(varpiC[-5:,2,:,:],axis=0)
    varsignal_end[3,:,:] = np.ma.average(varCO2[-5:,3,:,:],axis=0) - np.ma.average(varpiC[-5:,3,:,:],axis=0)

    # Reorganise i,j dims in single dimension data (speeds up loops)
    stdvarpiC_a = np.reshape(stdvarpiC[1,:,:], (levN*latN))
    stdvarpiC_p = np.reshape(stdvarpiC[2,:,:], (levN*latN))
    stdvarpiC_i = np.reshape(stdvarpiC[3,:,:], (levN*latN))
    varsignal_a = np.reshape(varCO2[:,1,:,:]-varpiC[:,1,:,:], (timN, levN*latN))
    varsignal_p = np.reshape(varCO2[:,2,:,:]-varpiC[:,2,:,:], (timN, levN*latN))
    varsignal_i = np.reshape(varCO2[:,3,:,:]-varpiC[:,3,:,:], (timN, levN*latN))

    # Initialize toe for each basin (density, lat)
    toe1_a = np.ma.masked_all((levN,latN))
    toe1_p = np.ma.masked_all((levN,latN))
    toe1_i = np.ma.masked_all((levN,latN))
    toe2_a = np.ma.masked_all((levN,latN))
    toe2_p = np.ma.masked_all((levN,latN))
    toe2_i = np.ma.masked_all((levN,latN))
    # Initialize output variable
    varToE1 = np.ma.masked_all((basinN,levN,latN)) # (>1std) (basin,density,latitude)
    varToE2 = np.ma.masked_all((basinN,levN,latN)) # (>2std) (basin,density,latitude)

    # Compute ToE as last date when diff 1%CO2 - PiControl is larger than mult * stddev
    toe1_a = np.reshape(findToE(varsignal_a, stdvarpiC_a, 1),(levN,latN))
    toe1_p = np.reshape(findToE(varsignal_p, stdvarpiC_p, 1),(levN,latN))
    toe1_i = np.reshape(findToE(varsignal_i, stdvarpiC_i, 1),(levN,latN))
    toe2_a = np.reshape(findToE(varsignal_a, stdvarpiC_a, multStd),(levN,latN))
    toe2_p = np.reshape(findToE(varsignal_p, stdvarpiC_p, multStd),(levN,latN))
    toe2_i = np.reshape(findToE(varsignal_i, stdvarpiC_i, multStd),(levN,latN))

    # Save in output variable
    varToE1[1,:,:] = toe1_a
    varToE1[2,:,:] = toe1_p
    varToE1[3,:,:] = toe1_i
    varToE2[1,:,:] = toe2_a
    varToE2[2,:,:] = toe2_p
    varToE2[3,:,:] = toe2_i

    fileName = 'cmip5.'+model['name']+'.toe_zonal_1pctCO2vsPiC.nc'
    dir = '/home/ysilvy/Density_bining/Yona_analysis/data/toe_zonal/toe_1pctCO2_piC/'
    description = 'Time of Emergence 1%CO2 vs. PiControl in latitude/density space for the three oceanic basins. ' \
    'The signal in each point is 1pctCO2-PiControl, and the noise is the standard deviation of the pre-industrial control. ' \
    'The ToE is computed by using once or twice the noise as the threshold.'

    fout = open_ncfile(dir+fileName,'w', format='NETCDF4')
    fout.description = description

    # dimensions
    fout.createDimension('basin', 4)
    fout.createDimension('density', levN)
    fout.createDimension('latitude', latN)

    # variables
    basin = fout.createVariable('basin', 'i4', ('basin',))
    density = fout.createVariable('density', 'd', ('density',))
    latitude = fout.createVariable('latitude', 'f4', ('latitude',))
    ToE2 = fout.createVariable(varname['var_zonal_w/bowl']+'ToE2', 'f4', ('basin','density','latitude',))
    ToE1 = fout.createVariable(varname['var_zonal_w/bowl']+'ToE1', 'f4', ('basin','density','latitude',))
    varchange = fout.createVariable(varname['var_zonal_w/bowl']+'_change', 'f4', ('basin','density','latitude',))

    # data
    basin[:] =  np.arange(0,basinN)
    density[:] = density
    latitude[:] = latitude
    ToE2[:,:,:] = varToE2
    ToE1[:,:,:] = varToE1
    varchange[:,:,:] = varsignal_end

    # units
    basin.units = 'basin index'
    density.units = 'kg.m-3'
    latitude.units = ''
    ToE2.units = 'Year'
    ToE2.long_name = 'ToE (>2std) for ' + legVar
    ToE1.units = 'Year'
    ToE1.long_name = 'ToE (>1std) for ' + legVar
    varchange.units = unit
    varchange.longName = legVar + ' signal averaged in the last five years'

    fout.close()
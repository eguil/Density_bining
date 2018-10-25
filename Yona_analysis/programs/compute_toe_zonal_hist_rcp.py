#!/bin/env python
# -*- coding: utf-8 -*-

"""
Compute ToE hist + RCP85 vs. histNat (or PiControl) in lat/rho domain for all runs of available models
Over the 2005-2100 period, use either RCP85-average(histNat) or RCP85-average(PiControl) as the signal
Save ToE in output files
For now only the first option is implemented
"""

import numpy as np
from netCDF4 import Dataset as open_ncfile
from maps_matplot_lib import defVarmme
from modelsDef import defModels
from libToE import findToE, findToE_2thresholds
import glob

# ----- Workspace ------

indir_histrcp85 = '/home/ysilvy/Density_bining/Yona_analysis/data/hist_rcp85/'
indir_histNat = '/data/ericglod/Density_binning/Prod_density_april15/mme_histNat/'
indir_piC = '/data/ericglod/Density_binning/Prod_density_april15/mme_piControl/'

models = defModels()

# ----- Work ------

varname = defVarmme('salinity'); v = 'S'

multStd = 2. # detect ToE at multStd std dev of histNat

use_piC = False # Over projection period, signal = RCP-average(histNat), noise = std(histNat)
# use_piC = True # Over projection period, signal = RCP-average(PiControl), noise = std(PiControl)

iniyear = 1860
finalyear = 2100
deltay = 10.

# Choose random file to read only the basic variables and properties common to all files
file = 'cmip5.' + models[1]['name'] + '.historicalNat.ensm.an.ocn.Omon.density.ver-' + \
       models[1]['file_end_histNat'] + '_zon2D.nc'
f = open_ncfile(indir_histNat + file,'r')

lat = f.variables['latitude'][:]; latN = lat.size
density = f.variables['lev'][:]; levN = density.size
timN = 240
var = varname['var_zonal_w/bowl']
basinN = 4

# Define variable properties
legVar = varname['legVar']

# ----- Compute zonal ToE for each simulation ------

nMembers = np.ma.zeros(len(models)) # Initialize array for keeping nb of members per model

for i, model in enumerate(models): # Loop on models

    print('Working on', model['name'])

    if (use_piC != True) or (use_piC == True and model['name'] != 'GISS-E2-R' and model['name'] != 'FGOALS-g2'
                             and model['name'] != 'MIROC-ESM'):

        # Read hist+rcp85 files
        listruns = glob.glob(indir_histrcp85 + 'cmip5.' + model['name'] + '.' + '*zon2D.nc')
        nruns = len(listruns)
        nMembers[i] = nruns
        if nruns != 0:
            # Index of common time interval
            tstart = model['props'][2]
            tend = model['props'][3]

            # Read histNat ensemble mean
            filehn = 'cmip5.' + model['name'] + '.historicalNat.ensm.an.ocn.Omon.density.ver-' \
                     + model['file_end_histNat'] + '_zon2D.nc'
            fhn = open_ncfile(indir_histNat + filehn,'r')

            # Read var histNat
            varhn_a = fhn.variables[var][:,1,:,:].squeeze()
            varhn_p = fhn.variables[var][:,2,:,:].squeeze()
            varhn_i = fhn.variables[var][:,3,:,:].squeeze()

            # Read std of histNat (max std of all runs for each model)
            stdvarhn_a = fhn.variables[var+'Std'][1,:,:].squeeze()
            stdvarhn_i = fhn.variables[var+'Std'][3,:,:].squeeze()
            stdvarhn_p = fhn.variables[var+'Std'][2,:,:].squeeze()

            # Compute time average of the whole histNat series (signal over projection = RCP - mean(histNat))
            meanvarhn_a = np.ma.mean(varhn_a, axis=0)
            meanvarhn_p = np.ma.mean(varhn_p, axis=0)
            meanvarhn_i = np.ma.mean(varhn_i, axis=0)

            # Reorganise i,j dims in single dimension data (speeds up loops)
            stdvarhn_a = np.reshape(stdvarhn_a, (levN*latN))
            stdvarhn_p = np.reshape(stdvarhn_p, (levN*latN))
            stdvarhn_i = np.reshape(stdvarhn_i, (levN*latN))

            if use_piC:
                # Read and Compute time average of PiControl over last 95 years + std of PiControl
                filepiC = glob.glob(indir_piC + 'cmip5.' + model['name'] + '.' + '*zon2D.nc')[0]
                fpiC = open_ncfile(filepiC,'r')
                varpiC_a = fpiC.variables[var][-95:,1,:,:].squeeze()
                varpiC_p = fpiC.variables[var][-95:,2,:,:].squeeze()
                varpiC_i = fpiC.variables[var][-95:,3,:,:].squeeze()
                meanvarpiC_a = np.ma.mean(varpiC_a, axis=0)
                meanvarpiC_p = np.ma.mean(varpiC_p, axis=0)
                meanvarpiC_i = np.ma.mean(varpiC_i, axis=0)
                stdvarpiC_a = np.ma.std(varpiC_a, axis=0)
                stdvarpiC_p = np.ma.std(varpiC_p, axis=0)
                stdvarpiC_i = np.ma.std(varpiC_i, axis=0)

                # Reorganise i,j dims in single dimension data (speeds up loops)
                stdvarpiC_a = np.reshape(stdvarpiC_a, (levN*latN))
                stdvarpiC_p = np.reshape(stdvarpiC_p, (levN*latN))
                stdvarpiC_i = np.reshape(stdvarpiC_i, (levN*latN))

            # Initialize toe for each basin (run, density, lat)
            toe_a = np.ma.masked_all((nruns,levN,latN))
            toe_p = np.ma.masked_all((nruns,levN,latN))
            toe_i = np.ma.masked_all((nruns,levN,latN))
            # Initialize output variable
            varToE = np.ma.masked_all((nruns,basinN,levN,latN)) # (members,basin,density,latitude)

            # Loop over number of runs
            for k in range(nruns):
                print('    . run number', k)

                # Read file
                fhrcp = open_ncfile(listruns[k],'r')
                # Read var hist + rcp85
                varh_a = fhrcp.variables[var][tstart:tend,1,:,:].squeeze()
                varh_p = fhrcp.variables[var][tstart:tend,2,:,:].squeeze()
                varh_i = fhrcp.variables[var][tstart:tend,3,:,:].squeeze()
                varrcp_a = fhrcp.variables[var][tend:tend+95,1,:,:].squeeze()
                varrcp_p = fhrcp.variables[var][tend:tend+95,2,:,:].squeeze()
                varrcp_i = fhrcp.variables[var][tend:tend+95,3,:,:].squeeze()

                # Initialize and fill var_signal for each basin (timN, density, latitude)
                varsignal_a = np.ma.masked_all((timN,levN,latN))
                varsignal_p = np.ma.masked_all((timN,levN,latN))
                varsignal_i = np.ma.masked_all((timN,levN,latN))

                varsignal_a[0:145,:,:] = varh_a-varhn_a
                varsignal_p[0:145,:,:] = varh_p-varhn_p
                varsignal_i[0:145,:,:] = varh_i-varhn_i

                if use_piC == False:
                    varsignal_a[145:,:,:] = varrcp_a-meanvarhn_a
                    varsignal_p[145:,:,:] = varrcp_p-meanvarhn_p
                    varsignal_i[145:,:,:] = varrcp_i-meanvarhn_i
                else:
                    varsignal_a[145:,:,:] = varrcp_a-meanvarpiC_a
                    varsignal_p[145:,:,:] = varrcp_p-meanvarpiC_p
                    varsignal_i[145:,:,:] = varrcp_i-meanvarpiC_i

                # Reorganise i,j dims in single dimension data (speeds up loops)
                varsignal_a = np.reshape(varsignal_a, (timN, levN*latN))
                varsignal_p = np.reshape(varsignal_p, (timN, levN*latN))
                varsignal_i = np.reshape(varsignal_i, (timN, levN*latN))

                # Compute ToE as last date when diff hist+RCP - histNat is larger than mult * stddev
                toe_a = np.reshape(findToE(varsignal_a, stdvarhn_a, multStd),(levN,latN))
                toe_p = np.reshape(findToE(varsignal_p, stdvarhn_p, multStd),(levN,latN))
                toe_i = np.reshape(findToE(varsignal_i, stdvarhn_i, multStd),(levN,latN))

                # Save in output variable
                varToE[k,1,:,:] = toe_a
                varToE[k,2,:,:] = toe_p
                varToE[k,3,:,:] = toe_i

            # Save in output file
            if use_piC == False:
                fileName = 'cmip5.'+model['name']+'.toe_zonal_rcp_histNat.nc'
                dir = '/home/ysilvy/Density_bining/Yona_analysis/data/toe_zonal/toe_rcp85_histNat/'
                description = 'Time of Emergence hist+rcp8.5 vs. histNat for each member. \n' \
                              'The historical runs are prolonged by the 95 years of RCP8.5. ' \
                              'The ensemble mean historicalNat is used here for all historical runs of the model. \n' \
                              'The signal in each point is hist-histNat over the historical period and ' \
                              'hist-timeaverage(histNat) over the projection period. ' \
                              'The noise in each point is the max standard deviation of the historicalNat runs \n' \
                              'The ToE is computed by using twice the noise as the threshold.'
            else:
                fileName = 'cmip5.'+model['name']+'.toe_zonal_rcp_PiControl.nc'
                dir = '/home/ysilvy/Density_bining/Yona_analysis/data/toe_zonal/toe_rcp85_PiControl/'
                description = 'Time of Emergence hist vs. histNat +rcp8.5 vs. PiControl for each member. \n' \
                              'The historical runs are prolonged by the 95 years of RCP8.5. ' \
                              'The ensemble mean historicalNat is used here for all historical runs of the model. \n' \
                              'The signal in each point is hist-histNat over the historical period and ' \
                               'hist-timeaverage(PiControl) over the projection period. ' \
                              'The noise in each point is the max standard deviation of the historicalNat runs over ' \
                              'the historical period and the standard deviation of the PiControl run over the projection period. \n' \
                              'The ToE is computed by using twice the noise as the threshold.'

            fout = open_ncfile(dir+fileName,'w', format='NETCDF4')
            fout.description = description

            # dimensions
            fout.createDimension('members', nruns)
            fout.createDimension('basin', 4)
            fout.createDimension('density', levN)
            fout.createDimension('latitude', latN)

            # variables
            members = fout.createVariable('members', 'f4', ('members',))
            basin = fout.createVariable('basin', 'i4', ('basin',))
            density = fout.createVariable('density', 'd', ('density',))
            latitude = fout.createVariable('latitude', 'f4', ('latitude',))
            varToE2 = fout.createVariable(varname['var_zonal_w/bowl']+'ToE2', 'f4', ('members','basin','density','latitude',))

            # data
            members[:] =  np.arange(0,nruns)
            basin[:] =  np.arange(0,basinN)
            density[:] = density
            latitude[:] = latitude
            varToE2[:,:,:] = varToE

            # units
            basin.units = 'basin index'
            density.units = 'kg.m-3'
            latitude.units = ''
            varToE2.units = 'Year'
            varToE2.long_name = 'ToE 2 for ' + legVar

            fout.close()


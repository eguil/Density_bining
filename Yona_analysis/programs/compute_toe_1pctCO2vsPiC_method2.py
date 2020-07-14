#!/bin/env python
# -*- coding: utf-8 -*-

"""
Average signal (1pctCO2 - PiControl) and noise (2*std(PiControl)) in the specified domains,
then compute ToE  1pctCO2 vs. Pi Control.
Save ToE in output files
(Method 2 = 'average_signal' compared to method 1 'average_toe' which first calculated the ToE everywhere 
then averaged it in the specified domains, making the calculations much longer)
"""

import numpy as np
from netCDF4 import Dataset as open_ncfile
from maps_matplot_lib import defVarmme, averageDom
from modelsDef import defModelsCO2piC
from libToE import findToE, ToEdomain1pctCO2vsPiC

# ----- Workspace ------

indir_1pctCO2 = '/data/ericglod/Density_binning/Prod_density_april15/mme_1pctCO2/'
indir_piC = '/data/ericglod/Density_binning/Prod_density_april15/mme_piControl/'

models = defModelsCO2piC()

# ----- Work ------

varname = defVarmme('salinity'); v = 'S'

method = 'average_signal'

# -- Choose which 'noise' to use for the ToE calculation
# method_noise = 'average_std' # Average the standard deviation of PiC in the specified domains
method_noise = 'average_piC' # Average PiC in the specified domains then determine the std of this averaged value

domains = ['Southern ST', 'SO', 'Northern ST', 'North Atlantic', 'North Pacific']
signal_domains = ['fresher','saltier','fresher','saltier','saltier']

multStd = 2. # detect ToE at multStd std dev of PiControl

iniyear = 0
finalyear = 140
deltay = 10.


# ----- Variables ------

# -- Choose random file to read only the basic variables and properties common to all files
file = 'cmip5.' + models[0]['name'] + '.1pctCO2.ensm.an.ocn.Omon.density.ver-' + models[0]['file_end_CO2'] + '_zon2D.nc'
f = open_ncfile(indir_1pctCO2 + file,'r')

lat = f.variables['latitude'][:]; latN = lat.size
density = f.variables['lev'][:]; levN = density.size
time = f.variables['time'][:]; timN = time.size
basinN = 4
var = varname['var_zonal']

# -- Define variable properties
legVar = varname['legVar']
unit = varname['unit']


# ----- Average signal and noise and compute ToE for each model ------

for i, model in enumerate(models):

    print('Working on', model['name'])

    # -- Read 1pctCO2 file and piControl file
    file_1pctCO2 = 'cmip5.' + model['name'] + '.1pctCO2.ensm.an.ocn.Omon.density.ver-' + model['file_end_CO2'] + '_zon2D.nc'
    file_piC = 'cmip5.' + model['name'] + '.piControl.ensm.an.ocn.Omon.density.ver-' + model['file_end_piC'] + '_zon2D.nc'
    fCO2 = open_ncfile(indir_1pctCO2 + file_1pctCO2,'r')
    fpiC = open_ncfile(indir_piC + file_piC,'r')

    # -- Read var 1pctCO2
    varCO2 = fCO2.variables[var][:]
    # -- Read var PiControl
    varpiC = fpiC.variables[var][-240:,:,:,:]

    # -- Compute std of PiControl
    varstd = np.ma.std(varpiC, axis=0)
    
    # -- Compute mean of PiControl
    meanvarpiC = np.ma.mean(varpiC,axis=0)

    # -- Initalize varnoise for each basin, containing averaged noise for each domain
    varnoise_a = np.ma.masked_all(len(domains))
    varnoise_p = np.ma.masked_all(len(domains))
    varnoise_i = np.ma.masked_all(len(domains))
    # -- Initialize varsignal for each basin, containing averaged signal for each domain
    varsignal_a = np.ma.masked_all((timN,len(domains)))
    varsignal_p = np.ma.masked_all((timN,len(domains)))
    varsignal_i = np.ma.masked_all((timN,len(domains)))
    # -- Initialize toe for each basin
    toe2_a = np.ma.masked_all(len(domains)) # >1std
    toe2_p = np.ma.masked_all(len(domains))
    toe2_i = np.ma.masked_all(len(domains))
    toe1_a = np.ma.masked_all(len(domains)) # >1std
    toe1_p = np.ma.masked_all(len(domains))
    toe1_i = np.ma.masked_all(len(domains))
    # -- Initialize ouput variable
    varToE1 = np.ma.masked_all((basinN, len(domains)))
    varToE2 = np.ma.masked_all((basinN, len(domains)))

    # -- Loop over domains
    for j, domain_name in enumerate(domains):
        print('- ' + domain_name)

        # Select domain to average
        domain = ToEdomain1pctCO2vsPiC(model['name'], domain_name)[0]

        # Average signal and noise
        if domain['Atlantic'] != None:
            varsignal_a[:,j] = averageDom(varCO2[:,1,:,:]-meanvarpiC[1,:,:], 3, domain['Atlantic'], lat, density)
            if method_noise == 'average_std':
                varnoise_a[j] = averageDom(varstd[1,:,:], 2, domain['Atlantic'], lat, density)
            else:
                varnoise_a[j] = np.ma.std(averageDom(varpiC[:,1,:,:], 3, domain['Atlantic'], lat, density), axis=0)
        if domain['Pacific'] != None:
            varsignal_p[:,j] = averageDom(varCO2[:,2,:,:]-meanvarpiC[2,:,:], 3, domain['Pacific'], lat, density)
            if method_noise == 'average_std':
                varnoise_p[j] = averageDom(varstd[2,:,:], 2, domain['Pacific'], lat, density)
            else:
                varnoise_p[j] = np.ma.std(averageDom(varpiC[:,2,:,:], 3, domain['Pacific'], lat, density), axis=0)
        if domain['Indian'] != None:
            varsignal_i[:,j] = averageDom(varCO2[:,3,:,:]-meanvarpiC[3,:,:], 3, domain['Indian'], lat, density)
            if method_noise == 'average_std':
                varnoise_i[j] = averageDom(varstd[3,:,:], 2, domain['Indian'], lat, density)
            else:
                varnoise_i[j] = np.ma.std(averageDom(varpiC[:,3,:,:], 3, domain['Indian'], lat, density), axis=0)

        # Compute ToE of averaged domain
        if domain['Atlantic'] != None and np.ma.is_masked(varnoise_a[j]) == False:
            toe1_a[j] = findToE(varsignal_a[:,j], varnoise_a[j], 1)
            toe2_a[j] = findToE(varsignal_a[:,j], varnoise_a[j], multStd)
            print(toe1_a[j],toe2_a[j])
        if domain['Pacific'] != None and np.ma.is_masked(varnoise_p[j]) == False:
            toe1_p[j] = findToE(varsignal_p[:,j], varnoise_p[j], 1)
            toe2_p[j] = findToE(varsignal_p[:,j], varnoise_p[j], multStd)
        if domain['Indian'] != None and np.ma.is_masked(varnoise_i[j]) == False:
            toe1_i[j] = findToE(varsignal_i[:,j], varnoise_i[j], 1)
            toe2_i[j] = findToE(varsignal_i[:,j], varnoise_i[j], multStd)

        # Take out runs where the signal is of opposite sign than expected
        if signal_domains[j] == 'fresher':
            if np.ma.mean(varsignal_a[-5:,j],axis=0) > multStd * varnoise_a[j]:
                toe2_a[j] = np.ma.masked
            if np.ma.mean(varsignal_a[-5:,j],axis=0) > 1 * varnoise_a[j]:
                toe1_a[j] = np.ma.masked
            if np.ma.mean(varsignal_p[-5:,j],axis=0) > multStd * varnoise_p[j]:
                toe2_p[j] = np.ma.masked
            if np.ma.mean(varsignal_p[-5:,j],axis=0) > 1 * varnoise_p[j]:
                toe1_p[j] = np.ma.masked
            if np.ma.mean(varsignal_i[-5:,j],axis=0) > multStd * varnoise_i[j]:
                toe2_i[j] = np.ma.masked
            if np.ma.mean(varsignal_i[-5:,j],axis=0) > 1 * varnoise_i[j]:
                toe1_i[j] = np.ma.masked
        else:
            if np.ma.mean(varsignal_a[-5:,j],axis=0) < -multStd * varnoise_a[j]:
                toe2_a[j] = np.ma.masked
            if np.ma.mean(varsignal_a[-5:,j],axis=0) < -1 * varnoise_a[j]:
                toe1_a[j] = np.ma.masked
            if np.ma.mean(varsignal_p[-5:,j],axis=0) < -multStd * varnoise_p[j]:
                toe2_p[j] = np.ma.masked
            if np.ma.mean(varsignal_p[-5:,j],axis=0) < -1 * varnoise_p[j]:
                toe1_p[j] = np.ma.masked
            if np.ma.mean(varsignal_i[-5:,j],axis=0) < -multStd * varnoise_i[j]:
                toe2_i[j] = np.ma.masked
            if np.ma.mean(varsignal_i[-5:,j],axis=0) < -1 * varnoise_i[j]:
                toe1_i[j] = np.ma.masked
        print(toe1_a[j], toe2_a[j])

    varToE1[1,:] = toe1_a
    varToE1[2,:] = toe1_p
    varToE1[3,:] = toe1_i
    varToE2[1,:] = toe2_a
    varToE2[2,:] = toe2_p
    varToE2[3,:] = toe2_i

    print(varToE1[1,:],varToE2[1,:])


    # Save in output file
    fileName = 'cmip5.'+model['name']+'.'+legVar+'_toe_1pctCO2vsPiControl_method2_'+method_noise+'.nc'
    if method_noise == 'average_std':
        dir = '/home/ysilvy/Density_bining/Yona_analysis/data/toe_1pctCO2vsPiC_average_signal/average_std/'
    else:
        dir = '/home/ysilvy/Density_bining/Yona_analysis/data/toe_1pctCO2vsPiC_average_signal/average_piC/'
    fout = open_ncfile(dir+fileName,'w', format='NETCDF4')
    if method_noise == 'average_std':
        noise_description = 'Noise is computed by averaging the standard deviation of the pre-industrial control runs ' \
                            'in the specified domains.'
    else :
        noise_description = 'Noise is computed by averaging the pre-industrial control runs in the specified domains,' \
                            ' and then taking the standard deviation of the average.'
    fout.description = 'ToE 1%CO2 vs. PiControl, in 5 domains : Southern Subtropics (0), Southern Ocean (1),' \
                        ' Northern Subtropics (2), North Atlantic (3), North Pacific (4). Signal is computed by averaging ' \
                       'the difference 1pctCO2 - timeaverage(PiControl_last240years) in those domains. ' + noise_description + ' Then ToE is computed ' \
                        'using once or twice the noise as the limit.'

    # dimensions
    fout.createDimension('basin', 4)
    fout.createDimension('domain', 5)

    # variables
    basin = fout.createVariable('basin', 'i4', ('basin',))
    domain = fout.createVariable('domain', 'f4', ('domain',))
    ToE2 = fout.createVariable(varname['var_zonal_w/bowl']+'ToE2', 'f4', ('basin','domain',))
    ToE1 = fout.createVariable(varname['var_zonal_w/bowl']+'ToE1', 'f4', ('basin','domain',))

    # data
    basin[:] =  np.arange(0,basinN)
    domain[:] = np.arange(0,len(domains))
    ToE2[:,:] = varToE2
    ToE1[:,:] = varToE1

    # units
    basin.units = 'basin index'
    domain.units = 'domain index'
    ToE2.units = 'Year'
    ToE2.long_name = 'ToE (>2std) for ' + legVar
    ToE1.units = 'Year'
    ToE1.long_name = 'ToE (>1std) for ' + legVar

    fout.close()

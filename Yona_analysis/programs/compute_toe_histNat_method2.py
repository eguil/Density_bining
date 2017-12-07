#!/bin/env python
# -*- coding: utf-8 -*-

"""
Average signal and noise for each run of each model in the specified domains, then compute ToE hist vs. histNat
Save ToE in output files
"""

import numpy as np
from netCDF4 import Dataset as open_ncfile
from maps_matplot_lib import defVarmme, averageDom
from modelsDef import defModels
from libToE import findToE, ToEdomainhistvshistNat
import glob
import os


# ----- Workspace ------

indir_hist = '/data/ericglod/Density_binning/Prod_density_april15/historical/'
indir_histNat = '/data/ericglod/Density_binning/Prod_density_april15/mme_histNat/'
indir_noise = '/home/ysilvy/Density_bining/Yona_analysis/data/noise_estimate/'

models = defModels()

# ----- Work ------

varname = defVarmme('salinity'); v = 'S'
#varname = defVarmme('temp'); v = 'T'
#varname = defVarmme('depth'); v = 'Z'

method = 'average_signal' # Average signal and noise in the box, then compute ToE

# -- Choose which 'noise' to use for the ToE calculation
method_noise = 'average_std' # Average the standard deviation of histNat in the specified domains
# method_noise = 'average_histNat' # Average histNat in the specified domains then determine the std of this averaged value
print('Noise computation: %s' %(method_noise,))

domains = ['Southern ST', 'SO', 'Northern ST', 'North Atlantic', 'North Pacific']

multStd = 2. # detect ToE at multStd std dev of histNat

iniyear = 1860
finalyear = 2005
deltay = 10.

# ----- Variables ------

# Choose random file to read only the basic variables and properties common to all files
file = 'cmip5.' + models[0]['name'] + '.historicalNat.ensm.an.ocn.Omon.density.ver-' + models[0]['file_end_histNat'] + '_zon2D.nc'
f = open_ncfile(indir_histNat + file,'r')

lat = f.variables['latitude'][:]; latN = lat.size
density = f.variables['lev'][:]; levN = density.size
time = f.variables['time'][:]; timN = time.size
var = varname['var_zonal_w/bowl']
basinN = 4

# Define variable properties
legVar = varname['legVar']


# ----- Average signal and noise and compute ToE for each simulation ------

nMembers = np.ma.empty(len(models)) # Initialize array for keeping nb of members per model

for i, model in enumerate(models):

    print('Working on', model['name'])

    # Read histNat ensemble mean
    filehn = 'cmip5.' + model['name'] + '.historicalNat.ensm.an.ocn.Omon.density.ver-' + model['file_end_histNat'] + '_zon2D.nc'
    fhn = open_ncfile(indir_histNat + filehn,'r')
    # Read var histNat
    varhn_a = fhn.variables[var][0:145,1,:,:].squeeze()
    varhn_p = fhn.variables[var][0:145,2,:,:].squeeze()
    varhn_i = fhn.variables[var][0:145,3,:,:].squeeze()
    # Read std of histNat (max std of all runs for each model)
    if method_noise == 'average_std':
        varstd_a = fhn.variables[var+'Std'][1,:,:].squeeze()
        varstd_p = fhn.variables[var+'Std'][2,:,:].squeeze()
        varstd_i = fhn.variables[var+'Std'][3,:,:].squeeze()

    # Read hist files
    listruns = glob.glob(indir_hist + 'cmip5.' + model['name'] + '.' + '*zon2D.nc')
    nruns = len(listruns)
    nMembers[i] = nruns
    # Index of common time interval
    tstart = model['props'][2]
    tend = model['props'][3]

    # Initialize varnoise for each basin, containing averaged noise for each domain
    varnoise_a = np.ma.masked_all(len(domains))
    varnoise_p = np.ma.masked_all(len(domains))
    varnoise_i = np.ma.masked_all(len(domains))
    # Initialize varsignal for each basin, containing averaged signal for each domain and each run
    varsignal_a = np.ma.masked_all((timN,nruns,len(domains)))
    varsignal_p = np.ma.masked_all((timN,nruns,len(domains)))
    varsignal_i = np.ma.masked_all((timN,nruns,len(domains)))
    # Initialize toe for each basin (run, domain)
    toe_a = np.ma.masked_all((nruns,len(domains)))
    toe_p = np.ma.masked_all((nruns,len(domains)))
    toe_i = np.ma.masked_all((nruns,len(domains)))
    # Initialize output variable
    varToE = np.ma.masked_all((nruns,basinN,len(domains))) # (members,basin,domain)
    # Initialize array to store run names
    run_names = np.empty(nruns).astype('S6')

    if method_noise == 'average_histNat':
        # Here we read the std of the averaged histNat in the 5 domains for all runs, then take the max as our noise
        filenoise = 'cmip5.' + model['name'] + '.noise_domains_hist_histNat.std_of_average.nc'
        fnoise = open_ncfile(indir_noise + filenoise,'r')
        varstd = fnoise.variables[var+'stdhn'][:]
        varnoise_a = np.ma.max(varstd[:,1,:], axis=0)
        varnoise_p = np.ma.max(varstd[:,2,:], axis=0)
        varnoise_i = np.ma.max(varstd[:,3,:], axis=0)

    # Loop over 5 domains
    for j, domain_name in enumerate(domains):
        print('- ', j, domains[j])

        # Select domain to average
        domain = ToEdomainhistvshistNat(model['name'], domain_name)[0]
        domain_char = ToEdomainhistvshistNat(model['name'], domain_name)[1]

        if method_noise == 'average_std':
            # Average noise
            if domain['Atlantic'] != None:
                varnoise_a[j] = averageDom(varstd_a, 2, domain['Atlantic'], lat, density)
            if domain['Pacific'] != None:
                varnoise_p[j] = averageDom(varstd_p, 2, domain['Pacific'], lat, density)
            if domain['Indian'] != None:
                varnoise_i[j] = averageDom(varstd_i, 2, domain['Indian'], lat, density)

        # Loop over number of runs
        for k in range(nruns):
            run_number = os.path.basename(listruns[k]).split('.')[3]
            run_names[k] = run_number
            print('    . run number', k, run_number)

            # Read file
            fh = open_ncfile(listruns[k],'r')
            # Read var hist
            varh_a = fh.variables[var][tstart:tend,1,:,:].squeeze()
            varh_p = fh.variables[var][tstart:tend,2,:,:].squeeze()
            varh_i = fh.variables[var][tstart:tend,3,:,:].squeeze()

            # Average signal var hist - var histNat
            if domain['Atlantic'] != None:
                varsignal_a[:,k,j] = averageDom(varh_a-varhn_a, 3, domain['Atlantic'], lat, density)
            if domain['Pacific'] != None:
                varsignal_p[:,k,j] = averageDom(varh_p-varhn_p, 3, domain['Pacific'], lat, density)
            if domain['Indian'] != None:
                varsignal_i[:,k,j] = averageDom(varh_i-varhn_i, 3, domain['Indian'], lat, density)

            print('      varsignal shape:', varsignal_a.shape, varsignal_p.shape, varsignal_i.shape)

            # Compute ToE of averaged domain for run k
            if domain['Atlantic'] != None and np.ma.is_masked(varnoise_a[j]) == False:
                toe_a[k,j] = findToE(varsignal_a[:,k,j], varnoise_a[j], multStd) + iniyear
            if domain['Pacific'] != None and np.ma.is_masked(varnoise_p[j]) == False:
                toe_p[k,j] = findToE(varsignal_p[:,k,j], varnoise_p[j], multStd) + iniyear
            if domain['Indian'] != None and np.ma.is_masked(varnoise_i[j]) == False:
                toe_i[k,j] = findToE(varsignal_i[:,k,j], varnoise_i[j], multStd) + iniyear

    varToE[:,1,:] = toe_a
    varToE[:,2,:] = toe_p
    varToE[:,3,:] = toe_i
    print('  varToE shape:', varToE.shape)
    print('  ', np.ma.around(np.ma.median(varToE[:,1,0])))
    print('')


    # Save in output file
    fileName = 'cmip5.'+model['name']+'.toe_histNat_method2_'+method_noise+'.nc'
    dir = '/home/ysilvy/Density_bining/Yona_analysis/data/toe_histNat_average_signal/'+method_noise+'/'
    fout = open_ncfile(dir+fileName,'w', format='NETCDF4')
    if method_noise == 'average_std':
        noise_description = 'Noise is computed by averaging the max standard deviation of the historical Nat runs ' \
                            'in the specified domains.'
    else :
        noise_description = 'Noise is computed by averaging the historical Nat runs in the specified domains,' \
                            ' then determining the standard deviation of the average, and then by keeping only the maximum ' \
                            'standard deviation of all runs.'
    fout.description = 'ToE hist vs. histNat for each member, in 5 domains : Southern Subtropics (0), Southern Ocean (1),' \
                        'Northern Subtropics (2), North Atlantic (3), North Pacific (4). Signal is computed by averaging ' \
                       'the difference historical - historicalNat in those domains. The ensemble mean historicalNat is ' \
                       'used here for all the historical runs of each model ' + noise_description + ' Then ToE is computed ' \
                        'using twice the noise as the limit.'

    # dimensions
    fout.createDimension('members', nruns)
    fout.createDimension('basin', 4)
    fout.createDimension('domain', 5)

    # variables
    members = fout.createVariable('members', 'f4', ('members',))
    run_name = fout.createVariable('run_name', 'str', ('members',))
    basin = fout.createVariable('basin', 'i4', ('basin',))
    domain = fout.createVariable('domain', 'f4', ('domain',))
    varToE2 = fout.createVariable(varname['var_zonal_w/bowl']+'ToE2', 'f4', ('members','basin','domain',))

    # data
    members[:] =  np.arange(0,nruns)
    run_name[:] = run_names
    basin[:] =  np.arange(0,basinN)
    domain[:] = np.arange(0,len(domains))
    varToE2[:,:,:] = varToE

    # units
    basin.units = 'basin index'
    domain.units = 'domain index'
    varToE2.units = 'Year'
    varToE2.long_name = 'ToE 2 for ' + legVar
    run_name.units = 'Name of the runs'

    fout.close()


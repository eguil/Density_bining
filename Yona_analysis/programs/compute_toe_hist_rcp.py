#!/bin/env python
# -*- coding: utf-8 -*-

"""
Average signal and noise for each run of each model in the specified domains, then compute ToE hist + RCP85
Save ToE in output files
"""

import numpy as np
from netCDF4 import Dataset as open_ncfile
from maps_matplot_lib import defVarmme, averageDom
from modelsDef import defModels
from libToE import findToE, ToEdomainhistvshistNat
import glob

# ----- Workspace ------

indir_histrcp85 = '/home/ysilvy/Density_bining/Yona_analysis/data/hist_rcp85/'
indir_histNat = '/data/ericglod/Density_binning/Prod_density_april15/mme_histNat/'

models = defModels()

# ----- Work ------

varname = defVarmme('salinity'); v = 'S'

method = 'average_signal' # Average signal and noise in the box, then compute ToE

# -- Choose which 'noise' to use for the ToE calculation
# method_noise = 'average_std' # Average the standard deviation of PiC in the specified domains
method_noise = 'average_histNat' # Average histNat in the specified domains then determine the std of this averaged value

domains = ['Southern ST', 'SO', 'Northern ST', 'North Atlantic', 'North Pacific']

multStd = 2. # detect ToE at multStd std dev of histNat

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

# ----- Average signal and noise and compute ToE for each simulation ------

nMembers = np.ma.empty(len(models)) # Initialize array for keeping nb of members per model
## Voir si on laisse l'initialisation ici ou dans la boucle au fur et à mesure du bon nombre
# de modèles rencontrés

for i, model in enumerate(models):

    print('Working on', model['name'])

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
        print(' varhn shape : ', varhn_a.shape)
        # Compute std of histNat (noise)
        varstda = np.ma.std(varhn_a, axis=0)
        varstdp = np.ma.std(varhn_p, axis=0)
        varstdi = np.ma.std(varhn_i, axis=0)
        # Compute time average of the whole histNat series (signal = hist+RCP - mean(histNat))
        meanvarhn_a = np.ma.mean(varhn_a, axis=0)
        meanvarhn_p = np.ma.mean(varhn_p, axis=0)
        meanvarhn_i = np.ma.mean(varhn_i, axis=0)

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

        # Loop over 5 domains
        for j, domain_name in enumerate(domains):
            print('- ', j, domains[j])

            # Select domain to average
            domain = ToEdomainhistvshistNat(model['name'], domain_name)[0]
            domain_char = ToEdomainhistvshistNat(model['name'], domain_name)[1]

            # Average noise
            if domain['Atlantic'] != None:
                if method_noise == 'average_std':
                    varnoise_a[j] = averageDom(varstda, 2, domain['Atlantic'], lat, density)
                else:
                    varnoise_a[j] = np.ma.std(averageDom(varhn_a, 3, domain['Atlantic'], lat, density), axis=0)
            if domain['Pacific'] != None:
                if method_noise == 'average_std':
                    varnoise_p[j] = averageDom(varstdp, 2, domain['Pacific'], lat, density)
                else:
                    varnoise_p[j] = np.ma.std(averageDom(varhn_p, 3, domain['Pacific'], lat, density), axis=0)
            if domain['Indian'] != None:
                if method_noise == 'average_std':
                    varnoise_i[j] = averageDom(varstdi, 2, domain['Indian'], lat, density)
                else:
                    varnoise_i[j] = np.ma.std(averageDom(varhn_i, 3, domain['Indian'], lat, density), axis=0)

            # Loop over number of runs
            for k in range(nruns):
                print('    . run number', k)

                # Read file
                fhrcp = open_ncfile(listruns[k],'r')
                # Read var hist + rcp85
                varh_a = fhrcp.variables[var][tstart:tend,1,:,:].squeeze()
                varh_p = fhrcp.variables[var][tstart:tend,2,:,:].squeeze()
                varh_i = fhrcp.variables[var][tstart:tend,3,:,:].squeeze()
                varrcp_a = fhrcp.variables[var][-95:,1,:,:].squeeze()
                varrcp_p = fhrcp.variables[var][-95:,2,:,:].squeeze()
                varrcp_i = fhrcp.variables[var][-95:,3,:,:].squeeze()
                # varhrcp_a = np.concatenate((varh_a,varrcp_a), axis=0)
                # varhrcp_p = np.concatenate((varh_p,varrcp_p), axis=0)
                # varhrcp_i = np.concatenate((varh_i,varrcp_i), axis=0)
                # print('       varhrcp shape : ', varhrcp_a.shape)


                # Average signal hist - histNat over historical period, rcp85 - mean(histNat) over projection period
                if domain['Atlantic'] != None:
                    varsignal_a[0:145,k,j] = averageDom(varh_a-varhn_a, 3, domain['Atlantic'], lat, density)
                    varsignal_a[145:,k,j] = averageDom(varrcp_a-meanvarhn_a, 3, domain['Atlantic'], lat, density)
                if domain['Pacific'] != None:
                    varsignal_p[0:145,k,j] = averageDom(varh_p-varhn_p, 3, domain['Pacific'], lat, density)
                    varsignal_p[145:,k,j] = averageDom(varrcp_p-meanvarhn_p, 3, domain['Pacific'], lat, density)
                if domain['Indian'] != None:
                    varsignal_i[0:145,k,j] = averageDom(varh_i-varhn_i, 3, domain['Indian'], lat, density)
                    varsignal_i[145:,k,j] = averageDom(varrcp_i-meanvarhn_i, 3, domain['Indian'], lat, density)

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
        fileName = 'cmip5.'+model['name']+'.toe_rcp_histNat_method2_'+method_noise+'.nc'
        dir = '/home/ysilvy/Density_bining/Yona_analysis/data/toe_rcp85_histNat_average_signal/'+method_noise+'/'
        fout = open_ncfile(dir+fileName,'w', format='NETCDF4')
        if method_noise == 'average_std':
            noise_description = 'Noise is computed by averaging the standard deviation of the historical Nat ensemble mean runs ' \
                                'in the specified domains, over the historical period only.'
        else :
            noise_description = 'Noise is computed by averaging the historical Nat ensemble mean runs in the specified domains ' \
                                'over the historical period only, and then taking the standard deviation of the average.'
        fout.description = 'ToE hist+rcp8.5 vs. histNat for each member, in 5 domains : Southern Subtropics (0), Southern Ocean (1),' \
                            'Northern Subtropics (2), North Atlantic (3), North Pacific (4). \n' \
                           'The historical runs are prolonged by the 95 years of RCP8.5. \n Signal is computed by averaging ' \
                           'the difference (historical+RCP8.5) - timeaverage(historicalNat) in those domains. The ensemble mean historicalNat is ' \
                           'used here for all runs of each model. ' + noise_description + ' Then ToE is computed ' \
                            'using twice the noise as the limit.'

        # dimensions
        fout.createDimension('members', nruns)
        fout.createDimension('basin', 4)
        fout.createDimension('domain', 5)

        # variables
        members = fout.createVariable('members', 'f4', ('members',))
        basin = fout.createVariable('basin', 'i4', ('basin',))
        domain = fout.createVariable('domain', 'f4', ('domain',))
        varToE2 = fout.createVariable(varname['var_zonal_w/bowl']+'ToE2', 'f4', ('members','basin','domain',))

        # data
        members[:] =  np.arange(0,nruns)
        basin[:] =  np.arange(0,basinN)
        domain[:] = np.arange(0,len(domains))
        varToE2[:,:,:] = varToE

        # units
        basin.units = 'basin index'
        domain.units = 'domain index'
        varToE2.units = 'Year'
        varToE2.long_name = 'ToE 2 for ' + legVar

        fout.close()



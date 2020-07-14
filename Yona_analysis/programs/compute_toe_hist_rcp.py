#!/bin/env python
# -*- coding: utf-8 -*-

"""
Average signal and noise for each run of each model in the specified domains, then compute ToE hist + RCP85
Over the 2005-2100 period, use either RCP85-average(histNat) or RCP85-average(PiControl) as the signal
Save ToE in output files

Dec 2018 : change use_piC criterion to signal = hist+RCP8.5 - mean(PiCtrl) and noise = mult*std(PiCtrl) over the whole period
(historical + projection)
"""

import numpy as np
from netCDF4 import Dataset as open_ncfile
from maps_matplot_lib import defVarmme, averageDom
from modelsDef import defModels
from libToE import findToE, findToE_2thresholds, ToEdomainrcp85vshistNat
import glob, os

# ----- Workspace ------

indir_histrcp85 = '/home/ysilvy/Density_bining/Yona_analysis/data/hist_rcp85/'
indir_histNat = '/data/ericglod/Density_binning/Prod_density_april15/mme_histNat/'
indir_piC = '/data/ericglod/Density_binning/Prod_density_april15/mme_piControl/'
indir_noise = '/home/ysilvy/Density_bining/Yona_analysis/data/noise_estimate/RCP85vshistNat_domains/'

models = defModels()

# ----- Work ------

varname = defVarmme('salinity'); v = 'S'
#varname = defVarmme('depth'); v = 'Z'

method = 'average_signal' # Average signal and noise in the box, then compute ToE

# -- Choose which 'noise' to use for the ToE calculation
# method_noise = 'average_std' # Average the standard deviation of histNat or PiC in the specified domains
method_noise = 'average_histNat' # Average histNat (or piC) in the specified domains then determine the std
# of this averaged value

domains = ['Southern ST', 'SO', 'Northern ST', 'North Atlantic', 'North Pacific']
signal_domains = ['fresher','saltier','fresher','saltier','saltier']

multStd = 2. # detect ToE at multStd std dev of histNat (or PiControl)
# In fact now we save both 1 std and multStd=2 std

# use_piC = False # Over projection period, signal = RCP-average(histNat), noise = std(histNat)
use_piC = True # Over projection period, signal = RCP-average(PiControl), noise = std(PiControl)

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

nMembers = np.ma.zeros(len(models)) # Initialize array for keeping nb of members per model

for i, model in enumerate(models):

    print('Working on', model['name'])

    if (use_piC != True) or (use_piC == True and model['name'] != 'GISS-E2-R' and model['name'] != 'FGOALS-g2'
                             and model['name'] != 'MIROC-ESM'):

        # Read hist+rcp85 files
        listruns = sorted(glob.glob(indir_histrcp85 + 'cmip5.' + model['name'] + '.' + '*zon2D.nc'))
        nruns = len(listruns)
        nMembers[i] = nruns
        if nruns != 0:
            # Index of common time interval
            tstart = model['props'][2]
            tend = model['props'][3]

            if use_piC == False :
                # Read histNat ensemble mean
                filehn = 'cmip5.' + model['name'] + '.historicalNat.ensm.an.ocn.Omon.density.ver-' \
                         + model['file_end_histNat'] + '_zon2D.nc'
                fhn = open_ncfile(indir_histNat + filehn,'r')
                # Read var histNat
                varhn_a = fhn.variables[var][:,1,:,:].squeeze()
                varhn_p = fhn.variables[var][:,2,:,:].squeeze()
                varhn_i = fhn.variables[var][:,3,:,:].squeeze()
                #print(' varhn shape : ', varhn_a.shape)
                # Read std of histNat (max std of all runs for each model)
                if method_noise == 'average_std':
                    varstd_a = fhn.variables[var+'Std'][1,:,:].squeeze()
                    varstd_p = fhn.variables[var+'Std'][2,:,:].squeeze()
                    varstd_i = fhn.variables[var+'Std'][3,:,:].squeeze()

                # Compute time average of the whole histNat series
                meanvarhn_a = np.ma.mean(varhn_a, axis=0)
                meanvarhn_p = np.ma.mean(varhn_p, axis=0)
                meanvarhn_i = np.ma.mean(varhn_i, axis=0)

            if use_piC == True:
                # Read and Compute time average of PiControl over 240 years + std of PiControl
                filepiC = glob.glob(indir_piC + 'cmip5.' + model['name'] + '.' + '*zon2D.nc')[0]
                fpiC = open_ncfile(filepiC,'r')
                varpiC_a = fpiC.variables[var][-240:,1,:,:].squeeze()
                varpiC_p = fpiC.variables[var][-240:,2,:,:].squeeze()
                varpiC_i = fpiC.variables[var][-240:,3,:,:].squeeze()
                meanvarpiC_a = np.ma.mean(varpiC_a, axis=0)
                meanvarpiC_p = np.ma.mean(varpiC_p, axis=0)
                meanvarpiC_i = np.ma.mean(varpiC_i, axis=0)
                if method_noise == 'average_std':
                    stdvarpiC_a = np.ma.std(varpiC_a, axis=0)
                    stdvarpiC_p = np.ma.std(varpiC_p, axis=0)
                    stdvarpiC_i = np.ma.std(varpiC_i, axis=0)

            # Initialize varnoise for each basin, containing averaged noise for each domain
            varnoise_a = np.ma.masked_all(len(domains))
            varnoise_p = np.ma.masked_all(len(domains))
            varnoise_i = np.ma.masked_all(len(domains))
            # Initialize varsignal for each basin, containing averaged signal for each domain and each run
            varsignal_a = np.ma.masked_all((timN,nruns,len(domains)))
            varsignal_p = np.ma.masked_all((timN,nruns,len(domains)))
            varsignal_i = np.ma.masked_all((timN,nruns,len(domains)))
            # Initialize toe 1 and 2 for each basin (run, domain)
            toe1_a = np.ma.masked_all((nruns,len(domains)))
            toe1_p = np.ma.masked_all((nruns,len(domains)))
            toe1_i = np.ma.masked_all((nruns,len(domains)))
            toe2_a = np.ma.masked_all((nruns,len(domains)))
            toe2_p = np.ma.masked_all((nruns,len(domains)))
            toe2_i = np.ma.masked_all((nruns,len(domains)))
            # Initialize output variable
            varToE1 = np.ma.masked_all((nruns,basinN,len(domains))) # (members,basin,domain) 1std
            varToE2 = np.ma.masked_all((nruns,basinN,len(domains))) # (members,basin,domain) 2std

            if method_noise == 'average_histNat':
                filenoise = 'cmip5.' + model['name'] + '.noise_domains_hist_histNat.std_of_average.nc'
                fnoise = open_ncfile(indir_noise + filenoise,'r')
                if use_piC == False:
                    # Here we read the std of the averaged histNat in the 5 domains for all runs, then take the max as our noise
                    varstd = fnoise.variables[var+'stdhn'][:]
                    varnoise_a = np.ma.max(varstd[:,1,:], axis=0)
                    varnoise_p = np.ma.max(varstd[:,2,:], axis=0)
                    varnoise_i = np.ma.max(varstd[:,3,:], axis=0)
                else: 
                    # Reading std of averaged PiControl in the 5 domains (boxes coordinates corresponding to RCP8.5 - histNat)
                    varstd = fnoise.variables[var+'stdpiC'][:]
                    varnoise_a = varstd[1,:]
                    varnoise_p = varstd[2,:]
                    varnoise_i = varstd[3,:]

            # Loop over 5 domains
            for j, domain_name in enumerate(domains):
                print('- ', j, domains[j])

                # Select domain to average
                domain = ToEdomainrcp85vshistNat(model['name'], domain_name)[0]
                domain_char = ToEdomainrcp85vshistNat(model['name'], domain_name)[1]


                if use_piC == False :
                    if method_noise == 'average_std':
                        # Average noise when using histNat
                        if domain['Atlantic'] != None:
                            varnoise_a[j] = averageDom(varstd_a, 2, domain['Atlantic'], lat, density)
                        if domain['Pacific'] != None:
                            varnoise_p[j] = averageDom(varstd_p, 2, domain['Pacific'], lat, density)
                        if domain['Indian'] != None:
                            varnoise_i[j] = averageDom(varstd_i, 2, domain['Indian'], lat, density)
                else :
                # Average noise when using PiControl
                    if domain['Atlantic'] != None:
                        if method_noise == 'average_std':
                            varnoise_a[j] = averageDom(stdvarpiC_a, 2, domain['Atlantic'], lat, density)
#                         else :
#                             varnoise_a[j] = np.ma.std(averageDom(varpiC_a, 3, domain['Atlantic'], lat, density), axis=0)
                    if domain['Pacific'] != None:
                        if method_noise == 'average_std':
                            varnoise_p[j] = averageDom(stdvarpiC_p, 2, domain['Pacific'], lat, density)
#                         else:
#                             varnoise_p[j] = np.ma.std(averageDom(varpiC_p, 3, domain['Pacific'], lat, density), axis=0)
                    if domain['Indian'] != None:
                        if method_noise == 'average_std':
                            varnoise_i[j] = averageDom(stdvarpiC_i, 2, domain['Indian'], lat, density)
#                         else:
#                             varnoise_i[j] = np.ma.std(averageDom(varpiC_i, 3, domain['Indian'], lat, density), axis=0)

                if j==0:
                    # Initialize list for run labels
                    run_names = ['']*nruns
                
                # Loop over number of runs
                for k in range(nruns):
            
                    namefile = os.path.basename(listruns[k])
                    run_nb = namefile.split('.')[3]
                    if j==0:
                        run_names[k] = run_nb # Save run label, i.e. 'r1i1p1'
                    print('    . run number', k, run_nb)

                    # Read file
                    fhrcp = open_ncfile(listruns[k],'r')
                    # Read var hist + rcp85
                    varhrcp_a = fhrcp.variables[var][tstart:tend+95,1,:,:].squeeze()
                    varhrcp_p = fhrcp.variables[var][tstart:tend+95,2,:,:].squeeze()
                    varhrcp_i = fhrcp.variables[var][tstart:tend+95,3,:,:].squeeze()

                    # Average signal hist+RCP8.5 - histNat or hist+RCP8.5 - PiControl
                    if use_piC == False:
                        if domain['Atlantic'] != None:
                            varsignal_a[:,k,j] = averageDom(varhrcp_a-meanvarhn_a, 3, domain['Atlantic'], lat, density)
                            #varsignal_a[145:,k,j] = averageDom(varhrcp_a[145:,:,:]-meanvarhn_a, 3, domain['Atlantic'], lat, density)
                        if domain['Pacific'] !=None:
                            varsignal_p[:,k,j] = averageDom(varhrcp_p-meanvarhn_p, 3, domain['Pacific'], lat, density)
                            #varsignal_p[145:,k,j] = averageDom(varhrcp_p[145:,:,:]-meanvarhn_p, 3, domain['Pacific'], lat, density)
                        if domain['Indian'] != None:
                            varsignal_i[:,k,j] = averageDom(varhrcp_i-meanvarhn_i, 3, domain['Indian'], lat, density)
                            #varsignal_i[145:,k,j] = averageDom(varhrcp_i[145:,:,:]-meanvarhn_i, 3, domain['Indian'], lat, density)
                    else:
                        if domain['Atlantic'] != None:
                            varsignal_a[:,k,j] = averageDom(varhrcp_a-meanvarpiC_a, 3, domain['Atlantic'], lat, density)
                        if domain['Pacific'] !=None:
                            varsignal_p[:,k,j] = averageDom(varhrcp_p-meanvarpiC_p, 3, domain['Pacific'], lat, density)
                        if domain['Indian'] != None:
                            varsignal_i[:,k,j] = averageDom(varhrcp_i-meanvarpiC_i, 3, domain['Indian'], lat, density)

                    # print('      varsignal shape:', varsignal_a.shape, varsignal_p.shape, varsignal_i.shape)

                    # Compute ToE of averaged domain for run k
                    if domain['Atlantic'] != None and np.ma.is_masked(varnoise_a[j]) == False:
                        toe2_a[k,j] = findToE(varsignal_a[:,k,j], varnoise_a[j], multStd) + iniyear
                        toe1_a[k,j] = findToE(varsignal_a[:,k,j], varnoise_a[j], 1) + iniyear
#                         print(toe1_a[k,j], toe2_a[k,j])
                    if domain['Pacific'] != None and np.ma.is_masked(varnoise_p[j]) == False:
                        toe2_p[k,j] = findToE(varsignal_p[:,k,j], varnoise_p[j], multStd) + iniyear
                        toe1_p[k,j] = findToE(varsignal_p[:,k,j], varnoise_p[j], 1) + iniyear
                    if domain['Indian'] != None and np.ma.is_masked(varnoise_i[j]) == False:
                        toe2_i[k,j] = findToE(varsignal_i[:,k,j], varnoise_i[j], multStd) + iniyear
                        toe1_i[k,j] = findToE(varsignal_i[:,k,j], varnoise_i[j], 1) + iniyear

                    # Take out runs where the signal is of opposite sign than expected
                    if v != 'Z':
                        if signal_domains[j] == 'fresher':
                            if np.ma.mean(varsignal_a[-5:,k,j],axis=0) > 2 * varnoise_a[j]:
                                toe2_a[k,j] = np.ma.masked
                            if np.ma.mean(varsignal_a[-5:,k,j],axis=0) > 1 * varnoise_a[j]:
                                toe1_a[k,j] = np.ma.masked
                            if np.ma.mean(varsignal_p[-5:,k,j],axis=0) > 2 * varnoise_p[j]:
                                toe2_p[k,j] = np.ma.masked
                            if np.ma.mean(varsignal_p[-5:,k,j],axis=0) > 1 * varnoise_p[j]:
                                toe1_p[k,j] = np.ma.masked
                            if np.ma.mean(varsignal_i[-5:,k,j],axis=0) > 2 * varnoise_i[j]:
                                toe2_i[k,j] = np.ma.masked
                            if np.ma.mean(varsignal_i[-5:,k,j],axis=0) > 1 * varnoise_i[j]:
                                toe1_i[k,j] = np.ma.masked
                        else:
                            if np.ma.mean(varsignal_a[-5:,k,j],axis=0) < -2 * varnoise_a[j]:
                                toe2_a[k,j] = np.ma.masked
                            if np.ma.mean(varsignal_a[-5:,k,j],axis=0) < -1 * varnoise_a[j]:
                                toe1_a[k,j] = np.ma.masked
                            if np.ma.mean(varsignal_p[-5:,k,j],axis=0) < -2 * varnoise_p[j]:
                                toe2_p[k,j] = np.ma.masked
                            if np.ma.mean(varsignal_p[-5:,k,j],axis=0) < -1 * varnoise_p[j]:
                                toe1_p[k,j] = np.ma.masked
                            if np.ma.mean(varsignal_i[-5:,k,j],axis=0) < -2 * varnoise_i[j]:
                                toe2_i[k,j] = np.ma.masked
                            if np.ma.mean(varsignal_i[-5:,k,j],axis=0) < -1 * varnoise_i[j]:
                                toe1_i[k,j] = np.ma.masked                 

            varToE1[:,1,:] = toe1_a
            varToE1[:,2,:] = toe1_p
            varToE1[:,3,:] = toe1_i
            varToE2[:,1,:] = toe2_a
            varToE2[:,2,:] = toe2_p
            varToE2[:,3,:] = toe2_i
            print('  varToE shape:', varToE2.shape)
            print('  ', np.ma.around(np.ma.median(varToE1[:,1,0])), np.ma.around(np.ma.median(varToE2[:,1,0])))
            print('')

            # Save in output file
            if use_piC != True:
                fileName = 'cmip5.'+model['name']+'.'+legVar+'_toe_rcp_histNat_'+method_noise+'.nc'
                dir = '/home/ysilvy/Density_bining/Yona_analysis/data/toe_rcp85_histNat_average_signal/'+method_noise+'/'
            else :
                if method_noise == 'average_histNat':
                    method = 'average_piC'
                else:
                    method = method_noise
                fileName = 'cmip5.'+model['name']+'.'+legVar+'_toe_rcp_PiControl_method2_'+method+'.nc'
                dir = '/home/ysilvy/Density_bining/Yona_analysis/data/toe_rcp85_PiControl_average_signal/'+method+'/'
            fout = open_ncfile(dir+fileName,'w', format='NETCDF4')
            if method_noise == 'average_std':
                if use_piC != True:
                    noise_description = 'Noise is computed by averaging the standard deviation of the historical Nat runs ' \
                                    'in the specified domains, then taking the max std, over the historical period only.'
                else :
                    noise_description = 'Noise is computed by averaging the standard deviation of the Pre-industrial Control ensemble mean.'
            else :
                if use_piC != True:
                    noise_description = 'Noise is computed by averaging the historical Nat runs in the specified domains, ' \
                                    'and then taking the standard deviation of the average,' \
                                        'and keeping the max std of all runs.'
                else:
                    noise_description = 'Noise is computed by averaging the Pre-industrial control ensemble mean in the specified domains ' \
                                    'over the last 240 years, and then taking the standard deviation of the average.'

            if use_piC != True:
                fout.description = 'ToE hist+rcp8.5 vs. histNat for each member, in 5 domains : Southern Subtropics (0), Southern Ocean (1),' \
                                'Northern Subtropics (2), North Atlantic (3), North Pacific (4). \n' \
                               'The historical runs are prolonged by the 95 years of RCP8.5. \n Signal is computed by averaging ' \
                               'the difference (historical+RCP8.5)-timeaverage(historicalNat) in those domains. ' \
                               'The ensemble mean historicalNat is used here for all historical runs of the model. ' \
                                   + noise_description + ' Then ToE is computed ' \
                                'using once or twice the noise as the limit.'
            else :
                fout.description = 'ToE hist + rcp8.5 vs. PiControl for each member, in 5 domains : Southern Subtropics (0), Southern Ocean (1),' \
                                'Northern Subtropics (2), North Atlantic (3), North Pacific (4). \n' \
                               'The historical runs are prolonged by the 95 years of RCP8.5. \n Signal is computed by averaging ' \
                               'the difference (historical + RCP8.5) - timeaverage(PiControl_last240years) in those domains. ' \
                                   'The ensemble mean PiControl is ' \
                               'used here for all historical runs of the model. ' + noise_description + ' Then ToE is computed ' \
                                'using once or twice the noise as the limit.'

            # dimensions
            fout.createDimension('members', nruns)
            fout.createDimension('basin', 4)
            fout.createDimension('domain', 5)

            # variables
            members = fout.createVariable('members', 'f4', ('members',))
            basin = fout.createVariable('basin', 'i4', ('basin',))
            domain = fout.createVariable('domain', 'f4', ('domain',))
            ToE1 = fout.createVariable(varname['var_zonal_w/bowl']+'ToE1', 'f4', ('members','basin','domain',))
            ToE2 = fout.createVariable(varname['var_zonal_w/bowl']+'ToE2', 'f4', ('members','basin','domain',))
            run_label = fout.createVariable('run_label', 'S6', ('members',))

            # data
            members[:] =  np.arange(0,nruns)
            basin[:] =  np.arange(0,basinN)
            domain[:] = np.arange(0,len(domains))
            ToE1[:,:,:] = varToE1
            ToE2[:,:,:] = varToE2
            for k in range(nruns):
                run_label[k] = run_names[k]

            # units
            basin.units = 'basin index'
            domain.units = 'domain index'
            ToE2.units = 'Year'
            ToE2.long_name = 'ToE (>2std) for ' + legVar
            ToE1.units = 'Year'
            ToE1.long_name = 'ToE (>1std) for ' + legVar
            run_label.long_name = 'Run number'

            fout.close()



#!/bin/env python
# -*- coding: utf-8 -*-


"""
Compute the noise estimate in the different domains for historical, historicalNat and PiControl simulations
and save in output file
"""


import numpy as np
from netCDF4 import Dataset as open_ncfile
from modelsDef import defModels, defModelsCO2piC
from maps_matplot_lib import defVarmme, averageDom
from libToE import ToEdomainhistvshistNat, ToEdomain1pctCO2vsPiC, ToEdomainrcp85vshistNat
import glob


# ----- Workspace ------

indir_hist = '/data/ericglod/Density_binning/Prod_density_april15/historical/'
indir_histNat = '/data/ericglod/Density_binning/Prod_density_april15/historicalNat/'
indir_piC = '/data/ericglod/Density_binning/Prod_density_april15/mme_piControl/'

models = defModels()
modelspiC = defModelsCO2piC()

# method_noise = 'average_std' # Compute std of the time series then average in the domains
method_noise = 'std_of_average' # Average time series in the domains then compute std

if method_noise == 'average_std':
    noise_description = 'The standard deviation is computed and then averaged in the domains.'
else:
    noise_description = 'The variables are averaged in the domains, then the standard deviation is computed ' \
                        'over those domains.'


# ----- Work ------

# varname = defVarmme('salinity'); v = 'S'
# varname = defVarmme('temp'); v = 'T'
varname = defVarmme('depth'); v = 'Z'

iniyear = 1860
finalyear = 2005


# ----- Variables ------

domains = ['Southern ST', 'SO', 'Northern ST', 'North Atlantic', 'North Pacific']

# Choose random file to read only the basic variables and properties common to all files
file = 'cmip5.' + models[0]['name'] + '.historicalNat.r1i1p1.an.ocn.Omon.density.ver-' + models[0]['file_end_histNat'] + '_zon2D.nc'
f = open_ncfile(indir_histNat + file,'r')

lat = f.variables['latitude'][:]; latN = lat.size
density = f.variables['lev'][:]; levN = density.size
time = f.variables['time'][:]; timN = time.size
basinN = 4
var = varname['var_zonal_w/bowl']
legVar = varname['legVar']
unit = varname['unit']


# ----- Compute noise and average in domains for each run, then save in output file for each model -----

# == Historical and historicalNat + PiControl in RCP8.5vshistNat boxes (for RCP8.5 ToE calculation) ==

nMembers_h = np.ma.empty(len(models)) # Initialize array for keeping nb of members per model
nMembers_hn = np.ma.empty(len(models))

for i, model in enumerate(models):

    print('Working on', model['name'])

    if model['hist-rcp85'] != [] :

        # Index of common time interval
        tstart = model['props'][2]
        tend = model['props'][3]

        # Read histNat files
        listruns_hn = glob.glob(indir_histNat + 'cmip5.'+model['name']+'*zon2D.nc')
        nruns_hn = len(listruns_hn)
        nMembers_hn[i] = nruns_hn

        # Read hist files
        listruns_h = glob.glob(indir_hist + 'cmip5.'+model['name']+'*zon2D.nc')
        nruns_h = len(listruns_h)
        nMembers_h[i] = nruns_h

        # Initialize varnoise for each basin, containing averaged noise for each domain and each run of the model
        varnoise_hna = np.ma.masked_all((nruns_hn,len(domains)))
        varnoise_hnp = np.ma.masked_all((nruns_hn,len(domains)))
        varnoise_hni = np.ma.masked_all((nruns_hn,len(domains)))
        varnoise_ha = np.ma.masked_all((nruns_h,len(domains)))
        varnoise_hp = np.ma.masked_all((nruns_h,len(domains)))
        varnoise_hi = np.ma.masked_all((nruns_h,len(domains)))
        # Initialize output variables
        varnoise_hn = np.ma.masked_all((nruns_hn,basinN,len(domains)))
        varnoise_h = np.ma.masked_all((nruns_h,basinN,len(domains)))

        # Initialize varnoise PiControl for each basin, containing averaged noise for each domain
        varnoise_piCa = np.ma.masked_all(len(domains))
        varnoise_piCp = np.ma.masked_all(len(domains))
        varnoise_piCi = np.ma.masked_all(len(domains))
        # Initialize output variable
        varnoise_piC = np.ma.masked_all((basinN,len(domains)))

        # Read PiC file
        filepiC = glob.glob(indir_piC + 'cmip5.' + model['name'] + '.' + '*zon2D.nc')
        if len(filepiC) !=0 :
            fpiC = open_ncfile(filepiC[0],'r')
            # Read var PiC over last 240 years
            varpiC = fpiC.variables[var][-240:,:,:,:]
            # Compute std of PiC
            varstdpiC = np.ma.std(varpiC, axis=0)

        for j, domain_name in enumerate(domains):

            print('- ', domain_name)

            # Select domain to average
            domain = ToEdomainrcp85vshistNat(model['name'], domain_name)[0]
            domain_char = ToEdomainrcp85vshistNat (model['name'], domain_name)[1]

            # Loop over number of histNat runs
            for khn in range(nruns_hn):
                print('      . histNat run number', khn)

                # Read file
                fhn = open_ncfile(listruns_hn[khn],'r')
                # Read var histNat
                varhn = fhn.variables[var][tstart:tend,:,:,:]
                # Compute std of histNat
                varstdhn = np.ma.std(varhn, axis=0)
                # Average std
                if domain['Atlantic'] != None:
                    if method_noise == 'average_std':
                        varnoise_hna[khn,j] = averageDom(varstdhn[1,:,:], 2, domain['Atlantic'], lat, density)
                    else:
                        varnoise_hna[khn,j] = np.ma.std(averageDom(varhn[:,1,:,:], 3, domain['Atlantic'], lat, density), axis=0)
                if domain['Pacific'] != None:
                    if method_noise == 'average_std':
                        varnoise_hnp[khn,j] = averageDom(varstdhn[2,:,:], 2, domain['Pacific'], lat, density)
                    else:
                        varnoise_hnp[khn,j] = np.ma.std(averageDom(varhn[:,2,:,:], 3, domain['Pacific'], lat, density), axis=0)
                if domain['Indian'] != None:
                    if method_noise == 'average_std':
                        varnoise_hni[khn,j] = averageDom(varstdhn[3,:,:], 2, domain['Indian'], lat, density)
                    else:
                        varnoise_hni[khn,j] = np.ma.std(averageDom(varhn[:,3,:,:], 3, domain['Indian'], lat, density), axis=0)

            # Loop over number of hist runs
            for kh in range(nruns_h):
                print('      . hist run number', kh)

                # Read file
                fh = open_ncfile(listruns_h[kh],'r')
                # Read var hist
                varh = fh.variables[var][tstart:tend,:,:,:]
                # Compute std of hist
                varstdh = np.ma.std(varh, axis=0)
                # Average
                if domain['Atlantic'] != None:
                    if method_noise == 'average_std':
                        varnoise_ha[kh,j] = averageDom(varstdh[1,:,:], 2, domain['Atlantic'], lat, density)
                    else:
                        varnoise_ha[kh,j] = np.ma.std(averageDom(varh[:,1,:,:], 3, domain['Atlantic'], lat, density), axis=0)
                if domain['Pacific'] != None:
                    if method_noise == 'average_std':
                        varnoise_hp[kh,j] = averageDom(varstdh[2,:,:], 2, domain['Pacific'], lat, density)
                    else:
                        varnoise_hp[kh,j] = np.ma.std(averageDom(varh[:,2,:,:], 2, domain['Pacific'], lat, density), axis=0)
                if domain['Indian'] != None:
                    if method_noise == 'average_std':
                        varnoise_hi[kh,j] = averageDom(varstdh[3,:,:], 2, domain['Indian'], lat, density)
                    else:
                        varnoise_hi[kh,j] = np.ma.std(averageDom(varh[:,3,:,:], 3, domain['Indian'], lat, density), axis=0)

            # Average PiControl
            if len(filepiC)!=0:
                print('      . PiControl')
                if domain['Atlantic'] != None:
                    if method_noise == 'average_std':
                        varnoise_piCa[j] = averageDom(varstdpiC[1,:,:], 2, domain['Atlantic'], lat, density)
                    else :
                        varnoise_piCa[j] = np.ma.std(averageDom(varpiC[:,1,:,:], 3, domain['Atlantic'], lat, density), axis=0)
                if domain['Pacific'] != None:
                    if method_noise == 'average_std':
                        varnoise_piCp[j] = averageDom(varstdpiC[2,:,:], 2, domain['Pacific'], lat, density)
                    else:
                        varnoise_piCp[j] = np.ma.std(averageDom(varpiC[:,2,:,:], 3, domain['Pacific'], lat, density), axis=0)
                if domain['Indian'] != None:
                    if method_noise == 'average_std':
                        varnoise_piCi[j] = averageDom(varstdpiC[3,:,:], 2, domain['Indian'], lat, density)
                    else:
                        varnoise_piCi[j] = np.ma.std(averageDom(varpiC[:,3,:,:], 3, domain['Indian'], lat, density), axis=0)

        varnoise_hn[:,1,:] = varnoise_hna
        varnoise_hn[:,2,:] = varnoise_hnp
        varnoise_hn[:,3,:] = varnoise_hni
        varnoise_h[:,1,:] = varnoise_ha
        varnoise_h[:,2,:] = varnoise_hp
        varnoise_h[:,3,:] = varnoise_hi
        varnoise_piC[1,:] = varnoise_piCa
        varnoise_piC[2,:] = varnoise_piCp
        varnoise_piC[3,:] = varnoise_piCi
        print('  varnoise_hn shape:', varnoise_hn.shape)
        print('  varnoise_h shape:', varnoise_h.shape)
        print('  varnoise_piC shape:', varnoise_piC.shape)

                
        # Save in output file
        fileName = 'cmip5.'+model['name']+'.'+legVar+'_noise_domains_hist_histNat.' + method_noise + '.nc'
        dir = '/home/ysilvy/Density_bining/Yona_analysis/data/noise_estimate/RCP85vshistNat_domains/'
        fout = open_ncfile(dir+fileName,'w', format='NETCDF4')
        fout.description = 'Standard deviation of historical, historicalNat and PiControl for each member, in 5 domains : ' \
                           'Southern Subtropics (0), Southern Ocean (1), Northern Subtropics (2), North Atlantic (3), ' \
                           'North Pacific (4) (taken in RCP8.5 vs. histNat boxes).' + noise_description

        # dimensions
        fout.createDimension('members_hist', nruns_h)
        fout.createDimension('members_histNat', nruns_hn)
        fout.createDimension('basin', 4)
        fout.createDimension('domain', 5)

        # variables
        members_hist = fout.createVariable('members_hist', 'f4', ('members_hist',))
        members_histNat = fout.createVariable('members_histNat', 'f4', ('members_histNat',))
        basin = fout.createVariable('basin', 'i4', ('basin',))
        domain = fout.createVariable('domain', 'f4', ('domain',))
        varstdh = fout.createVariable(varname['var_zonal_w/bowl']+'stdh', 'f4', ('members_hist','basin','domain',))
        varstdhn = fout.createVariable(varname['var_zonal_w/bowl']+'stdhn', 'f4', ('members_histNat','basin','domain',))
        varstdpiC = fout.createVariable(varname['var_zonal_w/bowl']+'stdpiC', 'f4', ('basin','domain',))

        # data
        members_hist[:] =  np.arange(0,nruns_h)
        members_histNat[:] =  np.arange(0,nruns_hn)
        basin[:] =  np.arange(0,basinN)
        domain[:] = np.arange(0,len(domains))
        varstdh[:,:,:] = varnoise_h
        varstdhn[:,:,:] = varnoise_hn
        varstdpiC[:,:] = varnoise_piC

        # units
        basin.units = 'basin index'
        domain.units = 'domain index'
        varstdh.units = unit
        varstdhn.units = unit
        varstdpiC.units = unit
        varstdh.long_name = 'Std for ' + legVar + ' of historical runs'
        varstdhn.long_name = 'Std for ' + legVar + ' of historicalNat runs'
        varstdpiC.long_name = 'Std for ' + legVar + ' of pre-industrial control run'

        fout.close()



# # == Pi Control ==
#
# for i, model in enumerate(modelspiC):
#
#     print('Working on', model['name'])
#
#     # Read PiC file
#     file_piC = 'cmip5.' + model['name'] + '.piControl.ensm.an.ocn.Omon.density.ver-' + model['file_end_piC'] + '_zon2D.nc'
#     fpiC = open_ncfile(indir_piC + file_piC,'r')
#
#     # Read var PiC
#     varpiC = fpiC.variables[var][-140:,:,:,:]
#     # Compute std of PiC
#     varstdpiC = np.ma.std(varpiC, axis=0)
#
#     # Initialize varnoise for each basin, containing averaged noise for each domain
#     varnoise_piCa = np.ma.masked_all(len(domains))
#     varnoise_piCp = np.ma.masked_all(len(domains))
#     varnoise_piCi = np.ma.masked_all(len(domains))
#     # Initialize output variable
#     varnoise_piC = np.ma.masked_all((basinN,len(domains)))
#
#     for j, domain_name in enumerate(domains):
#
#         print('- ', domain_name)
#
#         # Select domain to average
#         domain = ToEdomain1pctCO2vsPiC(model['name'], domain_name)[0]
#         domain_char = ToEdomain1pctCO2vsPiC(model['name'], domain_name)[1]
#
#         # Average
#         if domain['Atlantic'] != None:
#             if method_noise == 'average_std':
#                 varnoise_piCa[j] = averageDom(varstdpiC[1,:,:], 2, domain['Atlantic'], lat, density)
#             else :
#                 varnoise_piCa[j] = np.ma.std(averageDom(varpiC[:,1,:,:], 3, domain['Atlantic'], lat, density), axis=0)
#         if domain['Pacific'] != None:
#             if method_noise == 'average_std':
#                 varnoise_piCp[j] = averageDom(varstdpiC[2,:,:], 2, domain['Pacific'], lat, density)
#             else:
#                 varnoise_piCp[j] = np.ma.std(averageDom(varpiC[:,2,:,:], 3, domain['Pacific'], lat, density), axis=0)
#         if domain['Indian'] != None:
#             if method_noise == 'average_std':
#                 varnoise_piCi[j] = averageDom(varstdpiC[3,:,:], 2, domain['Indian'], lat, density)
#             else:
#                 varnoise_piCi[j] = np.ma.std(averageDom(varpiC[:,3,:,:], 3, domain['Indian'], lat, density), axis=0)
#
#
#     varnoise_piC[1,:] = varnoise_piCa
#     varnoise_piC[2,:] = varnoise_piCp
#     varnoise_piC[3,:] = varnoise_piCi
#     print('  varnoise_piC shape:', varnoise_piC.shape)
#
#
#     # Save in output file
#     fileName = 'cmip5.'+model['name']+'.noise_domains_PiControl.' + method_noise + '.nc'
#     dir = '/home/ysilvy/Density_bining/Yona_analysis/data/noise_estimate/'
#     fout = open_ncfile(dir+fileName,'w', format='NETCDF4')
#     fout.description = 'Standard deviation of pre-industrial control run in 5 domains : Southern Subtropics (0), ' \
#                        'Southern Ocean (1), Northern Subtropics (2), North Atlantic (3), North Pacific (4) .' + noise_description
#
#     # dimensions
#     fout.createDimension('basin', 4)
#     fout.createDimension('domain', 5)
#
#     # variables
#     basin = fout.createVariable('basin', 'i4', ('basin',))
#     domain = fout.createVariable('domain', 'f4', ('domain',))
#     varstdpiC = fout.createVariable(varname['var_zonal_w/bowl']+'stdpiC', 'f4', ('basin','domain',))
#
#     # data
#     basin[:] =  np.arange(0,basinN)
#     domain[:] = np.arange(0,len(domains))
#     varstdpiC[:,:] = varnoise_piC
#
#     # units
#     basin.units = 'basin index'
#     domain.units = 'domain index'
#     varstdpiC.units = unit
#     varstdpiC.long_name = 'Std for ' + legVar + ' of pre-industrial control run'
#
#     fout.close()

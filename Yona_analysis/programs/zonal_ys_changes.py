#!/bin/env python
# -*- coding: utf-8 -*-

"""
Python matplotlib
Plot density/latitude salinity changes for each basin (1950-2000)
Choose between :
- Durack and Wijffels data
- Historical simulations  (mme or ensemble means)
- hist vs. histNat (mme or ensemble means)
- 1%CO2 (mme or single models)
- 1%CO2 vs. piControl (mme or single models)
- RCP8.5 vs. histNat

"""

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset as open_ncfile
from maps_matplot_lib import defVarDurack, zonal_2D, defVarmme, custom_div_cmap, modelagree
from modelsDef import defModels, defModelsCO2piC
import glob
import datetime

# ----- Workspace ------

# name = 'Durack & Wijffels'
#name = 'mme_hist'
# name = 'mme_hist_histNat'
# name = 'ens_mean_hist'
#name = 'ens_mean_hist_histNat'
# name = 'mme_1pctCO2vsPiC'
# name = 'mme_1pctCO2'
# name = '1pctCO2'
# name = '1pctCO2vsPiC'
name = 'mme_rcp85_histNat'
# name = 'ens_mean_rcp85_histNat'

# -- Choose where to stop for 1%CO2 simulations : 2*CO2 (70 years) or 4*CO2 (140 years) or 1.4*CO2 (34 years)
focus_1pctCO2 = '2*CO2'  # 1.4 or 2*CO2 or 4*CO2

if focus_1pctCO2 == '1.4*CO2':
    y1 = 33; y2 = 38
elif focus_1pctCO2 == '2*CO2':
    y1 = 69; y2 = 74
else:
    y1 = 134; y2 = 140

# Model agreement level
agreelev = 0.5
modelAgree = True

imodel = 1 # Choose model index in model list (modelsDef.py)
# for imodel in range(16):

# -- Choose work files

if name == 'Durack & Wijffels':
    indir = '/home/ysilvy/Density_bining/Yona_analysis/data/' #'/data/ericglod/Density_binning/Obs_Prod_density_april16/'
    file = 'DurackandWijffels_GlobalOceanChanges-NeutralDensity_1950-2000_170224_20_48_22_beta.nc' #'DurackandWijffels_GlobalOceanChanges-NeutralDensity_1950-2000_120209_11_46_11_beta.nc'
    data = indir + file
    fh2d = open_ncfile(data, 'r')


if name == 'mme_hist':
    indir = '/data/ysilvy/Density_binning/mme_hist/'
    file_2d = 'cmip5.multimodel_Nat_rcp85.historical.ensm.an.ocn.Omon.density_zon2D.nc'
    file_1d = 'cmip5.multimodel_Nat_rcp85.historical.ensm.an.ocn.Omon.density_zon1D.nc'
    data_2d = indir + file_2d
    data_1d = indir + file_1d
    fh2d = open_ncfile(data_2d, 'r')
    fh1d = open_ncfile(data_1d, 'r')


if name == 'mme_hist_histNat':
    indirh = '/data/ericglod/Density_binning/Prod_density_april15/mme_hist/'
    fileh_2d = 'cmip5.multimodel_Nat.historical.ensm.an.ocn.Omon.density_zon2D.nc'
    fileh_1d = 'cmip5.multimodel_Nat.historical.ensm.an.ocn.Omon.density_zon1D.nc'
    datah_2d = indirh + fileh_2d
    datah_1d = indirh + fileh_1d
    indirhn = '/data/ericglod/Density_binning/Prod_density_april15/mme_histNat/'
    filehn_2d = 'cmip5.multimodel_Nat.historicalNat.ensm.an.ocn.Omon.density_zon2D.nc'
    filehn_1d = 'cmip5.multimodel_Nat.historicalNat.ensm.an.ocn.Omon.density_zon1D.nc'
    datahn_2d = indirhn + filehn_2d
    datahn_1d = indirhn + filehn_1d
    fh2d = open_ncfile(datah_2d,'r')
    fh1d = open_ncfile(datah_1d,'r')
    fhn2d = open_ncfile(datahn_2d,'r')
    fhn1d = open_ncfile(datahn_1d,'r')


if name == 'ens_mean_hist' or name == '1pctCO2' or name == 'ens_mean_hist_histNat' or name == '1pctCO2vsPiC'\
        or name == 'ens_mean_rcp85_histNat':
    if name == 'ens_mean_hist' or name == 'ens_mean_hist_histNat' or name == 'ens_mean_rcp85_histNat':
        models = defModels()
        model = models[imodel]
        if name == 'ens_mean_rcp85_histNat':
            nb_members = len(model['hist-rcp85'])
        else:
            nb_members = model['props'][0]
    else :
        models = defModelsCO2piC()
        model = models[imodel]

    indir = '/data/ericglod/Density_binning/Prod_density_april15/'
    if name !='ens_mean_rcp85_histNat':
        if name == 'ens_mean_hist' or name == 'ens_mean_hist_histNat':
            file_2d = 'mme_hist/cmip5.' + model['name'] + '.historical.ensm.an.ocn.Omon.density.ver-' + model['file_end_hist'] + '_zon2D.nc'
            file_1d = 'mme_hist/cmip5.' + model['name'] + '.historical.ensm.an.ocn.Omon.density.ver-' + model['file_end_hist'] + '_zon1D.nc'
        elif name == '1pctCO2' or name == '1pctCO2vsPiC':
            file_2d = 'mme_1pctCO2/cmip5.' + model['name'] + '.1pctCO2.ensm.an.ocn.Omon.density.ver-' + model['file_end_CO2'] + '_zon2D.nc'
            file_1d = 'mme_1pctCO2/cmip5.' + model['name'] + '.1pctCO2.ensm.an.ocn.Omon.density.ver-' + model['file_end_CO2'] + '_zon1D.nc'
        data_2d = indir + file_2d
        data_1d = indir + file_1d
        fh2d = open_ncfile(data_2d, 'r')
        fh1d = open_ncfile(data_1d, 'r')

    if name == 'ens_mean_hist_histNat' or name == 'ens_mean_rcp85_histNat':
        filehn_2d = 'mme_histNat/cmip5.' + model['name'] + '.historicalNat.ensm.an.ocn.Omon.density.ver-' + model['file_end_histNat'] + '_zon2D.nc'
        filehn_1d = 'mme_histNat/cmip5.' + model['name'] + '.historicalNat.ensm.an.ocn.Omon.density.ver-' + model['file_end_histNat'] + '_zon1D.nc'
        fhn2d = open_ncfile(indir + filehn_2d, 'r')
        fhn1d = open_ncfile(indir + filehn_1d, 'r')

    if name == '1pctCO2vsPiC':
        filehn_2d = 'mme_piControl/cmip5.' + model['name'] + '.piControl.ensm.an.ocn.Omon.density.ver-' + model['file_end_piC'] + '_zon2D.nc'
        filehn_1d = 'mme_piControl/cmip5.' + model['name'] + '.piControl.ensm.an.ocn.Omon.density.ver-' + model['file_end_piC'] + '_zon1D.nc'
        fhn2d = open_ncfile(indir + filehn_2d, 'r')
        fhn1d = open_ncfile(indir + filehn_1d, 'r')

    if name == 'ens_mean_rcp85_histNat':
        filercp85_2d = glob.glob(indir+'mme_rcp85/cmip5.' + model['name'] + '.' + '*zon2D.nc')[0]
        filercp85_1d = glob.glob(indir+'mme_rcp85/cmip5.' + model['name'] + '.' + '*zon1D.nc')[0]
        fh2d = open_ncfile(filercp85_2d,'r')
        fh1d = open_ncfile(filercp85_1d,'r')


if name == 'mme_1pctCO2vsPiC' or name == 'mme_1pctCO2':
    indir_1pctCO2 = '/data/ericglod/Density_binning/Prod_density_april15/mme_1pctCO2/'
    file_2d = 'cmip5.multimodel_piCtl.1pctCO2.ensm.an.ocn.Omon.density_zon2D.nc'
    file_1d = 'cmip5.multimodel_piCtl.1pctCO2.ensm.an.ocn.Omon.density_zon1D.nc'
    data_2d = indir_1pctCO2 + file_2d
    data_1d = indir_1pctCO2 + file_1d
    fh2d = open_ncfile(data_2d,'r')
    fh1d = open_ncfile(data_1d,'r')

    if name == 'mme_1pctCO2vsPiC':
        indir_piC = '/data/ericglod/Density_binning/Prod_density_april15/mme_piControl/'
        file_2d = 'cmip5.multimodel_1pct.piControl.ensm.an.ocn.Omon.density_zon2D.nc'
        file_1d = 'cmip5.multimodel_1pct.piControl.ensm.an.ocn.Omon.density_zon1D.nc'
        data_2d = indir_piC + file_2d
        data_1d = indir_piC + file_1d
        fhn2d = open_ncfile(data_2d,'r')
        fhn1d = open_ncfile(data_1d,'r')

if name == 'mme_rcp85_histNat':
    indir_rcp85 = '/data/ericglod/Density_binning/Prod_density_april15/mme_rcp85/'
    filercp85_2d = 'cmip5.multimodel_Nat.rcp85.ensm.an.ocn.Omon.density_zon2D.nc'
    filercp85_1d = 'cmip5.multimodel_Nat.rcp85.ensm.an.ocn.Omon.density_zon1D.nc'
    indirhn = '/data/ericglod/Density_binning/Prod_density_april15/mme_histNat/'
    filehn_2d = 'cmip5.multimodel_Nat.historicalNat.ensm.an.ocn.Omon.density_zon2D.nc'
    filehn_1d = 'cmip5.multimodel_Nat.historicalNat.ensm.an.ocn.Omon.density_zon1D.nc'
    fh2d = open_ncfile(indir_rcp85 + filercp85_2d, 'r')
    fh1d = open_ncfile(indir_rcp85 + filercp85_1d, 'r')
    fhn2d = open_ncfile(indirhn + filehn_2d, 'r')
    fhn1d = open_ncfile(indirhn + filehn_1d, 'r')

# ----- Variables ------

# == Read variables ==
lat = fh2d.variables['latitude'][:]

if name == 'Durack & Wijffels':
    # -- Choose wich variable to work on
    varname = defVarDurack('salinity'); v = 'S'
    # varname = defVarDurack('temp'); v = 'T'

    density = fh2d.variables['density'][:]
    var_mean = varname['var_mean_zonal']
    var_change = varname['var_change_zonal']
    var_change_er = varname['var_change_zonal_er']

    var_attributes = fh2d.variables[var_mean]
    var_mean = fh2d.variables[var_mean][:].squeeze()
    var_change = fh2d.variables[var_change][:].squeeze()
    var_change_er = fh2d.variables[var_change_er][:].squeeze()
    # -- Determine field in 1950 and 2000 from mean field
    var_1950 = var_mean - var_change/2
    var_2000 = var_mean + var_change/2

else:
    # -- Choose wich variable to work on
#     varname = defVarmme('salinity'); v = 'S'
    # varname = defVarmme('temp'); v = 'T'
    varname = defVarmme('depth'); v = 'Z'
    density = fh2d.variables['lev'][:]
    var = varname['var_zonal_w/bowl']

    if name == 'mme_hist' or name == 'ens_mean_hist':
        if modelAgree:
            var_agree = fh2d.variables[var + 'Agree'][-5:, :, :, :]
            var_agree_a = np.ma.average(var_agree[:,1,:,:].squeeze(), axis=0)
            var_agree_p = np.ma.average(var_agree[:,2,:,:].squeeze(), axis=0)
            var_agree_i = np.ma.average(var_agree[:,3,:,:].squeeze(), axis=0)
        var = fh2d.variables[var][88:,:,:,:]  # Index 88 = year 1950
        var_mean = np.ma.average(var[0:,:,:,:], axis=0)
        var_change = np.ma.average(var[-5:,:,:,:], axis=0) - np.ma.average(var[0:5,:,:,:], axis=0)
        bowl2 = fh1d.variables['ptopsigma'][-5:,:,:]
        bowl1 = fh1d.variables['ptopsigma'][88:93,:,:]
        bowl2 = np.ma.average(bowl2, axis=0)
        bowl1 = np.ma.average(bowl1, axis=0)
        labBowl = ['1950', '2010']


    if name == 'mme_hist_histNat' or name == 'ens_mean_hist_histNat':
        varh = fh2d.variables[var][-5:, :, :, :]
        varhn = fhn2d.variables[var][-5:, :, :, :]
        var_change = np.ma.average(varh, axis=0) - np.ma.average(varhn, axis=0)
        bowl2 = fh1d.variables['ptopsigma'][-5:,:,:]
        bowl1 = fhn1d.variables['ptopsigma'][-5:,:,:]
        bowl2 = np.ma.average(bowl2, axis=0)
        bowl1 = np.ma.average(bowl1, axis=0)
        labBowl = ['histNat', 'hist']

    if name == 'mme_rcp85_histNat' or name == 'ens_mean_rcp85_histNat':
        varh = fh2d.variables[var][-10:, :, :, :]
        varhn = fhn2d.variables[var][:, :, :, :]
        var_change = np.ma.average(varh, axis=0) - np.ma.average(varhn, axis=0)
        bowl2 = fh1d.variables['ptopsigma'][-5:,:,:]
        bowl1 = fhn1d.variables['ptopsigma'][:,:,:]
        bowl2 = np.ma.average(bowl2, axis=0)
        bowl1 = np.ma.average(bowl1, axis=0)
        labBowl = ['histNat', 'RCP8.5']

    if name == 'mme_1pctCO2' or name == '1pctCO2':
        if focus_1pctCO2 == '2*CO2':
            var_end = fh2d.variables[var][69:74,:,:,:]
            var_start = fh2d.variables[var][0:5,:,:,:]
            var_change = np.ma.average(var_end, axis=0) - np.ma.average(var_start, axis=0)
            bowl2 = fh1d.variables['ptopsigma'][69:74,:,:]
            bowl1 = fh1d.variables['ptopsigma'][0:5,:,:]
            bowl2 = np.ma.average(bowl2, axis=0)
            bowl1 = np.ma.average(bowl1, axis=0)
            labBowl = ['Beginning', '2*CO2']
        if focus_1pctCO2 == '4*CO2':
            var_end = fh2d.variables[var][-5:,:,:,:]
            var_start = fh2d.variables[var][0:5,:,:,:]
            var_change = np.ma.average(var_end, axis=0) - np.ma.average(var_start, axis=0)
            bowl2 = fh1d.variables['ptopsigma'][-5:,:,:]
            bowl1 = fh1d.variables['ptopsigma'][0:5,:,:]
            bowl2 = np.ma.average(bowl2, axis=0)
            bowl1 = np.ma.average(bowl1, axis=0)
            labBowl = ['Beginning', '4*CO2']
        if focus_1pctCO2 == '1.4*CO2':
            var_end = fh2d.variables[var][33:38,:,:,:]
            var_start = fh2d.variables[var][0:5,:,:,:]
            var_change = np.ma.average(var_end, axis=0) - np.ma.average(var_start, axis=0)
            bowl2 = fh1d.variables['ptopsigma'][33:38,:,:]
            bowl1 = fh1d.variables['ptopsigma'][0:5,:,:]
            bowl2 = np.ma.average(bowl2, axis=0)
            bowl1 = np.ma.average(bowl1, axis=0)
            labBowl = ['Beginning', '1.4*CO2']

    if name == 'mme_1pctCO2vsPiC' or name == '1pctCO2vsPiC':
        if name == 'mme_1pctCO2vsPiC':
            ypiC = -20
        else:
            ypiC = -50
            y1 = y1-5
        # TODO : Simplify code below (no need to have multiple cases since y1 and y2 are defined above)
        if focus_1pctCO2 == '2*CO2':
            varCO2 = fh2d.variables[var][y1:y2,:,:,:]
            varPiC = fhn2d.variables[var][ypiC:,:,:,:]
            var_change = np.ma.average(varCO2, axis=0) - np.ma.average(varPiC, axis=0)
            bowl2 = fh1d.variables['ptopsigma'][y1:y2,:,:]
            bowl1 = fhn1d.variables['ptopsigma'][-10:,:,:]
            bowl2 = np.ma.average(bowl2, axis=0)
            bowl1 = np.ma.average(bowl1, axis=0)
            labBowl = ['PiControl', '2*CO2']
        if focus_1pctCO2 == '4*CO2':
            varCO2 = fh2d.variables[var][y1:y2,:,:,:]
            varPiC = fhn2d.variables[var][ypiC:,:,:,:]
            var_change = np.ma.average(varCO2, axis=0) - np.ma.average(varPiC, axis=0)
            bowl2 = fh1d.variables['ptopsigma'][y1:y2,:,:]
            bowl1 = fhn1d.variables['ptopsigma'][-10:,:,:]
            bowl2 = np.ma.average(bowl2, axis=0)
            bowl1 = np.ma.average(bowl1, axis=0)
            labBowl = ['PiControl', '4*CO2']
        if focus_1pctCO2 == '1.4*CO2':
            varCO2 = fh2d.variables[var][y1:y2,:,:,:]
            varPiC = fhn2d.variables[var][ypiC:,:,:,:]
            var_change = np.ma.average(varCO2, axis=0) - np.ma.average(varPiC, axis=0)
            bowl2 = fh1d.variables['ptopsigma'][y1:y2,:,:]
            bowl1 = fhn1d.variables['ptopsigma'][-10:,:,:]
            bowl2 = np.ma.average(bowl2, axis=0)
            bowl1 = np.ma.average(bowl1, axis=0)
            labBowl = ['PiControl', '1.4*CO2']

# == Define variable properties ==
minmax = varname['minmax_zonal']
clevsm = varname['clevsm_zonal']
clevsm_bold = varname['clevsm_bold']
legVar = varname['legVar']
unit = varname['unit']
if name == 'mme_rcp85_histNat' or name == 'ens_mean_rcp85_histNat':
    minmax = varname['minmax_zonal_rcp85']

# == density domain ==
rhomin = 21
rhomid = 26
rhomax = 28
domrho = [rhomin, rhomid, rhomax]

valmask = 1e+20

# == Build plot variables ==

if name == 'Durack & Wijffels':
    var_change_p = var_change[:,:,1].squeeze()
    var_change_a = var_change[:,:,2].squeeze()
    var_change_i = var_change[:,:,3].squeeze()
    var_change_er_p = var_change_er[:,:,1].squeeze()
    var_change_er_a = var_change_er[:,:,2].squeeze()
    var_change_er_i = var_change_er[:,:,3].squeeze()
    var_mean_p = var_mean[:,:,1].squeeze()
    var_mean_a = var_mean[:,:,2].squeeze()
    var_mean_i = var_mean[:,:,3].squeeze()

#     # Get bowl position
#     bowl_1950 = np.ma.masked_all((len(lat),4))
#     bowl_2000 = np.ma.masked_all((len(lat),4))
#     for ibasin in range (1,4):
#         for ilat in range(len(lat)):
#             nomask_1950 = np.ma.flatnotmasked_edges(var_1950[:, ilat, ibasin]) # returns indices of the 1st and last unmasked values
#             nomask_2000 = np.ma.flatnotmasked_edges(var_2000[:, ilat, ibasin])
#             if nomask_1950 != None:
#                 bowl_1950[ilat, ibasin] = density[nomask_1950[0]]
#             else:
#                 bowl_1950[ilat, ibasin] = np.ma.masked
#             if nomask_2000 != None:
#                 bowl_2000[ilat, ibasin] = density[nomask_2000[0]]
#             else:
#                 bowl_2000[ilat, ibasin] = np.ma.masked
#     # Bowl legend
#     labBowl = ['1950', '2000']

    # -- Create variable bundles
    varPac = {'name': 'Pacific', 'var_change': var_change_p, 'var_mean': var_mean_p, 'var_error': var_change_er_p,
              'bowl1': None, 'bowl2': None, 'labBowl': [None,None]}
    varAtl = {'name': 'Atlantic', 'var_change': var_change_a, 'var_mean': var_mean_a, 'var_error': var_change_er_a,
              'bowl1': None, 'bowl2': None, 'labBowl': [None,None]}
    varInd = {'name': 'Indian', 'var_change': var_change_i, 'var_mean': var_mean_i, 'var_error': var_change_er_i,
              'bowl1': None, 'bowl2': None, 'labBowl': [None,None]}


# if name == 'mme_hist' or name == 'ens_mean_hist':
else:
    var_change_p = var_change[2,:,:].squeeze()
    var_change_a = var_change[1,:,:].squeeze()
    var_change_i = var_change[3,:,:].squeeze()
    bowl2_p = bowl2[2,:].squeeze(); bowl1_p = bowl1[2,:].squeeze()
    bowl2_a = bowl2[1,:].squeeze(); bowl1_a = bowl1[1,:].squeeze()
    bowl2_i = bowl2[3,:].squeeze(); bowl1_i = bowl1[3,:].squeeze()
    
    if name == 'mme_hist' or name == 'ens_mean_hist':
        var_mean_p = var_mean[2,:,:].squeeze()
        var_mean_a = var_mean[1,:,:].squeeze()
        var_mean_i = var_mean[3,:,:].squeeze()
    else :
        var_mean_a = None
        var_mean_p = None
        var_mean_i = None

    # In mme 1%CO2 vs. Pi Control (and others) problem below the bowl, so take variable with bowl, and mask the data above the bowl
    # Pb with masked values in the bottom
    # Also to read model agreement, need to have variable with bowl
    var_change_a[np.ma.nonzero(var_change_a>valmask/10)] = np.ma.masked
    var_change_p[np.ma.nonzero(var_change_p>valmask/10)] = np.ma.masked
    var_change_i[np.ma.nonzero(var_change_i>valmask/10)] = np.ma.masked
    # Now mask points in the bowl
    for ilat in range(len(lat)):
        if np.ma.is_masked(bowl2_a[ilat]) == False :
            inda = np.ma.nonzero(bowl2_a[ilat]>=density)
            var_change_a[inda,ilat] = np.ma.masked
            if modelAgree and name =='mme_hist':
                var_agree_a[inda,ilat] = np.ma.masked
        if np.ma.is_masked(bowl2_p[ilat]) == False :
            indp = np.ma.nonzero(bowl2_p[ilat]>=density)
            var_change_p[indp,ilat] = np.ma.masked
            if modelAgree and name =='mme_hist':
                var_agree_p[indp,ilat] = np.ma.masked
        if np.ma.is_masked(bowl2_i[ilat]) == False :
            indi = np.ma.nonzero(bowl2_i[ilat]>=density)
            var_change_i[indi,ilat] = np.ma.masked
            if modelAgree and name =='mme_hist':
                var_agree_i[indi,ilat] = np.ma.masked

    # -- Create variable bundles
    varPac = {'name': 'Pacific', 'var_change': var_change_p, 'var_mean': var_mean_p,
              'bowl1': bowl1_p, 'bowl2': bowl2_p, 'labBowl': labBowl}
    varAtl = {'name': 'Atlantic', 'var_change': var_change_a, 'var_mean': var_mean_a,
              'bowl1': bowl1_a, 'bowl2': bowl2_a, 'labBowl': labBowl}
    varInd = {'name': 'Indian', 'var_change': var_change_i, 'var_mean': var_mean_i,
              'bowl1': bowl1_i, 'bowl2': bowl2_i, 'labBowl': labBowl}


# if name == 'mme_hist_histNat' or name == 'mme_1pctCO2vsPiC' or name == 'mme_1pctCO2' or name == '1pctCO2' \
#         or name == 'ens_mean_hist_histNat' or name == '1pctCO2vsPiC' or name == 'mme_rcp85_histNat' \
#         or name == 'ens_mean_rcp85_histNat':
#     var_change_p = var_change[2, :, :].squeeze()
#     var_change_a = var_change[1, :, :].squeeze()
#     var_change_i = var_change[3, :, :].squeeze()
#     bowl2_p = bowl2[2,:].squeeze(); bowl1_p = bowl1[2,:].squeeze()
#     bowl2_a = bowl2[1,:].squeeze(); bowl1_a = bowl1[1,:].squeeze()
#     bowl2_i = bowl2[3,:].squeeze(); bowl1_i = bowl1[3,:].squeeze()
#     var_mean_a = None
#     var_mean_p = None
#     var_mean_i = None

    # # In mme 1%CO2 vs. Pi Control (and others) problem below the bowl, so take variable with bowl, and mask the data above the bowl
    # # Pb with masked values in the bottom
    # var_change_a[np.ma.nonzero(var_change_a>valmask/10)] = np.ma.masked
    # var_change_p[np.ma.nonzero(var_change_p>valmask/10)] = np.ma.masked
    # var_change_i[np.ma.nonzero(var_change_i>valmask/10)] = np.ma.masked
    # # Now mask points in the bowl
    # for ilat in range(len(lat)):
    #     if np.ma.is_masked(bowl2_a[ilat]) == False :
    #         inda = np.ma.nonzero(bowl2_a[ilat]>=density)
    #         var_change_a[inda,ilat] = np.ma.masked
    #     if np.ma.is_masked(bowl2_p[ilat]) == False :
    #         indp = np.ma.nonzero(bowl2_p[ilat]>=density)
    #         var_change_p[indp,ilat] = np.ma.masked
    #     if np.ma.is_masked(bowl2_i[ilat]) == False :
    #         indi = np.ma.nonzero(bowl2_i[ilat]>=density)
    #         var_change_i[indi,ilat] = np.ma.masked
    #

    # # -- Create variable bundles
    # varPac = {'name': 'Pacific', 'var_change': var_change_p, 'var_mean': None,
    #           'bowl1': bowl1_p, 'bowl2': bowl2_p, 'labBowl': labBowl}
    # varAtl = {'name': 'Atlantic', 'var_change': var_change_a, 'var_mean': None,
    #           'bowl1': bowl1_a, 'bowl2': bowl2_a, 'labBowl': labBowl}
    # varInd = {'name': 'Indian', 'var_change': var_change_i, 'var_mean': None,
    #           'bowl1': bowl1_i, 'bowl2': bowl2_i, 'labBowl': labBowl}

# ------------------------------------
#               Plot
# ------------------------------------

fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(17,5))


if name == 'Durack & Wijffels':
    levels = np.linspace(minmax[0], minmax[1], minmax[2])
    cmap = custom_div_cmap()

    cnplot = zonal_2D(plt, 'total', axes[0,0], axes[1,0], 'left', lat, density, varAtl, domrho, cmap, levels, clevsm, clevsm_bold)

    cnplot = zonal_2D(plt, 'total', axes[0,1], axes[1,1], 'mid', lat, density, varPac, domrho, cmap, levels, clevsm, clevsm_bold)

    cnplot = zonal_2D(plt, 'total', axes[0,2], axes[1,2], 'right', lat, density, varInd, domrho, cmap, levels, clevsm, clevsm_bold)

else:
    levels = np.linspace(minmax[0], minmax[1], minmax[2])
    cmap = custom_div_cmap() # plt.get_cmap('bwr')

    cnplot = zonal_2D(plt, 'total_mme', axes[0, 0], axes[1, 0], 'left', lat, density, varAtl, domrho, cmap, levels, clevsm, clevsm_bold)

    cnplot = zonal_2D(plt, 'total_mme', axes[0, 1], axes[1, 1], 'mid', lat, density, varPac, domrho, cmap, levels, clevsm, clevsm_bold)

    cnplot = zonal_2D(plt, 'total_mme', axes[0, 2], axes[1, 2], 'right', lat, density, varInd, domrho, cmap, levels, clevsm, clevsm_bold)

if modelAgree and name=='mme_hist' :
    modelagree(axes[0,0],axes[1,0],agreelev,lat,density,var_agree_a)
    modelagree(axes[0,1],axes[1,1],agreelev,lat,density,var_agree_p)
    modelagree(axes[0,2],axes[1,2],agreelev,lat,density,var_agree_i)

plt.subplots_adjust(hspace=.012, wspace=0.05, left=0.04, right=0.86)

cb = plt.colorbar(cnplot[0], ax=axes.ravel().tolist(), ticks=levels[::3], fraction=0.015, shrink=2.0, pad=0.05)
cb.set_label('%s (%s)' % (legVar, unit), fontweight='bold')

#plt.figtext(.006,.96,'b',fontweight='bold',fontsize=16)

if name == 'mme_hist' or name == 'Durack & Wijffels':
    plotTitle = '%s changes (2000-1950), %s' %(legVar, name)
    plotName = name + '_' + v + 'changes'
    if modelAgree and name == 'mme_hist':
        plotName = name + '_' + v + 'changes_modelAgree'
    if name == 'mme_hist':
        figureDir = 'models/zonal_ys/hist/'
    else:
        figureDir = 'obs/zonal_ys/'

elif name == 'mme_hist_histNat':
    plotTitle = '%s changes %s' %(legVar, name)
    plotName = name + '_' + legVar
    figureDir = 'models/zonal_ys/hist-histNat/'

elif name == 'mme_rcp85_histNat':
    plotTitle = '%s changes %s' %(legVar, name)
    plotName = name + '_' + legVar + '_meanhistNat'
    figureDir = 'models/zonal_ys/rcp85-histNat/'

elif name == 'mme_1pctCO2vsPiC':
    plotTitle = '%s changes %s (%s)' %(legVar, name, focus_1pctCO2)
    plotName = name + '_' + focus_1pctCO2 + '_' + legVar
    figureDir = 'models/zonal_ys/1pctCO2vsPiC/'

elif name == '1pctCO2vsPiC':
    plotTitle = '%s changes 1pctCO2 vs. pi Control (%s, %s)' %(legVar, model['name'], focus_1pctCO2)
    plotName = model['name'] + '_' + name + '_' + focus_1pctCO2 + '_' + legVar
    figureDir = 'models/zonal_ys/1pctCO2vsPiC/With_grid/'

elif name == 'mme_1pctCO2':
    plotTitle = '%s changes (%s, %s)' %(legVar, name, focus_1pctCO2)
    plotName = name + '_' + focus_1pctCO2 + '_' + v + 'changes'
    figureDir = 'models/zonal_ys/1pctCO2/'

elif name == '1pctCO2':
    plotTitle = '%s changes (%s, %s)' %(legVar, model['name'], focus_1pctCO2)
    plotName = model['name'] + '_' + focus_1pctCO2 + '_' + v + 'changes'
    figureDir = 'models/zonal_ys/1pctCO2/'

elif name == 'ens_mean_hist_histNat':
    plotTitle = '%s changes hist vs. histNat, %s ensemble mean (%d members)' %(legVar, model['name'], nb_members)
    plotName = model['name'] + '_hist_histNat' + legVar
    figureDir = 'models/zonal_ys/hist-histNat/'

elif name == 'ens_mean_rcp85_histNat':
    if nb_members>1:
        nb = '%d members'%(nb_members,)
    else:
        nb = '%d member'%(nb_members,)
    plotTitle = '%s changes RCP8.5 - histNat, %s ensemble mean (%s)' %(legVar, model['name'], nb)
    plotName = model['name'] + '_rcp85_histNat' + legVar
    figureDir = 'models/zonal_ys/rcp85-histNat/With_grid/'

else:
    plotTitle = '%s changes (2000-1950), %s ensemble mean (%d members)' %(legVar, model['name'], nb_members)
    plotName = model['name'] + '_' + v + 'changes'
    figureDir = 'models/zonal_ys/hist/'

# Date
now = datetime.datetime.now()
date = now.strftime("%Y-%m-%d")

#plt.suptitle(plotTitle, fontweight='bold', fontsize=14, verticalalignment='top')
plt.figtext(.006,.45,'Density (kg.m-3)',rotation='vertical',horizontalalignment='left',fontweight='bold')
#plt.figtext(.5,.015,'Computed by : zonal_ys_changes.py,  '+date,fontsize=9,ha='center')
#if modelAgree and name=='mme_hist':
    #plt.figtext(.1,.015,'Model agreement level : ' + str(agreelev),fontsize=9,ha='center')

# plt.show()
plt.savefig('/home/ysilvy/figures/'+figureDir+plotName+'.png', bbox_inches='tight')

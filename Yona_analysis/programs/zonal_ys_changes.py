#!/bin/env python
# -*- coding: utf-8 -*-

"""
Python matplotlib
Plot density/latitude salinity changes for each basin
Use for mme rather than Durack and Wijffels data
(for the latter, see the script "isopycnal_migration_2D.py")


"""

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset as open_ncfile
from maps_matplot_lib import defVarDurack, zonal_2D, defVarmme, custom_div_cmap
from modelsDef import defModels, defModelsCO2piC

# ----- Workspace ------

#name = 'Durack & Wijffels'
#name = 'mme_hist'
#name = 'mme_hist_histNat'
#name = 'ens_mean_hist'
#name = 'mme_1pctCO2vsPiC'
name = 'mme_1pctCO2'
#name = '1pctCO2'

focus_1pctCO2 = '2*CO2' # 2 or 4*CO2

imodel = 5 # Choose model index in model list (modelsDef.py)

# -- Choose work files

if name == 'Durack & Wijffels':
    indir = '/data/ericglod/Density_binning/Obs_Prod_density_april16/'
    file = 'DurackandWijffels_GlobalOceanChanges-NeutralDensity_1950-2000_120209_11_46_11_beta.nc'
    data = indir + file
    fh2d = open_ncfile(data, 'r')


if name == 'mme_hist':
    indir = '/data/ericglod/Density_binning/Prod_density_april15/mme_hist/'
    file_2d = 'cmip5.multimodel_All.historical.ensm.an.ocn.Omon.density_zon2D.nc'
    file_1d = 'cmip5.multimodel_All.historical.ensm.an.ocn.Omon.density_zon1D.nc'
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
    filehn_2d = 'cmip5.multimodel_All.historicalNat.ensm.an.ocn.Omon.density_zon2D.nc'
    filehn_1d = 'cmip5.multimodel_All.historicalNat.ensm.an.ocn.Omon.density_zon1D.nc'
    datahn_2d = indirhn + filehn_2d
    datahn_1d = indirhn + filehn_1d
    fh2d = open_ncfile(datah_2d,'r')
    fh1d = open_ncfile(datah_1d,'r')
    fhn_2d = open_ncfile(datahn_2d,'r')
    fhn_1d = open_ncfile(datahn_1d,'r')


if name == 'ens_mean_hist' or name == '1pctCO2':
    if name == 'ens_mean_hist':
        models = defModels()
        model = models[imodel]
        nb_members = model['props'][0]
    else :
        models = defModelsCO2piC()
        model = models[imodel]

    indir = '/data/ericglod/Density_binning/Prod_density_april15/'
    if name == 'ens_mean_hist':
        file_2d = 'mme_hist/cmip5.' + model['name'] + '.historical.ensm.an.ocn.Omon.density.ver-' + model['file_end'] + '_zon2D.nc'
        file_1d = 'mme_hist/cmip5.' + model['name'] + '.historical.ensm.an.ocn.Omon.density.ver-' + model['file_end'] + '_zon1D.nc'
    else:
        file_2d = 'mme_1pctCO2/cmip5.' + model['name'] + '.1pctCO2.ensm.an.ocn.Omon.density.ver-' + model['file_end_CO2'] + '_zon2D.nc'
        file_1d = 'mme_1pctCO2/cmip5.' + model['name'] + '.1pctCO2.ensm.an.ocn.Omon.density.ver-' + model['file_end_CO2'] + '_zon1D.nc'
    data_2d = indir + file_2d
    data_1d = indir + file_1d
    fh2d = open_ncfile(data_2d, 'r')
    fh1d = open_ncfile(data_1d, 'r')


if name == 'mme_1pctCO2vsPiC' or name == 'mme_1pctCO2':
    indir_1pctCO2 = '/data/ericglod/Density_binning/Prod_density_april15/mme_1pctCO2/'
    file_2d = 'zon2D_mean_bowl/cmip5.multimodel_piCtl.1pctCO2.ensm.an.ocn.Omon.density_zon2D.nc'
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



# ----- Variables ------

# == Read variables ==
lat = fh2d.variables['latitude'][:]

if name == 'Durack & Wijffels':
    # -- Choose wich variable to work on
    #varname = defVarDurack('salinity')
    varname = defVarDurack('temp')

    density = fh2d.variables['density'][:]
    var_mean = varname['var_mean_zonal']
    var_change = varname['var_change_zonal']
    var_change_er = varname['var_change_zonal_er']

    var_attributes = fh2d.variables[var_mean]
    var_mean = fh2d.variables[var_mean][:].squeeze()
    var_change = fh2d.variables[var_change][:].squeeze()
    var_change_er = fh2d.variables[var_change_er][:].squeeze()

else:
    # -- Choose wich variable to work on
    varname = defVarmme('salinity'); v = 'S'
    #varname = defVarmme('temp'); v = 'T'
    #varname= defVarmme('depth'); v = 'Z'
    density = fh2d.variables['lev'][:]
    var = varname['var_zonal']

    if name == 'mme_hist' or name == 'ens_mean_hist':
        var = fh2d.variables[var][88:,:,:,:] # Index 88 = year 1950
        var_mean = np.ma.average(var[0:,:,:,:], axis=0)
        var_change = np.ma.average(var[-5:,:,:,:], axis=0) - np.ma.average(var[0:5,:,:,:], axis=0)
        bowl2 = fh1d.variables['ptopsigma'][-5:,:,:]
        bowl1 = fh1d.variables['ptopsigma'][88:93,:,]
        bowl2 = np.ma.average(bowl2, axis=0)
        bowl1 = np.ma.average(bowl1, axis=0)
        labBowl = ['1950', '2010']

    if name == 'mme_hist_histNat' or name == 'mme_1pctCO2vsPiC':
        varh = fh2d.variables[var][-5:, :, :, :]
        varhn = fhn2d.variables[var][-5:, :, :, :]
        var_change = np.ma.average(varh, axis=0) - np.ma.average(varhn, axis=0)
        bowl2 = fh1d.variables['ptopsigma'][-5:,:,:]
        bowl1 = fhn1d.variables['ptopsigma'][-5:,:,:]
        bowl2 = np.ma.average(bowl2, axis=0)
        bowl1 = np.ma.average(bowl1, axis=0)
        if name == 'mme_hist_histNat':
            labBowl = ['histNat', 'hist']
        else:
            labBowl = ['piControl', '1pctCO2']

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


# == Define variable properties ==
minmax = varname['minmax_zonal']
clevsm = varname['clevsm_zonal']
clevsm_bold = varname['clevsm_bold']
legVar = varname['legVar']
unit = varname['unit']


# == density domain ==
rhomin = 21
rhomid = 26
rhomax = 28
domrho = [rhomin, rhomid, rhomax]


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
    # -- Create variable bundles
    varPac = {'name': 'Pacific', 'var_change': var_change_p, 'var_mean': var_mean_p, 'var_error': var_change_er_p}
    varAtl = {'name': 'Atlantic', 'var_change': var_change_a, 'var_mean': var_mean_a, 'var_error': var_change_er_a}
    varInd = {'name': 'Indian', 'var_change': var_change_i, 'var_mean': var_mean_i, 'var_error': var_change_er_i}

if name == 'mme_hist' or name == 'ens_mean_hist':
    var_change_p = var_change[2,:,:].squeeze()
    var_change_a = var_change[1,:,:].squeeze()
    var_change_i = var_change[3,:,:].squeeze()
    var_mean_p = var_mean[2,:,:].squeeze()
    var_mean_a = var_mean[1,:,:].squeeze()
    var_mean_i = var_mean[3,:,:].squeeze()
    bowl2_p = bowl2[2,:].squeeze(); bowl1_p = bowl1[2,:].squeeze()
    bowl2_a = bowl2[1,:].squeeze(); bowl1_a = bowl1[1,:].squeeze()
    bowl2_i = bowl2[3,:].squeeze(); bowl1_i = bowl1[3,:].squeeze()
    # -- Create variable bundles
    varPac = {'name': 'Pacific', 'var_change': var_change_p, 'var_mean': var_mean_p,
              'bowl1': bowl1_p, 'bowl2': bowl2_p, 'labBowl': labBowl}
    varAtl = {'name': 'Atlantic', 'var_change': var_change_a, 'var_mean': var_mean_a,
              'bowl1': bowl1_p, 'bowl2': bowl2_p, 'labBowl': labBowl}
    varInd = {'name': 'Indian', 'var_change': var_change_i, 'var_mean': var_mean_i,
              'bowl1': bowl1_p, 'bowl2': bowl2_p, 'labBowl': labBowl}

if name == 'mme_hist_histNat' or name == 'mme_1pctCO2vsPiC' or name == 'mme_1pctCO2' or name == '1pctCO2':
    var_change_p = var_change[2, :, :].squeeze()
    var_change_a = var_change[1, :, :].squeeze()
    var_change_i = var_change[3, :, :].squeeze()
    bowl2_p = bowl2[2,:].squeeze(); bowl1_p = bowl1[2,:].squeeze()
    bowl2_a = bowl2[1,:].squeeze(); bowl1_a = bowl1[1,:].squeeze()
    bowl2_i = bowl2[3,:].squeeze(); bowl1_i = bowl1[3,:].squeeze()
    # -- Create variable bundles
    varPac = {'name': 'Pacific', 'var_change': var_change_p, 'var_mean': None,
              'bowl1': bowl1_p, 'bowl2': bowl2_p, 'labBowl': labBowl}
    varAtl = {'name': 'Atlantic', 'var_change': var_change_a, 'var_mean': None,
              'bowl1': bowl1_a, 'bowl2': bowl2_a, 'labBowl': labBowl}
    varInd = {'name': 'Indian', 'var_change': var_change_i, 'var_mean': None,
              'bowl1': bowl1_i, 'bowl2': bowl2_i, 'labBowl': labBowl}

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
    cmap = custom_div_cmap() #plt.get_cmap('bwr')

    cnplot = zonal_2D(plt, 'total_mme', axes[0, 0], axes[1, 0], 'left', lat, density, varAtl, domrho, cmap, levels, clevsm, clevsm_bold)

    cnplot = zonal_2D(plt, 'total_mme', axes[0, 1], axes[1, 1], 'mid', lat, density, varPac, domrho, cmap, levels, clevsm, clevsm_bold)

    cnplot = zonal_2D(plt, 'total_mme', axes[0, 2], axes[1, 2], 'right', lat, density, varInd, domrho, cmap, levels, clevsm, clevsm_bold)


plt.subplots_adjust(hspace=.0001, wspace=0.05, left=0.04, right=0.86)

cb = plt.colorbar(cnplot, ax=axes.ravel().tolist(), ticks=levels[::3], fraction=0.015, shrink=2.0, pad=0.05)
cb.set_label('%s (%s)' % (legVar, unit), fontweight='bold')


if name == 'mme_hist' or name == 'Durack & Wijffels':
    plotTitle = '%s changes (2000-1950), %s' %(legVar, name)
    plotName = name + '_' + v + 'changes'
    figureDir = 'hist/' # save only models here

elif name == 'mme_hist_histNat':
    plotTitle = '%s changes %s (last 5 years), ' %(legVar, name)

elif name == 'mme_1pctCO2vsPiC':
    plotTitle = '%s changes 1pctCO2 vs pi Control (last 5 years)' %(legVar, )
    plotName = name + '_' + legVar
    figureDir = '1pctCO2vsPiC/'

elif name == 'mme_1pctCO2':
    plotTitle = '%s changes (%s, %s)' %(legVar, name, focus_1pctCO2)
    plotName = name + '_' + focus_1pctCO2 + '_' + v + 'changes'
    figureDir = '1pctCO2/'

elif name == '1pctCO2':
    plotTitle = '%s changes (%s, %s)' %(legVar, model['name'], focus_1pctCO2)
    plotName = model['name'] + '_' + focus_1pctCO2 + '_' + v + 'changes'
    figureDir = '1pctCO2/With_grid'

else:
    plotTitle = '%s changes (2000-1950), %s ensemble mean (%d members)' %(legVar, model['name'], nb_members)
    plotName = model['name'] + '_' + v + 'changes'
    figureDir = 'hist/'

plt.suptitle(plotTitle, fontweight='bold', fontsize=14, verticalalignment='top')

plt.show()

#plt.savefig('/home/ysilvy/Density_bining/Yona_analysis/figures/models/zonal_ys/'+figureDir+plotName+'.png', bbox_inches='tight')
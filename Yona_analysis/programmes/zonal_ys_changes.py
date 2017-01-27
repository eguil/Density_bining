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
from modelsDef import defModels

# ----- Workspace ------

#name = 'Durack & Wijffels'
#name = 'mme_hist'
#name = 'mme_hist_histNat'
#name = 'ens_mean_hist'
#name = '1pctCO2vsPiC'
#name = 'mme_1pctCO2'
name = 'ens_mean_1pctCO2'

focus_1pctCO2 = '2*CO2' # 2 or 4*CO2


if name == 'Durack & Wijffels':
    indir = '/data/ericglod/Density_binning/Obs_Prod_density_april16/'
    file = 'DurackandWijffels_GlobalOceanChanges-NeutralDensity_1950-2000_120209_11_46_11_beta.nc'
    data = indir + file
    fh = open_ncfile(data, 'r')

if name == 'mme_hist':
    indir = '/data/ericglod/Density_binning/Prod_density_april15/mme_hist/'
    file = 'cmip5.multimodel_All.historical.ensm.an.ocn.Omon.density_zon2D.nc'
    data = indir + file
    fh = open_ncfile(data, 'r')

if name == 'mme_hist_histNat':
    indirh = '/data/ericglod/Density_binning/Prod_density_april15/mme_hist/'
    fileh = 'cmip5.multimodel_Nat.historical.ensm.an.ocn.Omon.density_zon2D.nc'
    datah = indirh + fileh
    indirhn = '/data/ericglod/Density_binning/Prod_density_april15/mme_histNat/'
    filehn = 'cmip5.multimodel_All.historicalNat.ensm.an.ocn.Omon.density_zon2D.nc'
    datahn = indirhn + filehn
    fh = open_ncfile(datah,'r')
    fhn = open_ncfile(datahn,'r')

if name == 'ens_mean_hist' or name == 'ens_mean_1pctCO2':
    models = defModels()
    model = models[14] #Iterate
    nb_members = model['props'][0]

    indir = '/data/ericglod/Density_binning/Prod_density_april15/'
    if name == 'ens_mean_hist':
        file = 'mme_hist/cmip5.' + model['name'] + '.historical.ensm.an.ocn.Omon.density.ver-' + model['file_end'] + '_zon2D.nc'
    else:
        file = 'mme_1pctCO2/cmip5.' + model['name'] + '.1pctCO2.ensm.an.ocn.Omon.density.ver-' + model['file_end'] + '_zon2D.nc'
    data = indir + file
    fh = open_ncfile(data, 'r')

if name == '1pctCO2vsPiC':
    indir_1pctCO2 = '/data/ericglod/Density_binning/Prod_density_april15/mme_1pctCO2/'
    file_1pctCO2 = 'cmip5.multimodel_All.1pctCO2.ensm.an.ocn.Omon.density_zon2D.nc'
    data_1pctCO2 = indir_1pctCO2 + file_1pctCO2
    indir_piC = '/data/ericglod/Density_binning/Prod_density_april15/mme_piControl/'
    file_piC = 'cmip5.' + model['name'] + '.piControl.ensm.an.ocn.Omon.density.ver-' + model['file_end'] + '_zon2D.nc'
    data_piC = indir_piC + file_piC
    fh = open_ncfile(data_1pctCO2,'r')
    fhn = open_ncfile(data_piC,'r')

if name == 'mme_1pctCO2':
    indir_1pctCO2 = '/data/ericglod/Density_binning/Prod_density_april15/mme_1pctCO2/'
    file_1pctCO2 = 'cmip5.GFDL-ESM2M.1pctCO2.ensm.an.ocn.Omon.density.ver-v20130226_zon2D.nc' #'cmip5.multimodel_All.1pctCO2.ensm.an.ocn.Omon.density_zon2D.nc'
    data_1pctCO2 = indir_1pctCO2 + file_1pctCO2
    fh = open_ncfile(data_1pctCO2,'r')


# ----- Variables ------

# Read variables
lat = fh.variables['latitude'][:]

if name == 'Durack & Wijffels':
    #varname = defVarDurack('salinity')
    varname = defVarDurack('temp')
    density = fh.variables['density'][:]
    var_mean = varname['var_mean_zonal']
    var_change = varname['var_change_zonal']
    var_change_er = varname['var_change_zonal_er']

    var_attributes = fh.variables[var_mean]
    var_mean = fh.variables[var_mean][:].squeeze()
    var_change = fh.variables[var_change][:].squeeze()
    var_change_er = fh.variables[var_change_er][:].squeeze()

else:
    varname = defVarmme('salinity'); v = 'S'
    #varname = defVarmme('temp'); v = 'T'
    density = fh.variables['lev'][:]
    var = varname['var_zonal']

    if name == 'mme_hist' or name == 'ens_mean_hist':
        var = fh.variables[var][88:,:,:,:] # Index 88 = year 1950
        var_mean = np.ma.average(var[0:,:,:,:], axis=0)
        var_change = np.ma.average(var[-5:,:,:,:], axis=0) - np.ma.average(var[0:5,:,:,:], axis=0)
    if name == 'mme_hist_histNat' or name == '1pctCO2vsPiC':
        varh = fh.variables[var][-5:, :, :, :]
        varhn = fhn.variables[var][-5:, :, :, :]
        var_change = np.ma.average(varh, axis=0) - np.ma.average(varhn, axis=0)
    if name == 'mme_1pctCO2' or name == 'ens_mean_1pctCO2':
        if focus_1pctCO2 == '2*CO2':
            var_end = fh.variables[var][69:75,:,:,:]
            var_start = fh.variables[var][0:5,:,:,:]
            var_change = np.ma.average(var_end, axis=0) - np.ma.average(var_start, axis=0)
        if focus_1pctCO2 == '4*CO2':
            var_end = fh.variables[var][-5:,:,:,:]
            var_start = fh.variables[var][0:5,:,:,:]
            var_change = np.ma.average(var_end, axis=0) - np.ma.average(var_start, axis=0)


# Define variable properties
minmax = varname['minmax_zonal']
clevsm = varname['clevsm_zonal']
clevsm_bold = varname['clevsm_bold']
legVar = varname['legVar']
unit = varname['unit']


# density domain
rhomin = 21
rhomid = 26
rhomax = 28
domrho = [rhomin, rhomid, rhomax]


# Build plot variables

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
    # -- Create variable bundles
    varPac = {'name': 'Pacific', 'var_change': var_change_p, 'var_mean': var_mean_p}
    varAtl = {'name': 'Atlantic', 'var_change': var_change_a, 'var_mean': var_mean_a}
    varInd = {'name': 'Indian', 'var_change': var_change_i, 'var_mean': var_mean_i}

if name == 'mme_hist_histNat' or name == '1pctCO2vsPiC' or name == 'mme_1pctCO2' or name == 'ens_mean_1pctCO2':
    var_change_p = var_change[2, :, :].squeeze()
    var_change_a = var_change[1, :, :].squeeze()
    var_change_i = var_change[3, :, :].squeeze()
    # -- Create variable bundles
    varPac = {'name': 'Pacific', 'var_change': var_change_p, 'var_mean': None}
    varAtl = {'name': 'Atlantic', 'var_change': var_change_a, 'var_mean': None}
    varInd = {'name': 'Indian', 'var_change': var_change_i, 'var_mean': None}

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
    cmap = custom_div_cmap() #plt.get_cmap('bwr') #

    cnplot = zonal_2D(plt, 'total_mme', axes[0, 0], axes[1, 0], 'left', lat, density, varAtl, domrho, cmap, levels, clevsm, clevsm_bold)

    cnplot = zonal_2D(plt, 'total_mme', axes[0, 1], axes[1, 1], 'mid', lat, density, varPac, domrho, cmap, levels, clevsm, clevsm_bold)

    cnplot = zonal_2D(plt, 'total_mme', axes[0, 2], axes[1, 2], 'right', lat, density, varInd, domrho, cmap, levels, clevsm, clevsm_bold)


plt.subplots_adjust(hspace=.0001, wspace=0.05, left=0.04, right=0.86)

cb = plt.colorbar(cnplot, ax=axes.ravel().tolist(), ticks=levels[::3], fraction=0.015, shrink=2.0, pad=0.05)
cb.set_label('%s (%s)' % (legVar, unit), fontweight='bold')


if name == 'mme_hist' or name == 'Durack & Wijffels':
    plotTitle = '%s changes (2000-1950), %s' %(legVar, name)

elif name == 'mme_hist_histNat':
    plotTitle = '%s changes %s (last 5 years), ' %(legVar, name)

elif name == '1pctCO2vsPiC':
    plotTitle = '%s changes (1pctCO2 vs pi Control), %s ensemble mean (%d members)' %(legVar, model['name'], nb_members)
    plotName = model['name'] + name + '_' + legVar
    figureDir = '1pctCO2vsPiC/'

elif name == 'mme_1pctCO2':
    plotTitle = '%s changes (%s, %s)' %(legVar, name, focus_1pctCO2)

elif name == 'ens_mean_1pctCO2':
    plotTitle = '%s changes (%s, %s), %s ensemble mean (%d members)' %(legVar, name, focus_1pctCO2, model['name'], nb_members)
    plotName = model['name'] + '_' + focus_1pctCO2 + '_' + v + 'changes'
    figureDir = '1pctCO2/'

else:
    plotTitle = '%s changes (2000-1950), %s ensemble mean (%d members)' %(legVar, model['name'], nb_members)
    plotName = model['name'] + '_' + v + 'changes'
    figureDir = 'hist/'

plt.suptitle(plotTitle, fontweight='bold', fontsize=14, verticalalignment='top')

#plt.show()


plt.savefig('/home/ysilvy/Density_bining/Yona_analysis/figures/models/zonal_ys/'+figureDir+plotName+'.png', bbox_inches='tight')
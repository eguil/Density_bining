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
from maps_matplot_lib import defVarDurack, zonal_2D, defVarmme
from modelsDef import defModels

# ----- Workspace ------

#name = 'Durack & Wijffels'
#name = 'mme_hist'
#name = 'mme_hist_histNat'
name = 'ens_mean_hist'

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

if name == 'ens_mean_hist':
    models = defModels()
    model = models[1] #Iterate
    nb_members = model['props'][0]

    indir = '/data/ericglod/Density_binning/Prod_density_april15/mme_hist/'
    file = 'cmip5.' + model['name'] + '.historical.ensm.an.ocn.Omon.density.ver-' + model['file_end'] + '_zon2D.nc'
    data = indir + file
    fh = open_ncfile(data, 'r')


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
    if name == 'mme_hist_histNat':
        varh = fh.variables[var][-5:, :, :, :]
        varhn = fhn.variables[var][-5:, :, :, :]
        var_change = np.ma.average(varh, axis=0) - np.ma.average(varhn, axis=0)

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

if name == 'mme_hist_histNat':
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
    cnplot = zonal_2D(plt, 'total', axes[0,0], axes[1,0], 'left', lat, density, varAtl, minmax, domrho, clevsm, clevsm_bold)

    cnplot = zonal_2D(plt, 'total', axes[0,1], axes[1,1], 'mid', lat, density, varPac, minmax, domrho, clevsm, clevsm_bold)

    cnplot = zonal_2D(plt, 'total', axes[0,2], axes[1,2], 'right', lat, density, varInd, minmax, domrho, clevsm, clevsm_bold)

else:
    cnplot = zonal_2D(plt, 'total_mme', axes[0, 0], axes[1, 0], 'left', lat, density, varAtl, minmax, domrho, clevsm, clevsm_bold)

    cnplot = zonal_2D(plt, 'total_mme', axes[0, 1], axes[1, 1], 'mid', lat, density, varPac, minmax, domrho, clevsm, clevsm_bold)

    cnplot = zonal_2D(plt, 'total_mme', axes[0, 2], axes[1, 2], 'right', lat, density, varInd, minmax, domrho, clevsm, clevsm_bold)

    if name == 'ens_mean_hist':
        name = model['name']
    plotName = name + '_' + v + 'changes'


plt.subplots_adjust(hspace=.0001, wspace=0.05, left=0.04, right=0.86)

cb = plt.colorbar(cnplot[0], ax=axes.ravel().tolist(), ticks=cnplot[1][::3])
cb.set_label('%s (%s)' % (legVar, unit), fontweight='bold')

if name == 'mme_hist' or name == 'Durack & Wijffels':
    plt.suptitle('%s changes (2000-1950), %s' %(legVar, name),
          fontweight='bold', fontsize=14, verticalalignment='top')

elif name == 'mme_hist_histNat':
    plt.suptitle('%s changes %s (last 5 years), ' %(legVar, name),
          fontweight='bold', fontsize=14, verticalalignment='top')

else:
    plt.suptitle('%s changes (2000-1950), %s ensemble mean (%d members)' %(legVar, name, nb_members), fontweight='bold',
                 fontsize=14, va='top')

#plt.show()

plt.savefig('/home/ysilvy/Density_bining/Yona_analysis/figures/models/zonal_ys/hist/'+plotName+'.png', bbox_inches='tight')
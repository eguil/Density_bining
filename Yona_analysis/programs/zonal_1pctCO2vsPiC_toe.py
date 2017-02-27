#!/bin/env python
# -*- coding: utf-8 -*-

"""
Python matplotlib
Plot density/latitude Time of Emergence for several variables

"""

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset as open_ncfile
from maps_matplot_lib import zonal_2D, defVarmme, custom_div_cmap
from modelsDef import defModelsCO2piC
from libToE import findToE
import colormaps as cmaps

# ----- Workspace ------

#imodel = 6
for imodel in range(13,17):

    #name = 'mme'
    name = 'single_model'

    indir_1pctCO2 = '/data/ericglod/Density_binning/Prod_density_april15/mme_1pctCO2/'
    indir_piC = '/data/ericglod/Density_binning/Prod_density_april15/mme_piControl/'


    if name == 'mme':
        file2d_1pctCO2 = 'zon2D_mean_bowl/cmip5.multimodel_piCtl.1pctCO2.ensm.an.ocn.Omon.density_zon2D.nc'
        file1d_1pctCO2 = 'cmip5.multimodel_piCtl.1pctCO2.ensm.an.ocn.Omon.density_zon1D.nc'
        file2d_piC = 'cmip5.multimodel_1pct.piControl.ensm.an.ocn.Omon.density_zon2D.nc'
        file1d_piC = 'cmip5.multimodel_1pct.piControl.ensm.an.ocn.Omon.density_zon1D.nc'

    else :
        models = defModelsCO2piC()
        model = models[imodel ] # Iterate

        file2d_1pctCO2 = 'cmip5.' + model['name'] + '.1pctCO2.ensm.an.ocn.Omon.density.ver-' + model['file_end_CO2'] + '_zon2D.nc'
        file1d_1pctCO2 = 'cmip5.' + model['name'] + '.1pctCO2.ensm.an.ocn.Omon.density.ver-' + model['file_end_CO2'] + '_zon1D.nc'
        file2d_piC = 'cmip5.' + model['name'] + '.piControl.ensm.an.ocn.Omon.density.ver-' + model['file_end_piC'] + '_zon2D.nc'
        file1d_piC = 'cmip5.' + model['name'] + '.piControl.ensm.an.ocn.Omon.density.ver-' + model['file_end_piC'] + '_zon1D.nc'

    f2dCO2 = open_ncfile(indir_1pctCO2 + file2d_1pctCO2,'r')
    f1dCO2 = open_ncfile(indir_1pctCO2 + file1d_1pctCO2,'r')
    f2dpiC = open_ncfile(indir_piC + file2d_piC,'r')
    f1dpiC = open_ncfile(indir_piC + file1d_piC,'r')


    # ----- Work ------

    varname = defVarmme('salinity'); v = 'S'
    #varname = defVarmme('temp'); v = 'T'
    #varname = defVarmme('depth'); v = 'Z'

    multStd = 2. # detect ToE at multStd std dev of piControl

    labBowl = ['piControl', '2*CO2']

    valmask = 1.e20

    # Years for bowl : average around year 70 when CO2 has doubled
    y1 = 69
    y2 = 74

    # ----- Variables ------

    lat = f2dCO2.variables['latitude'][:]; latN = lat.size
    density = f2dCO2.variables['lev'][:]; levN = density.size
    time = f2dCO2.variables['time'][:]; timN = time.size

    var = varname['var_zonal']

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

    # -- Read var 1pctCO2
    varCO2_a = f2dCO2.variables[var][:,1,:,:].squeeze()
    varCO2_p = f2dCO2.variables[var][:,2,:,:].squeeze()
    varCO2_i = f2dCO2.variables[var][:,3,:,:].squeeze()
    # -- Read var piControl
    if name == 'mme':
        varpiC_a = f2dpiC.variables[var][:,1,:,:].squeeze()
        varpiC_p = f2dpiC.variables[var][:,2,:,:].squeeze()
        varpiC_i = f2dpiC.variables[var][:,3,:,:].squeeze()
    else :
        varpiC_a = f2dpiC.variables[var][-140:,1,:,:].squeeze()
        varpiC_p = f2dpiC.variables[var][-140:,2,:,:].squeeze()
        varpiC_i = f2dpiC.variables[var][-140:,3,:,:].squeeze()

    # -- Read lightest density of persistent ocean (ptopsigma)
    bowlCO2_a = f1dCO2.variables['ptopsigma'][y1:y2,1,:].squeeze()
    bowlCO2_p = f1dCO2.variables['ptopsigma'][y1:y2,2,:].squeeze()
    bowlCO2_i = f1dCO2.variables['ptopsigma'][y1:y2,3,:].squeeze()
    bowlpiC_a = f1dpiC.variables['ptopsigma'][y1:y2,1,:].squeeze()
    bowlpiC_p = f1dpiC.variables['ptopsigma'][y1:y2,2,:].squeeze()
    bowlpiC_i = f1dpiC.variables['ptopsigma'][y1:y2,3,:].squeeze()

    # -- Check mean fields
    varCO2mean_a = np.ma.average(varCO2_a[0:5,:,:], axis=0)
    varCO2mean_p = np.ma.average(varCO2_p[0:5,:,:], axis=0)
    varCO2mean_i = np.ma.average(varCO2_i[0:5,:,:], axis=0)

    # ----- Build plot variables ------

    # -- Compute significance of difference when diff within 1 stddev of piControl variability (in the MME sense)
    varams = np.ma.std(varpiC_a, axis=0)
    varpms = np.ma.std(varpiC_p, axis=0)
    varims = np.ma.std(varpiC_i, axis=0)

    # -- reorganise i,j dims in single dimension data (speeds up loops)
    varCO2_a  = np.reshape(varCO2_a, (timN,levN*latN))
    varpiC_a = np.reshape(varpiC_a,(timN,levN*latN))
    varams  = np.reshape(varams, (levN*latN))
    varCO2_p  = np.reshape(varCO2_p, (timN,levN*latN))
    varpiC_p = np.reshape(varpiC_p,(timN,levN*latN))
    varpms  = np.reshape(varpms, (levN*latN))
    varCO2_i  = np.reshape(varCO2_i, (timN,levN*latN))
    varpiC_i = np.reshape(varpiC_i,(timN,levN*latN))
    varims  = np.reshape(varims, (levN*latN))

    # -- Compute ToE as last date when diff 1pctCO2 - piControl is larger than mult * stddev
    toe_a = np.reshape(findToE(varCO2_a-varpiC_a, varams, multStd),(levN,latN))
    toe_p = np.reshape(findToE(varCO2_p-varpiC_p, varpms, multStd),(levN,latN))
    toe_i = np.reshape(findToE(varCO2_i-varpiC_i, varims, multStd),(levN,latN))

    # -- Average bowl position
    bowlCO2_a = np.ma.average(bowlCO2_a, axis=0)
    bowlCO2_p = np.ma.average(bowlCO2_p, axis=0)
    bowlCO2_i = np.ma.average(bowlCO2_i, axis=0)
    bowlpiC_a = np.ma.average(bowlpiC_a, axis=0)
    bowlpiC_p = np.ma.average(bowlpiC_p, axis=0)
    bowlpiC_i = np.ma.average(bowlpiC_i, axis=0)

    # -- Mask
    var_mask = np.ma.getmask(np.ma.average(f2dCO2.variables[var][:], axis=0))
    toe_a = np.ma.array(toe_a, mask=var_mask[1,:,:])
    toe_p = np.ma.array(toe_p, mask=var_mask[2,:,:])
    toe_i = np.ma.array(toe_i, mask=var_mask[3,:,:])

    # -- Create variable bundles
    varAtl = {'name': 'Atlantic', 'ToE': toe_a, 'bowl2': bowlCO2_a, 'bowl1': bowlpiC_a, 'labBowl': labBowl,
              'var_mean': varCO2mean_a, 'bowl': bowlCO2_a}
    varPac = {'name': 'Pacific', 'ToE': toe_p, 'bowl2': bowlCO2_p, 'bowl1': bowlpiC_p, 'labBowl': labBowl,
              'var_mean': varCO2mean_p, 'bowl': bowlCO2_p}
    varInd = {'name': 'Indian', 'ToE': toe_i, 'bowl2': bowlCO2_i, 'bowl1': bowlpiC_i, 'labBowl': labBowl,
              'var_mean': varCO2mean_i, 'bowl': bowlCO2_i}


    # ------------------------------------
    #               Plot
    # ------------------------------------

    # # == Check mean fields ==
    #
    # fig0, axes = plt.subplots(nrows=2, ncols=3, figsize=(17,5))
    #
    # cmap = plt.get_cmap('jet')
    # levels = None
    #
    # cnplot = zonal_2D(plt, 'var_mean', axes[0, 0], axes[1, 0], 'left', lat, density, varAtl, domrho, cmap, levels, clevsm, clevsm_bold)
    #
    # cnplot = zonal_2D(plt, 'var_mean', axes[0, 1], axes[1, 1], 'mid', lat, density, varPac, domrho, cmap, levels, clevsm, clevsm_bold)
    #
    # cnplot = zonal_2D(plt, 'var_mean', axes[0, 2], axes[1, 2], 'right', lat, density, varInd, domrho, cmap, levels, clevsm, clevsm_bold)
    #
    #
    # plt.subplots_adjust(hspace=.0001, wspace=0.05, left=0.04, right=0.86)
    #
    # cb = plt.colorbar(cnplot, ax=axes.ravel().tolist(), fraction=0.015, shrink=2.0, pad=0.05)
    # cb.set_label('%s (%s)' % (legVar,unit), fontweight='bold')
    #
    # if name == 'ensemble_mean':
    #     name = model['name']
    # plotTitle = 'Mean ' + legVar + ' field (' + name + ' 1pctCO2, first 5 years)'
    #
    # plt.suptitle(plotTitle, fontweight='bold', fontsize=14, verticalalignment='top')


    # == Plot ToE ==
    fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(17,5))

    minmax = [0, 141, 10]
    unit = 'ToE'
    cmap = cmaps.viridis
    levels = np.arange(minmax[0], minmax[1], minmax[2])

    cnplot = zonal_2D(plt, 'ToE', axes[0, 0], axes[1, 0], 'left', lat, density, varAtl, domrho, cmap, levels, clevsm, clevsm_bold)

    cnplot = zonal_2D(plt, 'ToE', axes[0, 1], axes[1, 1], 'mid', lat, density, varPac, domrho, cmap, levels, clevsm, clevsm_bold)

    cnplot = zonal_2D(plt, 'ToE', axes[0, 2], axes[1, 2], 'right', lat, density, varInd, domrho, cmap, levels, clevsm, clevsm_bold)


    plt.subplots_adjust(hspace=.0001, wspace=0.05, left=0.04, right=0.86)

    cb = plt.colorbar(cnplot, ax=axes.ravel().tolist(), ticks = levels[::2], fraction=0.015, shrink=2.0, pad=0.05)
    cb.set_label('%s' % (unit,), fontweight='bold')

    if name == 'single_model':
        name = model['name']
    plotTitle = name + ' ' + legVar + ' ToE 1pctCO2 vs. piControl [> ' + str(multStd) + ' std]'
    plotName = name + '_ToE_' + legVar + '_1pctCO2vsPiC'

    plt.suptitle(plotTitle, fontweight='bold', fontsize=14, verticalalignment='top')

    plt.figtext(.5,.02,'Computed by : zonal_1pctCO2vsPiC_toe.py',fontsize=9,ha='center')

    #plt.show()

    plt.savefig('/home/ysilvy/Density_bining/Yona_analysis/figures/models/zonal_ys/ToE/'+plotName+'.png', bbox_inches='tight')

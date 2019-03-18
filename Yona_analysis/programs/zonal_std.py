#!/bin/env python
# -*- coding: utf-8 -*-

"""
Python matplotlib
Plot density/latitude standard deviation of salinity for each basin, each model
Choose between :
- HistoricalNat (ensemble mean, max std)
- PiControl

"""

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset as open_ncfile
from maps_matplot_lib import zonal_2D, defVarmme, custom_div_cmap
import glob, os
import datetime
import colormaps as cmaps

# ===== Define Work ======

# work  = 'histNat'
work = 'piC'

indir_histNat = '/data/ericglod/Density_binning/Prod_density_april15/mme_histNat/'
indir_piC = '/data/ericglod/Density_binning/Prod_density_april15/mme_piControl/'
fig_dir = '/home/ysilvy/Density_bining/Yona_analysis/figures/models/zonal_ys/'

varname = defVarmme('salinity'); v = 'S'
# varname = defVarmme('temp'); v = 'T'

var = varname['var_zonal_w/bowl']
legVar = varname['legVar']
unit = varname['unit']

# Density domain
rhomin = 21
rhomid = 26
rhomax = 28
domrho = [rhomin, rhomid, rhomax]

if work == 'histNat':
    indir = indir_histNat
    longName = 'historicalNat'
    dir = 'histNat_std/'
else:
    indir = indir_piC
    longName = 'piControl'
    dir = 'piControl_std/'

# ==== Read std and plot, for each model in list ====

# -- Read files in directory
listfiles = sorted(glob.glob(indir + '/*2D.nc'))
nmodels = len(listfiles)

# -- Loop on models
for i in range(nmodels):
# i=18
    file = os.path.basename(listfiles[i])
    f = open_ncfile(indir+file,'r')
    name = file.split('.')[1] # Read model name

    print('Reading '+work+' for '+name)

    if work == 'histNat':
        if name != 'multimodel_Nat':
            # Read varstd of histNat (max std of all runs for each model)
            stdvar_a = f.variables[var+'Std'][1,:,:].squeeze()
            stdvar_i = f.variables[var+'Std'][3,:,:].squeeze()
            stdvar_p = f.variables[var+'Std'][2,:,:].squeeze()
        else:
            # Read var of histNat for the multimodel
            varmme = f.variables[var][:]
            stdvar_a = np.ma.std(varmme[:,1,:,:],axis=0)
            stdvar_p = np.ma.std(varmme[:,2,:,:],axis=0)
            stdvar_i = np.ma.std(varmme[:,3,:,:],axis=0)

    else:
        # Read PiControl over 240 years + compute std of PiControl
        varpiC = f.variables[var][-240:,:,:,:]
        stdvar_a = np.ma.std(varpiC[:,1,:,:], axis=0)
        stdvar_p = np.ma.std(varpiC[:,2,:,:], axis=0)
        stdvar_i = np.ma.std(varpiC[:,3,:,:], axis=0)

    density = f.variables['lev'][:]
    lat = f.variables['latitude'][:]

    # Read bowl and take the average
    file2 = glob.glob(indir+'/*'+name+'*1D.nc')[0]
    f2 = open_ncfile(file2,'r')
    if work == 'histNat':
        bowl = f2.variables['ptopsigma'][:]
    else:
        bowl = f2.variables['ptopsigma'][-240:,:,:]
    bowl = np.ma.average(bowl,axis=0)

    # == Plot ==

    # Create variable bundles
    varPac = {'name': 'Pacific', 'var_std': stdvar_p, 'bowl':bowl[2,:]}
    varAtl = {'name': 'Atlantic','var_std': stdvar_a, 'bowl':bowl[1,:]}
    varInd = {'name': 'Indian', 'var_std': stdvar_i, 'bowl':bowl[3,:]}

    fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(17,5))

    cmap = cmaps.viridis
    levels = np.arange(0,0.1201,0.01)

    cnplot = zonal_2D(plt, 'var_std', axes[0,0], axes[1,0], 'left', lat, density, varAtl, domrho, cmap, levels)
    cnplot = zonal_2D(plt, 'var_std', axes[0,1], axes[1,1], 'mid', lat, density, varPac, domrho, cmap, levels)
    cnplot = zonal_2D(plt, 'var_std', axes[0,2], axes[1,2], 'right', lat, density, varInd, domrho, cmap, levels)


    plt.subplots_adjust(hspace=.0001, wspace=0.05, left=0.04, right=0.86)

    cb = plt.colorbar(cnplot[0], ax=axes.ravel().tolist(), ticks=levels[::2], fraction=0.015, shrink=2.0, pad=0.05)
    cb.set_label('Standard deviation (%s)' % (unit,), fontweight='bold')

    plotTitle = 'Standard deviation of '+legVar+', '+name+' '+longName
    plotName = name+'_'+legVar+'_std_'+work

    # Date
    now = datetime.datetime.now()
    date = now.strftime("%Y-%m-%d")

    plt.suptitle(plotTitle, fontweight='bold', fontsize=14, verticalalignment='top')
    plt.figtext(.5,.015,'Computed by : zonal_std.py  '+date,fontsize=8,ha='center')

    plt.savefig(fig_dir+dir+plotName+'.png', bbox_inches='tight')

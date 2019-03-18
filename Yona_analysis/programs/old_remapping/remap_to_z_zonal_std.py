#!/bin/env python
# -*- coding: utf-8 -*-

"""
Python matplotlib
Plot pseudo-depth/latitude standard deviation of salinity for each basin, each model
(read from saved file with remapped field)
Choose between :
- HistoricalNat (ensemble mean, max std)
- PiControl

"""

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset as open_ncfile
from maps_matplot_lib import defVarmme, remapToZ, zon_2Dz
import glob, os
import datetime
import colormaps as cmaps

# ===== Define Work ======

work  = 'histNat'
# work = 'piC'

indir_histNat = '/data/ericglod/Density_binning/Prod_density_april15/mme_histNat/'
indir_piC = '/data/ericglod/Density_binning/Prod_density_april15/mme_piControl/'
indir_histNat_remap = '/home/ysilvy/Density_bining/Yona_analysis/data/remaptoz/histNat_std/'
indir_piC_remap = '/home/ysilvy/Density_bining/Yona_analysis/data/remaptoz/piControl_std/'
fig_dir = '/home/ysilvy/Density_bining/Yona_analysis/figures/models/zonal_remaptoz/'

varname = defVarmme('salinity'); v = 'S'
# varname = defVarmme('temp'); v = 'T'

var = varname['var_zonal_w/bowl']
legVar = varname['legVar']
unit = varname['unit']

if work == 'histNat':
    indir_z = indir_histNat_remap
    indir = indir_histNat
    longName = 'historicalNat'
    dir = 'histNat_std/'
else:
    indir_z = indir_piC_remap
    indir = indir_piC
    longName = 'piControl'
    dir = 'piControl_std/'

# ==== Read remapped std from file ====

listfiles = sorted(glob.glob(indir_z+'*.nc'))
# -- Loop on models
for i in range(len(listfiles)):
# i=0
    file = os.path.basename(listfiles[i])
    f = open_ncfile(indir_z+file,'r')
    name = file.split('.')[1] # Read model name

    print('Reading '+work+' for '+name)

    stdvar_z = f.variables[var+'Std'][:]
    lat = f.variables['latitude'][:]
    pseudo_depth = f.variables['pseudo_depth'][:]

    # Read bowl and take the average
    file2 = glob.glob(indir+'/*'+name+'*1D.nc')[0]
    f2 = open_ncfile(file2,'r')
    if work == 'histNat':
        bowl = f2.variables['ptopdepth'][:]
    else:
        bowl = f2.variables['ptopdepth'][-240:,:,:]
    bowl = np.ma.average(bowl,axis=0)

    # ==== Plot ====

    # Create variable bundles
    varAtl = {'name': 'Atlantic','var_change': stdvar_z[1,:,:], 'bowl1':bowl[1,:], 'bowl2':None, 'labBowl': ['bowl',None], 'density':None}
    varPac = {'name': 'Pacific', 'var_change': stdvar_z[2,:,:], 'bowl1':bowl[2,:], 'bowl2':None, 'labBowl': ['bowl',None], 'density':None}
    varInd = {'name': 'Indian', 'var_change': stdvar_z[3,:,:], 'bowl1':bowl[3,:], 'bowl2':None, 'labBowl': ['bowl',None], 'density':None}

    domzed = [0,500,5000]

    fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(17, 5))

    cmap = cmaps.viridis
    levels = np.arange(0,0.1201,0.01)
    ext_cmap = 'max'

    contourDict = {'cmap':cmap, 'levels':levels, 'ext_cmap':ext_cmap}

    cnplot = zon_2Dz(plt, axes[0,0], axes[1,0], 'left', lat, pseudo_depth, varAtl,
                 contourDict, domzed)
    cnplot = zon_2Dz(plt, axes[0,1], axes[1,1], 'mid', lat, pseudo_depth, varPac,
                     contourDict, domzed)
    cnplot = zon_2Dz(plt, axes[0,2], axes[1,2], 'right', lat, pseudo_depth, varInd,
                     contourDict, domzed)

    plt.subplots_adjust(hspace=.0001, wspace=0.05, left=0.04, right=0.86)

    # Add colorbar
    cb = plt.colorbar(cnplot[0], ax=axes.ravel().tolist(), ticks=levels[::2], fraction=0.015, shrink=2.0, pad=0.05)
    cb.set_label('Standard deviation (%s)' % (unit,), fontweight='bold')

    plotTitle = 'Standard deviation of '+legVar+', '+name+' '+longName
    plotName = name+'_'+legVar+'_std_'+work+'_remapping'

    # Date
    now = datetime.datetime.now()
    date = now.strftime("%Y-%m-%d")

    plt.suptitle(plotTitle, fontweight='bold', fontsize=14, verticalalignment='top')
    plt.figtext(.5,.015,'Computed by : remap_to_z_zonal_std.py  '+date,fontsize=8,ha='center')
    plt.figtext(.004,.65,'Pseudo-depth (m)',rotation='vertical',fontweight='bold')

    plt.savefig(fig_dir+dir+plotName+'.png', bbox_inches='tight')

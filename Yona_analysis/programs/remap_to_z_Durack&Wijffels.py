#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 24 14:51:58 2018

@author: Yona
"""

from scipy.interpolate import griddata
from matplotlib.ticker import AutoMinorLocator
from maps_matplot_lib import zon_2Dz

#TODO : début du script puis transférer sur Ciclad

pressure = fh2d.variables['pressure_mean_basin_zonal'][:].squeeze()

# WOA13 grid
targetz = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80,
   85, 90, 95, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 375,
   400, 425, 450, 475, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950,
   1000, 1050, 1100, 1150, 1200, 1250, 1300, 1350, 1400, 1450, 1500, 1550,
   1600, 1650, 1700, 1750, 1800, 1850, 1900, 1950, 2000, 2100, 2200, 2300,
   2400, 2500, 2600, 2700, 2800, 2900, 3000, 3100, 3200, 3300, 3400, 3500,
   3600, 3700, 3800, 3900, 4000, 4100, 4200, 4300, 4400, 4500, 4600, 4700,
   4800, 4900, 5000, 5100, 5200, 5300, 5400, 5500]

var_change_z = np.ma.masked_all((len(targetz),141,4))

for ibasin in range(4):
    for ilat in range(len(lat)):
        Psig = pressure[:,ilat,ibasin]
        var_change_sig = var_change[:,ilat,ibasin]
        Psort = np.ma.compressed(np.sort(Psig))
        isort = np.argsort(Psig)
        var_change_sort = np.ma.compressed(var_change_sig[isort])
        if len(Psort) > 1:
            var_change_z[:,ilat,ibasin] = griddata(Psort,var_change_sort,targetz)
        else :
            var_change_z[:,ilat,ibasin] = np.ma.masked
            

domzed = [0,500,2000]

# -- Make variable bundles for each basin
varAtl = {'name': 'Atlantic', 'var_change': var_change_z[:,:,2], 'bowl1': None, 'bowl2': None, 'labBowl': None} #Rajouter mean et err si code okay
varPac = {'name': 'Pacific', 'var_change': var_change_z[:,:,1], 'bowl1': None, 'bowl2': None, 'labBowl': None}
varInd = {'name': 'Indian', 'var_change': var_change_z[:,:,3], 'bowl1': None, 'bowl2': None, 'labBowl': None}

# == PLOT ==

# -- Create figure and axes instances
fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(17, 5))

# -- color map
cmap = custom_div_cmap()

# -- levels
levels = np.linspace(minmax[0],minmax[1],minmax[2])
    
# -- Contourf of signal
cnplot = zon_2Dz(plt,axes[0,0], axes[1,0], 'left', lat, targetz, varAtl,
                 cmap, levels, domzed)
cnplot = zon_2Dz(plt,axes[0,1], axes[1,1], 'mid', lat, targetz, varPac,
                 cmap, levels, domzed)
cnplot = zon_2Dz(plt,axes[0,2], axes[1,2], 'right', lat, targetz, varInd,
                 cmap, levels, domzed)


plt.subplots_adjust(hspace=.0001, wspace=0.05, left=0.04, right=0.86)

# -- Add colorbar
cb = plt.colorbar(cnplot[0], ax=axes.ravel().tolist(), ticks=levels[::3], fraction=0.015, shrink=2.0, pad=0.05)
cb.set_label('%s (%s)' % (legVar, unit), fontweight='bold')


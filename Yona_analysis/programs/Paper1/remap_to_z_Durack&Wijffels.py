#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 24 14:51:58 2018

@author: Yona Silvy

Data : Durack & Wijffels 1950-2000 zonal changes, analyzed in density framework
Grid back to a pressure coordinate and plot
"""
import sys
#sys.path.append('/home/ysilvy/Density_bining/Yona_analysis/programs/')
sys.path.append('/Users/Yona 1/Documents/Thèse/Density_bining/Yona_analysis/programs/') # From local
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset as open_ncfile
from scipy.interpolate import griddata
from maps_matplot_lib import defVarDurack, zon_2Dz, custom_div_cmap
import datetime
import gsw

# ===== Workspace =====

#indir = '/home/ysilvy/Density_bining/Yona_analysis/data/'
indir = '/Users/Yona 1/Documents/Thèse/Data/' #Local
file = 'DurackandWijffels_GlobalOceanChanges-NeutralDensity_1950-2000_170224_20_48_22_beta.nc'
data = indir + file
f2d = open_ncfile(data, 'r')

out = 'save' # View or save output figure
#out = 'view'

show_isopyc = True #False

# ===== Read variables =====

lat = f2d.variables['latitude'][:]
density = f2d.variables['density'][:]
pressure = f2d.variables['pressure_mean_basin_zonal'][:].squeeze()

# -- Choose wich variable to work on
varname = defVarDurack('salinity'); v = 'S'
# varname = defVarDurack('temp'); v = 'T'

var_mean = varname['var_mean_zonal']
var_change = varname['var_change_zonal']
var_change_er = varname['var_change_zonal_er']

var_attributes = f2d.variables[var_mean]
var_mean = f2d.variables[var_mean][:].squeeze()
var_change = f2d.variables[var_change][:].squeeze()
var_change_er = f2d.variables[var_change_er][:].squeeze()

# -- Determine field in 1950 and 2000 from mean field
var_1950 = var_mean - var_change/2
var_2000 = var_mean + var_change/2

# -- Define variable properties
minmax = varname['minmax_zonal']
clevsm = varname['clevsm_zonal']
clevsm_bold = varname['clevsm_bold']
legVar = varname['legVar']
unit = varname['unit']

# -- Density domain
rhomin = 21
rhomid = 26
rhomax = 28
domrho = [rhomin, rhomid, rhomax]

valmask = 1e+20

# ===== REMAP ======

# -- First convert pressure to depth using gsw package
depth_a = -gsw.z_from_p(pressure[:,:,2],lat)
depth_p = -gsw.z_from_p(pressure[:,:,1],lat)
depth_i = -gsw.z_from_p(pressure[:,:,3],lat)

depth = np.ma.masked_all_like(pressure)
depth[:,:,2] = depth_a
depth[:,:,1] = depth_p
depth[:,:,3] = depth_i

# -- WOA13 grid for target grid
targetz = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80,
   85, 90, 95, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 375,
   400, 425, 450, 475, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950,
   1000, 1050, 1100, 1150, 1200, 1250, 1300, 1350, 1400, 1450, 1500, 1550,
   1600, 1650, 1700, 1750, 1800, 1850, 1900, 1950, 2000, 2100, 2200, 2300,
   2400, 2500, 2600, 2700, 2800, 2900, 3000, 3100, 3200, 3300, 3400, 3500,
   3600, 3700, 3800, 3900, 4000, 4100, 4200, 4300, 4400, 4500, 4600, 4700,
   4800, 4900, 5000, 5100, 5200, 5300, 5400, 5500]

# -- Meshgrid density for remapping isopycnals

lat2d, density2d = np.meshgrid(lat,density)

# -- Initialize fields in pseudo-z
var_change_z = np.ma.masked_all((len(targetz),141,4))
var_mean_z = np.ma.masked_all((len(targetz),141,4))
var_change_er_z = np.ma.masked_all((len(targetz),141,4))
density_z = np.ma.masked_all((len(targetz),141,4))

# -- Grid fields to target depth grid, using the depth field in density/lat framework
for ibasin in range(4):
    for ilat in range(len(lat)):

        zsig = depth[:,ilat,ibasin]
        var_change_sig = var_change[:,ilat,ibasin]
        var_mean_sig = var_mean[:,ilat,ibasin]
        var_change_er_sig = var_change_er[:,ilat,ibasin]
        density_sig = density2d[:,ilat]

        zsort = np.ma.compressed(np.sort(zsig)) # Sort depth and get rid of nans
        isort = np.argsort(zsig)
        var_change_sort = np.ma.compressed(var_change_sig[isort]) # Change order according to corresponding depth vector
        var_mean_sort = np.ma.compressed(var_mean_sig[isort])
        var_change_er_sort = np.ma.compressed(var_change_er_sig[isort])
        density_sort = np.ma.compressed(density_sig[isort])

        if len(zsort) > 1:
            var_change_z[:,ilat,ibasin] = griddata(zsort,var_change_sort,targetz) # Grid field with target depth grid
            var_mean_z[:,ilat,ibasin] = griddata(zsort,var_mean_sort,targetz)
            var_change_er_z[:,ilat,ibasin] = griddata(zsort,var_change_er_sort,targetz)
            density_z[:,ilat,ibasin] = griddata(zsort,density_sort,targetz)
        else :
            var_change_z[:,ilat,ibasin] = np.ma.masked
            var_mean_z[:,ilat,ibasin] = np.ma.masked
            var_change_er_z[:,ilat,ibasin] = np.ma.masked
            density_z[:,ilat,ibasin] = np.ma.masked

# -- Depth domain for the plot
domzed = [0,500,5000]

# -- Make variable bundles for each basin
varAtl = {'name': 'Atlantic', 'var_change': var_change_z[:,:,2], 'var_mean':var_mean_z[:,:,2],
          'bowl1': None, 'bowl2': None, 'labBowl': None, 'density':density_z[:,:,2]}
varPac = {'name': 'Pacific', 'var_change': var_change_z[:,:,1], 'var_mean':var_mean_z[:,:,1],
          'bowl1': None, 'bowl2': None, 'labBowl': None, 'density':density_z[:,:,1]}
varInd = {'name': 'Indian', 'var_change': var_change_z[:,:,3], 'var_mean':var_mean_z[:,:,3],
          'bowl1': None, 'bowl2': None, 'labBowl': None, 'density':density_z[:,:,3]}

# ===== PLOT ======

# -- Create figure and axes instances
fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(17, 5))

# -- color map
cmap = custom_div_cmap()

# -- levels
levels = np.linspace(minmax[0],minmax[1],minmax[2])

# Dictionary
ext_cmap = 'both'
contourDict = {'cmap':cmap, 'levels':levels, 'levels2':levels, 'ext_cmap':ext_cmap, 'isopyc':show_isopyc}

# -- Error field
var_change_er_z = var_change_er_z * 1.1  # to account for a potential underestimation of the error determined by a bootstrap analysis
var_change_er_z = var_change_er_z * 1.64  # 90% confidence level
not_signif_change = np.ma.where(np.absolute(var_change_z) < var_change_er_z, 1, 0)

# -- Contourf of signal
cnplot = zon_2Dz(plt,axes[0,0], axes[1,0], 'left', lat, targetz, varAtl,
                 contourDict, domzed)
cnplot = zon_2Dz(plt,axes[0,1], axes[1,1], 'mid', lat, targetz, varPac,
                 contourDict, domzed)
cnplot = zon_2Dz(plt,axes[0,2], axes[1,2], 'right', lat, targetz, varInd,
                 contourDict, domzed)

# -- Draw areas where signal is not significant
lat2d, lev2d = np.meshgrid(lat, targetz)
for i in range(2):
    axes[i,0].contourf(lat2d, lev2d, not_signif_change[:,:,2], levels=[0.25,0.5,1.5], colors='None',
                       hatches=['','....'])
    axes[i,1].contourf(lat2d, lev2d, not_signif_change[:,:,1], levels=[0.25,0.5,1.5], colors='None',
                       hatches=['','....'])
    axes[i,2].contourf(lat2d, lev2d, not_signif_change[:,:,3], levels=[0.25,0.5,1.5], colors='None',
                       hatches=['','....'])

axes[0,1].tick_params(axis='both',labelleft=False,which='both',top=False,bottom=False)
axes[1,1].tick_params(axis='both',labelleft=False,which='both',top=False)
axes[0,2].tick_params(axis='both',labelleft=False,which='both',top=False,bottom=False)
axes[1,2].tick_params(axis='both',labelleft=False,which='both',top=False)
axes[0,0].tick_params(axis='x',which='both',top=False,bottom=False)
axes[1,0].tick_params(axis='x',which='both',top=False)

plt.subplots_adjust(hspace=.012, wspace=0.05, left=0.04, right=0.86)
plt.figtext(.005,.38,'Pseudo-depth (m)',rotation='vertical',fontweight='bold',ha='center')

plt.figtext(.006,.96,'a',fontweight='bold',fontsize=16)

# -- Add colorbar
cb = plt.colorbar(cnplot[0], ax=axes.ravel().tolist(), ticks=levels[::3], fraction=0.015, shrink=2.0, pad=0.05)
cb.set_label('%s (%s)' % (legVar, unit), fontweight='bold')


name = 'Durack & Wijffels (2010)'
plotTitle = '%s changes 1950-2000, %s' %(legVar, name)
plotName = 'remapping_depth_Durack&Wijffels_'+ v + 'changes_isopyccontours_5000_paper'

# Date
now = datetime.datetime.now()
date = now.strftime("%Y-%m-%d")

#plt.suptitle(plotTitle, fontweight='bold', fontsize=14, verticalalignment='top')
#axes[0,1].set_title(plotTitle, y=1.25, fontweight='bold', fontsize=15, verticalalignment='top')
#plt.figtext(.5,.01,'Computed by : remap_to_z_Durack&Wijffels.py  '+date,fontsize=8,ha='center')

figureDir = 'obs/zonal_remaptoz/'

if out == 'view':
    plt.show()
else:
    #plt.savefig('/home/ysilvy/figures/'+figureDir+plotName+'.png', bbox_inches='tight') 
    plt.savefig('/Users/Yona 1/Documents/Thèse/Yona_analysis/figures/'+figureDir+plotName+'.png', bbox_inches='tight',dpi=150) 

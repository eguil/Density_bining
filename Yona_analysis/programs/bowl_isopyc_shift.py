#!/bin/env python
# -*- coding: utf-8 -*-

"""
Python matplotlib
Make a map of the shift in isopycnals on the bowl between 1945 and 2010

"""

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset as open_ncfile
from maps_matplot_lib import custom_div_cmap
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1 import make_axes_locatable

# ----- Workspace ------

#name = 'Ishii'
name = 'mme'

if name == 'Ishii' :
    indir = '/data/ericglod/Density_binning/Obs_Prod_density_april16/'
    file = 'obs.Ishii.historical.r0i0p0.an.ocn.Omon.density.ver-1.latestXCorr.nc'
    ### Ishii : 1945 to 2013
else :
    indir = '/data/ericglod/Density_binning/Prod_density_april15/Raw/mme_hist/'
    file = 'cmip5.multimodel_Nat.historical.ensm.an.ocn.Omon.density_2D.nc'
    ### mme : 1861 to 2006

data = indir + file

f = open_ncfile(data,'r')

# Read grid for coloring continents
f2 = open_ncfile('/home/ysilvy/Density_bining/Yona_analysis/data/140807_WOD13_masks.nc', 'r')


# ----- Variables ------

# Read variables
lat = f.variables['latitude'][:]
lon = f.variables['longitude'][:]
ptopsigma = f.variables['ptopsigmaxy'][:]
ptopsalinity = f.variables['ptopsoxy'][:]


landsea = f2.variables['landsea'][:]


# ----- Build plot variables -----

mean_salinity = np.ma.average(ptopsalinity, axis=0)
if name == 'Ishii' :
    density_begin = np.ma.average(ptopsigma[0:5,:,:], axis=0)
    density_end = np.ma.average(ptopsigma[-5:,:,:], axis=0)
if name == 'mme' :
    density_begin = np.ma.average(ptopsigma[84:89, :, :], axis=0)
    density_end = np.ma.average(ptopsigma[-5:, :, :], axis=0)



lon2d, lat2d = np.meshgrid(lon, lat)

# Create mask for continents
sea_mask = landsea != 1
landsea = np.ma.array(landsea, mask=sea_mask)
landsea[landsea == 1] = 0.2


# ----- Plot density changes on a map -----

fig = plt.figure(figsize=(13,5))

# Basemap
map = Basemap(projection='cyl', llcrnrlon=20, llcrnrlat=-70., urcrnrlon=380, urcrnrlat=70)
map.drawmapboundary(fill_color='1')
map.drawparallels(np.arange(-60, 61, 20.), labels=[1, 0, 1, 0], linewidth=0.5)
map.drawmeridians(np.arange(-180, 180, 60), labels=[0, 0, 0, 1], linewidth=0.5)

# Draw continents
pc2 = map.pcolormesh(lon2d, lat2d, landsea, shading='flat', cmap=plt.cm.gray, latlon=True)

# Draw filled contour of mean salinity field
levels = np.arange(33,37.5,0.5)
cnplot = map.contourf(lon2d, lat2d, mean_salinity, levels=levels,
                      cmap=plt.get_cmap('RdYlBu_r'), latlon=True, extend='both')

# Draw 1945 isopycnals
plot1945 = map.contour(lon2d, lat2d, density_begin, levels=np.arange(21, 29, 1), colors='black',
            linewidths=2, latlon=True)
plt.clabel(plot1945, inline=1, fontsize=10, fmt='%d', fontweight='bold')
# Draw 21st century isopycnals
plot2013 = map.contour(lon2d, lat2d, density_end, levels=np.arange(22, 29, 1), colors='white',
            linewidths=2, latlon=True)
#plt.clabel(plot2013, plot2013.levels[1::2], inline=1, fontsize=10, fmt='%d', fontweight='bold')

ax = plt.gca()
ax.set_title('Isopycnal shift on the bowl (%s)' % (name,),
          fontsize=16, va='bottom')

# Make colorbar same height as map
ax = plt.gca()
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="3%", pad=0.15)
cb = plt.colorbar(cnplot, cax=cax)


plt.show()

#!/bin/env python
# -*- coding: utf-8 -*-
"""
Python matplotlib 
Make density/latitude section for Atl/Pac/Ind for a number of variables

(c) Eric Guilyardi Feb 2016

"""
import numpy as np
from   netCDF4 import Dataset as open_ncfile
import matplotlib as mpl
from   mpl_toolkits.basemap import Basemap, cm
from   mpl_toolkits.axes_grid1 import Grid
import matplotlib.pyplot as plt
from   matplotlib.colors import LinearSegmentedColormap

from   densitlib import zon_2dom, gmtColormap

# -------------------------------------------------------------------------------
#                               Define work
# -------------------------------------------------------------------------------

indir = "/Users/ericg/Projets/Density_bining/Prod_density_april15/mme_hist/r1i1p1"
file  = "cmip5.multimodel.historical.ensm.an.ocn.Omon.density_zon2D.nc"

salinity = {
    'var'    : 'isonso',
    'minmax' : [-.2,.2],            # for diff shading
    'clevsm' : np.arange(30,40,.2), # for mean contours
    'legVar' : "Salinity",
    'unit'   : "PSU",               # TODO: could be read from file
    }

temp = {
    'var'    : 'isonthetao',
    'minmax' : [-.4,.4],
    'clevsm' : np.arange(-2,30,1),
    'legVar' : "Temperature",
    'unit'   : "C" ,
    }

# Define variable
varname = salinity
varname = temp

# years for diff
y11 = 140-1
y12 = 146-1
y21 = 2-1
y22 = 30-1

# density domain
domrho = [21.,26.,28.]               # min/mid/max

#
# -------------------------------------------------------------------------------

# Define variable properties

var    = varname['var']
minmax = varname['minmax']
clevsm = varname['clevsm']
legVar = varname['legVar']
unit = varname['unit']

#-- open netcdf file
nc = open_ncfile(indir+'/'+file)

#-- read variable
tvara = nc.variables[var+'Bowl'][:,1,:,:].squeeze()
tvarp = nc.variables[var+'Bowl'][:,2,:,:].squeeze()
tvari = nc.variables[var+'Bowl'][:,3,:,:].squeeze()
lev = nc.variables['lev'][:]
lat = nc.variables['latitude'][:]

#-- Build plot variables
# difference
vara = np.ma.average(tvara[y11:y12], axis=0)-np.ma.average(tvara[y21:y22], axis=0)
varp = np.ma.average(tvarp[y11:y12], axis=0)-np.ma.average(tvarp[y21:y22], axis=0)
vari = np.ma.average(tvari[y11:y12], axis=0)-np.ma.average(tvari[y21:y22], axis=0)
# mean
varam = np.ma.average(tvara, axis=0)
varap = np.ma.average(tvarp, axis=0)
varai = np.ma.average(tvari, axis=0)

#-- create figure and axes instances
fig, axes = plt.subplots(nrows=2,ncols=3,figsize=(17,5))

#-- color map
#fileGMT='palet_diff_6'
#cdict = gmtColormap(fileGMT,GMTPath = '/Users/ericg/POST_IT/config/palettes')
#cmap = LinearSegmentedColormap('campgmt',cdict)
cmap = plt.get_cmap('bwr') # red/white/blue difference map

#
# -------- Make plot ----------------
#
cnplot=zon_2dom(axes[0,0],axes[1,0],lat,lev,vara,varam,unit,minmax,clevsm,cmap,domrho,legVar+" Atl.",'F')

cnplot=zon_2dom(axes[0,1],axes[1,1],lat,lev,varp,varap,unit,minmax,clevsm,cmap,domrho,legVar+" Pac.",'T')

cnplot=zon_2dom(axes[0,2],axes[1,2],lat,lev,vari,varai,unit,minmax,clevsm,cmap,domrho,legVar+" Ind.",'R')

plt.subplots_adjust(hspace = .00001, wspace=0.05, left=0.04, right=0.86)

#-- Add colorbar
cbar = fig.colorbar(cnplot, ax=axes.ravel().tolist(),fraction=0.015, shrink=2.0,pad=0.03)
cbar.set_label(unit)       

plt.show()
#plt.savefig('test.pdf', bbox_inches='tight')

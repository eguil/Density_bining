#!/bin/env python
# -*- coding: utf-8 -*-

import os, glob
from netCDF4 import Dataset as open_ncfile
import matplotlib.pyplot as plt
from densit_matplot_lib import defVar, remapToZ, zon_2domz
import numpy as np


# -------------------------------------------------------------------------------
#                                Define work
# -------------------------------------------------------------------------------

inDir = '/Users/ericg/Projets/Density_bining/'
work = 'Prod_density_april15/mme_hist'
inDir = inDir + work

# output format
outfmt = 'view'
#outfmt = 'save'

# Define variable  TODO: read as argument
varname = defVar('salinity')
#varname = defVar('temp')
minmax = varname['minmax']
minmax=[30,40,16]
clevsm = varname['clevsm']
legVar = varname['legVar']
unit = varname['unit']

iniyear = 1860
finalyear = 2005

plotName = 'cmip5_remap_test' + varname['var']

# -------------------------------------------------------------------------------
# file inits

os.chdir(inDir)
file='cmip5.multimodel_All.historical.ensm.an.ocn.Omon.density_zon2D.nc'
var = varname['var']

# find dimensions
fi = open_ncfile(inDir+'/'+file)
print inDir+'/'+file

fieldr = fi.variables[var][:]
depthr = fi.variables['isondepth'][:]
volumr = fi.variables['isonvol'][:]
lat = fi.variables['latitude'][:]

valmask = 1.e20

# Target grid

targetz = [0.,5.,15.,25.,35.,45.,55.,65.,75.,85.,95.,105.,116.,128.,142.,158.,181.,216.,272.,364.,511.,732.,1033.,1405.,1830.,2289.,2768.,3257.,3752.,4350.,4749.,5250.]

# Remap
fieldz = remapToZ(fieldr.data,depthr.data,volumr.data, valmask, targetz)
idx = np.argwhere(fieldz < valmask/10.)
print idx.shape

# plot

domzed = [0,1000,5000]

# -- Create figure and axes instances
fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(10, 10))

# -- color map
cmap = plt.get_cmap('bwr')  # red/white/blue difference map

#print lat,targetz, fieldz[0,:]
print len(lat),len(targetz),fieldz.shape
cnplot = zon_2domz(plt, axes[0], axes[1], lat, targetz, fieldz[0,0,:,:], 'Test remap', minmax, clevsm, cmap, domzed, 'F')

plt.subplots_adjust(hspace=.00001, wspace=0.05, left=0.04, right=0.86)

# -- Add colorbar
#cbar = fig.colorbar(cnplot, ax=axes.ravel().tolist(), fraction=0.015, shrink=2.0, pad=0.05)
#cbar.set_label(unit)

# add Title text
#ttxt = fig.suptitle(legVar + ' for ' + work+'/ '+dow, fontsize=14, fontweight='bold')

# -- Output  # TODO read as argument

if outfmt == 'view':
    plt.show()
else:
    print 'Save',plotName+'.pdf'
    plt.savefig(plotName+'.pdf', bbox_inches='tight')

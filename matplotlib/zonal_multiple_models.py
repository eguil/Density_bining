#!/bin/env python
# -*- coding: utf-8 -*-
"""
Python matplotlib 
Make density/latitude section for Atl/Pac/Ind for a number of variables

(c) Eric Guilyardi Feb 2016

TODO: - add arguments for variable and output type
      - read unit from file
      - add mesh/dots for agreement zone

"""
import numpy as np
from   netCDF4 import Dataset as open_ncfile
import matplotlib as mpl
from   mpl_toolkits.basemap import Basemap, cm
from   mpl_toolkits.axes_grid1 import Grid
import matplotlib.pyplot as plt
from   matplotlib.colors import LinearSegmentedColormap

from   densitlib import zon_2dom, defVar

# -------------------------------------------------------------------------------
#                               Define work
# -------------------------------------------------------------------------------

indirHist = '/Users/ericg/Projets/Density_bining/'
workHist = 'Prod_density_april15/historical/r1i1p1'
indirHist = indirHist+workHist
indirHistN = '/Users/ericg/Projets/Density_bining/'
workHistN = 'Prod_density_april15/historicalNat/r1i1p1'
indirHistN = indirHistN+workHistN

# model  name, version Hist, version HistNat

model1 = ['CCSM4',      'v20121128', 'v20121128']
model2 = ['CESM1-CAM5', 'v20130302', 'v20130909']
model3 = ['CNRM-CM5',   'v20130401', 'v20130101']
model1 = ['CESM1-CAM5', 'v20130302', 'v20130909']
model1 = ['CNRM-CM5',   'v20130101', 'v20130101']
#model3 = ['','']
#model3 = ['','']
#model3 = ['','']


file1H2d  = 'cmip5.'+model1[0]+'.historical.r1i1p1.an.ocn.Omon.density.ver-'+model1[1]+'_zon2D.nc'
file1H1d  = 'cmip5.'+model1[0]+'.historical.r1i1p1.an.ocn.Omon.density.ver-'+model1[1]+'_zon1D.nc'
file2H2d  = 'cmip5.'+model2[0]+'.historical.r1i1p1.an.ocn.Omon.density.ver-'+model2[1]+'_zon2D.nc'
file2H1d  = 'cmip5.'+model2[0]+'.historical.r1i1p1.an.ocn.Omon.density.ver-'+model2[1]+'_zon1D.nc'
file3H2d  = 'cmip5.'+model3[0]+'.historical.r1i1p1.an.ocn.Omon.density.ver-'+model3[1]+'_zon2D.nc'
file3H1d  = 'cmip5.'+model3[0]+'.historical.r1i1p1.an.ocn.Omon.density.ver-'+model3[1]+'_zon1D.nc'

file1HN2d  = 'cmip5.'+model1[0]+'.historicalNat.r1i1p1.an.ocn.Omon.density.ver-'+model1[2]+'_zon2D.nc'
file1HN1d  = 'cmip5.'+model1[0]+'.historicalNat.r1i1p1.an.ocn.Omon.density.ver-'+model1[2]+'_zon1D.nc'
file2HN2d  = 'cmip5.'+model2[0]+'.historicalNat.r1i1p1.an.ocn.Omon.density.ver-'+model2[2]+'_zon2D.nc'
file2HN1d  = 'cmip5.'+model2[0]+'.historicalNat.r1i1p1.an.ocn.Omon.density.ver-'+model2[2]+'_zon1D.nc'
file3HN2d  = 'cmip5.'+model3[0]+'.historicalNat.r1i1p1.an.ocn.Omon.density.ver-'+model3[2]+'_zon2D.nc'
file3HN1d  = 'cmip5.'+model3[0]+'.historicalNat.r1i1p1.an.ocn.Omon.density.ver-'+model3[2]+'_zon1D.nc'

# Define variable  TODO: read as argument
varname = defVar('salinity')
#varname = defVar('temp')
#varname = defVar('depth')
#varname = defVar('volume')
#varname = defVar('persist')

# Define plot name
plotName = 'cmip5_compare_hist_histNat_r1i1p1_'+varname['var']

# years for difference
y11 = 0
y12 = 2

# density domain
domrho = [21.,26.,28.]   # min/mid/max
#
# -------------------------------------------------------------------------------

#-- Define variable properties

var    = varname['var']
minmax = varname['minmax']
clevsm = varname['clevsm']
legVar = varname['legVar']
unit = varname['unit']

#-- Open netcdf files model 1
print
print file1H2d
nchist2d  = open_ncfile(indirHist+'/'+file1H2d)
nchist1d  = open_ncfile(indirHist+'/'+file1H1d)
print file1HN2d
nchistn2d  = open_ncfile(indirHistN+'/'+file1HN2d)
nchistn1d  = open_ncfile(indirHistN+'/'+file1HN1d)
print

agreelev=0.6 # not used

#-- Read variables
# Restrict variables to bowl
tvar  = nchist2d.variables[var]
tvarn = nchistn2d.variables[var]
lev = nchist2d.variables['lev'][:]
lat = nchist2d.variables['latitude'][:]


# Read lightest density of persistent ocean (ptopsigma)
ptopsig  = nchist1d.variables['ptopsigma']
ptopsign = nchistn1d.variables['ptopsigma']

#-- Build plot variables
# difference
vara = np.ma.average(tvar[y11:y12,0,:,:], axis=0)-np.ma.average(tvarn[y11:y12,0,:,:], axis=0)
varp = np.ma.average(tvar[y11:y12,2,:,:], axis=0)-np.ma.average(tvarn[y11:y12,2,:,:], axis=0)
vari = np.ma.average(tvar[y11:y12,3,:,:], axis=0)-np.ma.average(tvarn[y11:y12,3,:,:], axis=0)
# mean
varam = np.ma.average(tvar[:,0,:,:], axis=0)
varpm = np.ma.average(tvar[:,2,:,:], axis=0)
varim = np.ma.average(tvar[:,3,:,:], axis=0)
# Periods ptopsigma
ptopsigha = np.ma.average(ptopsig[y11:y12,0,], axis=0)
ptopsighp = np.ma.average(ptopsig[y11:y12,2,], axis=0)
ptopsighi = np.ma.average(ptopsig[y11:y12,3,], axis=0)
ptopsighna = np.ma.average(ptopsign[y11:y12,0,], axis=0)
ptopsighnp = np.ma.average(ptopsign[y11:y12,2,], axis=0)
ptopsighni = np.ma.average(ptopsign[y11:y12,3,], axis=0)

#-- Create variable bundles
varAtl = {'name': 'Global','diffBowl': vara, 'meanBowl': varam, 'agree': 'none'}
varPac = {'name': 'Pacific','diffBowl': varp, 'meanBowl': varpm, 'agree': 'none'}
varInd = {'name': 'Indian','diffBowl': vari, 'meanBowl': varim, 'agree': 'none'}
vartsiga = {'yr1': ptopsigha, 'yr2': ptopsighna}
vartsigp = {'yr1': ptopsighp, 'yr2': ptopsighnp}
vartsigi = {'yr1': ptopsighi, 'yr2': ptopsighni}
labBowl = ['Hist','HistNat']
 
#-- Create figure and axes instances
fig, axes = plt.subplots(nrows=2,ncols=3,figsize=(17,5))

#-- color map
#fileGMT='palet_diff_6'
#cdict = gmtColormap(fileGMT,GMTPath = '/Users/ericg/POST_IT/config/palettes')
#cmap = LinearSegmentedColormap('campgmt',cdict)
cmap = plt.get_cmap('bwr') # red/white/blue difference map

#
# -------- Make plot ----------------
#
cnplot=zon_2dom(plt,axes[0,0],axes[1,0],lat,lev,varAtl,vartsiga,unit,minmax,clevsm,cmap,domrho,agreelev,False,'F',labBowl)

cnplot=zon_2dom(plt,axes[0,1],axes[1,1],lat,lev,varPac,vartsigp,unit,minmax,clevsm,cmap,domrho,agreelev,False,'T',labBowl)

cnplot=zon_2dom(plt,axes[0,2],axes[1,2],lat,lev,varInd,vartsigi,unit,minmax,clevsm,cmap,domrho,agreelev,False,'R',labBowl)

plt.subplots_adjust(hspace = .00001, wspace=0.05, left=0.04, right=0.86)

#-- Add colorbar
cbar = fig.colorbar(cnplot[0], ax=axes.ravel().tolist(),fraction=0.015, shrink=2.0,pad=0.05)
cbar.set_label(unit)

# add Title text
ttxt = fig.suptitle(legVar+' for diff between Hist and HistNat for model '+model1[0]+' [yrs '+str(y11+1850)+'-'+str(y12+1850)+']', fontsize=14, fontweight='bold')

#-- Output  # TODO read as argument

plt.show()
#plt.savefig(plotName+'.pdf', bbox_inches='tight')

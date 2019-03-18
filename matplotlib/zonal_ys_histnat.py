#!/bin/env python
# -*- coding: utf-8 -*-
"""
Python matplotlib 
Make density/latitude section diff between hist and histNat
Per model and for mme of each
for Atl/Pac/Ind for a number of variables

(c) Eric Guilyardi Feb 2016

TODO: - add arguments for variable and output type
      - read variable unit from file

"""

import os.path
import sys
from inspect import getsourcefile

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MaxNLocator
from netCDF4 import Dataset as open_ncfile

from densit_matplot_lib import zon_2dom, defVar
import colormaps as cmaps

current_path = os.path.abspath(getsourcefile(lambda:0))
current_dir = os.path.dirname(current_path)
parent_dir = current_dir[:current_dir.rfind(os.path.sep)]

sys.path.insert(0, parent_dir)

from libToE import findToE

#import matplotlib as mpl
#from   mpl_toolkits.basemap import Basemap, cm
#from   mpl_toolkits.axes_grid1 import Grid
#from   matplotlib.colors import LinearSegmentedColormap

# -------------------------------------------------------------------------------
#                               Define work
# -------------------------------------------------------------------------------

indir = '/data/ericglod/Density_binning/'
workh = 'Prod_density_april15/mme_hist/'
workhn = 'Prod_density_april15/mme_histNat/'
indirh = indir + workh
indirhn = indir + workhn
file2dh = 'cmip5.multimodel_Nat.historical.ensm.an.ocn.Omon.density_zon2D.nc'
file1dh = 'cmip5.multimodel_Nat.historical.ensm.an.ocn.Omon.density_zon1D.nc'
file2dhn = 'cmip5.multimodel_All.historicalNat.ensm.an.ocn.Omon.density_zon2D.nc'
file1dhn = 'cmip5.multimodel_All.historicalNat.ensm.an.ocn.Omon.density_zon1D.nc'


# Define variable  TODO: read as argument
varname = defVar('salinity')
#varname = defVar('temp')
#varname = defVar('depth')
varname = defVar('volume')
#varname = defVar('persist')
#varname = defVar('heatcontent')

ToE = True
#ToE = False
multStd = 2. # detect ToE at multStd std dev of histNat
# Define plot name
plotName = 'cmip5_hist_vs_histNat_stddev_' + varname['var']
if ToE:
    plotName = 'cmip5_ToE_hist_vs_histNat_stddev_' + varname['var']
valmask = 1.e20
# years for difference (last 10 years)
nyearsComp = 5
y11 = -nyearsComp
y12 = -1
iniyear = -40
labBowl = ['histNat', 'hist']

# density domain
domrho = [21., 26., 28.]  # min/mid/max
delrho = [.5, .2]
#
# -------------------------------------------------------------------------------

# -- Define variable properties

var = varname['var']
minmax = varname['minmax']
clevsm = varname['clevsmdif']
#clevsm = varname['clevsm']
legVar = varname['legVar']
unit = varname['unit']

agreelev = 0. # not used

# -- Open netcdf files
nc1dh = open_ncfile(indirh + '/' + file1dh)
nc2dh = open_ncfile(indirh + '/' + file2dh)
nc1dhn = open_ncfile(indirhn + '/' + file1dhn)
nc2dhn = open_ncfile(indirhn + '/' + file2dhn)

# -- Read variables
# Restrict variables to bowl (hist) TODO add switch for bowl
tvarha = nc2dh.variables[var + 'Bowl'][:, 1, :, :].squeeze()
tvarhp = nc2dh.variables[var + 'Bowl'][:, 2, :, :].squeeze()
tvarhi = nc2dh.variables[var + 'Bowl'][:, 3, :, :].squeeze()
#tvarha = nc2dh.variables[var][:, 1, :, :].squeeze()
#tvarhp = nc2dh.variables[var][:, 2, :, :].squeeze()
#tvarhi = nc2dh.variables[var][:, 3, :, :].squeeze()
lev = nc2dh.variables['lev'][:]
levN = lev.size
lat = nc2dh.variables['latitude'][:]
latN = lat.size
time = nc2dh.variables['time'][:]
timN = time.size
# Restrict variables to bowl (histNat) TODO add switch for bowl
#tvarhna = nc2dhn.variables[var + 'Bowl'][:, 1, :, :].squeeze()
#tvarhnp = nc2dhn.variables[var + 'Bowl'][:, 2, :, :].squeeze()
#tvarhni = nc2dhn.variables[var + 'Bowl'][:, 3, :, :].squeeze()
tvarhna = nc2dhn.variables[var ][:, 1, :, :].squeeze()
tvarhnp = nc2dhn.variables[var ][:, 2, :, :].squeeze()
tvarhni = nc2dhn.variables[var ][:, 3, :, :].squeeze()

# Read lightest density of persistent ocean (ptopsigma)
ptopsigha = nc1dh.variables['ptopsigma'][:, 1, :].squeeze()
ptopsighp = nc1dh.variables['ptopsigma'][:, 2, :].squeeze()
ptopsighi = nc1dh.variables['ptopsigma'][:, 3, :].squeeze()
ptopsighna = nc1dhn.variables['ptopsigma'][:, 1, :].squeeze()
ptopsighnp = nc1dhn.variables['ptopsigma'][:, 2, :].squeeze()
ptopsighni = nc1dhn.variables['ptopsigma'][:, 3, :].squeeze()

# -- Build plot variables
# difference
vara = np.ma.average(tvarha[y11:y12], axis=0) - np.ma.average(tvarhna[y11:y12], axis=0)
varp = np.ma.average(tvarhp[y11:y12], axis=0) - np.ma.average(tvarhnp[y11:y12], axis=0)
vari = np.ma.average(tvarhi[y11:y12], axis=0) - np.ma.average(tvarhni[y11:y12], axis=0)
# Compute significance of difference when diff within 1 stddev of histNat variability (in the MME sense)
varams = np.ma.std(tvarhna, axis=0)
varpms = np.ma.std(tvarhnp, axis=0)
varims = np.ma.std(tvarhni, axis=0)

if ToE:
    # reorganise i,j dims in single dimension data (speeds up loops)
    tvarha  = np.reshape(tvarha, (timN,levN*latN))
    tvarhna = np.reshape(tvarhna,(timN,levN*latN))
    varams  = np.reshape(varams, (levN*latN))
    tvarhp  = np.reshape(tvarhp, (timN,levN*latN))
    tvarhnp = np.reshape(tvarhnp,(timN,levN*latN))
    varpms  = np.reshape(varpms, (levN*latN))
    tvarhi  = np.reshape(tvarhi, (timN,levN*latN))
    tvarhni = np.reshape(tvarhni,(timN,levN*latN))
    varims  = np.reshape(varims, (levN*latN))

    # Compute ToE as last date when diff hist-histNat is larger than mult * stddev
    varam = np.reshape(findToE(tvarha-tvarhna, varams, multStd)+iniyear,(levN,latN))
    varpm = np.reshape(findToE(tvarhp-tvarhnp, varpms, multStd)+iniyear,(levN,latN))
    varim = np.reshape(findToE(tvarhi-tvarhni, varims, multStd)+iniyear,(levN,latN))

    # shade ToE and contour diff hist-histNat
    tmpa = vara
    vara = varam+1900
    varam = tmpa
    tmpp = varp
    varp = varpm+1900
    varpm = tmpp
    tmpi = vari
    vari = varim+1900
    varim = tmpi
else:
    # std dev TODO: add switch
    varam = varams
    varpm = varpms
    varim = varims
    # mean field
    varam = np.ma.average(tvarha[y11:y12], axis=0)
    varpm = np.ma.average(tvarhp[y11:y12], axis=0)
    varim = np.ma.average(tvarhi[y11:y12], axis=0)


# Not used
varaa=0.
varap=0.
varai=0.
# Periods ptopsigma
ptopsigyr1a = np.ma.average(ptopsighna[y11:y12], axis=0)
ptopsigyr1p = np.ma.average(ptopsighnp[y11:y12], axis=0)
ptopsigyr1i = np.ma.average(ptopsighni[y11:y12], axis=0)
ptopsigyr2a = np.ma.average(ptopsigha[y11:y12], axis=0)
ptopsigyr2p = np.ma.average(ptopsighp[y11:y12], axis=0)
ptopsigyr2i = np.ma.average(ptopsighi[y11:y12], axis=0)

# -- Create variable bundles
varAtl = {'name': 'Atlantic', 'diffBowl': vara, 'meanBowl': varam, 'agree': varaa}
varPac = {'name': 'Pacific', 'diffBowl': varp, 'meanBowl': varpm, 'agree': varap}
varInd = {'name': 'Indian', 'diffBowl': vari, 'meanBowl': varim, 'agree': varai}
vartsiga = {'yr1': ptopsigyr1a, 'yr2': ptopsigyr2a} # first is dashed , second solid
vartsigp = {'yr1': ptopsigyr1p, 'yr2': ptopsigyr2p}
vartsigi = {'yr1': ptopsigyr1i, 'yr2': ptopsigyr2i}

# -- Create figure and axes instances
fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(17, 5))

# -- color map
cmap = plt.get_cmap('bwr')  # red/white/blue difference map
if ToE:
    plt.register_cmap(name='viridis', cmap=cmaps.viridis)
    cmap=cmaps.viridis
    #cmap = plt.get_cmap('viridis')
    minmax = [1900, 2000, 20]
    unit = 'ToE'

#
# -------- Make plot ----------------
#
cnplot = zon_2dom(plt, axes[0, 0], axes[1, 0], lat, lev, varAtl, vartsiga, unit, minmax, clevsm, cmap, domrho, agreelev,
                  False, 'F', labBowl)

cnplot = zon_2dom(plt, axes[0, 1], axes[1, 1], lat, lev, varPac, vartsigp, unit, minmax, clevsm, cmap, domrho, agreelev,
                  False, 'T', labBowl)

cnplot = zon_2dom(plt, axes[0, 2], axes[1, 2], lat, lev, varInd, vartsigi, unit, minmax, clevsm, cmap, domrho, agreelev,
                  False, 'R', labBowl)

plt.subplots_adjust(hspace=.00001, wspace=0.05, left=0.04, right=0.86)

# -- Add colorbar
levels = MaxNLocator(nbins=minmax[2]).tick_values(minmax[0], minmax[1])
cbar = fig.colorbar(cnplot[0], ax=axes.ravel().tolist(), fraction=0.015, shrink=2.0, pad=0.05)
#extend='both',boundaries=[1880] + levels + [2010], extendfrac='auto')
cbar.set_label(unit)

# add Title text
titleText=legVar + ' [' + workh +'] - Hist minus HistNat '
if ToE:
    titleText = titleText+'(contour) - & MME ToE [> '+str(multStd)+' std] (shading)'
else:
    titleText = titleText+' (last '+str(nyearsComp)+' yrs)'

ttxt = fig.suptitle(titleText, fontsize=14, fontweight='bold')
# -- Output  # TODO read as argument

plt.show()
#plt.savefig(plotName+'.pdf', bbox_inches='tight')

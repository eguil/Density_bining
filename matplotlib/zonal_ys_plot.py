#!/bin/env python
# -*- coding: utf-8 -*-
"""
Python matplotlib 
Make density/latitude section for Atl/Pac/Ind for a number of variables

(c) Eric Guilyardi Feb 2016

TODO: - add arguments for variable and output type
      - read variable unit from file

"""
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset as open_ncfile
from densit_matplot_lib import zon_2dom, defVar


# -------------------------------------------------------------------------------
#                               Define work
# -------------------------------------------------------------------------------

indir = '/Users/ericg/Projets/Density_bining/'

# description of work (dow)
dow = 'model'
dow = 'EN4'
dow = 'ishii'

# output format
outfmt = 'view'
#outfmt = 'save'

# models
if dow == 'model':
    work = 'Prod_density_april15/mme_hist'
    file2d = 'cmip5.multimodel_All.historical.ensm.an.ocn.Omon.density_zon2D.nc'
    file1d = 'cmip5.multimodel_All.historical.ensm.an.ocn.Omon.density_zon1D.nc'

# observations
if dow == 'EN4':
    work = 'Prod_density_obs_april16/mme_obs'
    file2d = 'obs.EN4.historical.ensm.an.ocn.Omon.density.ver-1.latestX_zon2D.nc'
    file1d = 'obs.EN4.historical.ensm.an.ocn.Omon.density.ver-1.latestX_zon1D.nc'
if dow == 'ishii':
    work = 'Prod_density_obs_april16/mme_obs'
    file2d = 'obs.Ishii.historical.ensm.an.ocn.Omon.density.ver-1.latestX_zon2D.nc'
    file1d = 'obs.Ishii.historical.ensm.an.ocn.Omon.density.ver-1.latestX_zon1D.nc'

indir = indir + work

# Model agreement level
agreelev = 0.6

# Define variable  TODO: read as argument
varname = defVar('salinity')
#varname = defVar('temp')
varname = defVar('depth')
#varname = defVar('volume')
# varname = defVar('persist')

# Define plot name, years for difference, bowl label and bowl and model agreement plot options

if dow == 'model': # MME for 1861-2005 (146 time steps)
    plotName = 'cmip5_mme_hist_42models_' + varname['var']
    y11 = 140 ; y12 = 145
    y21 = 0 ; y22 = 50
    labBowl = ['<1950', '2000']
    restrictBowl = True
    modelAgree = True

if dow == 'EN4': # 1900.01 - 2015.04 (115 time steps, ignore last year) Good et al.
    plotName = 'obs_EN4_' + varname['var']
    y11 = 108 ; y12 = 113
    y21 = 0 ; y22 = 30
    labBowl = ['<1930', '2010']
    restrictBowl = True
    modelAgree = False

if dow == 'ishii': # 1945.01 - 2012.12 (68 time steps)
    plotName = 'obs_ishii_' + varname['var']
    y11 = 62 ; y12 = 67
    y21 = 0 ; y22 = 5
    labBowl = ['<1950', '2010']
    restrictBowl = True
    modelAgree = False
#(ignore for now JAMSTEC: 2001.01 - 2014.12
#(ignore for now UCSD: 2004.01 - 2015.03
# SmithAndMurphy2007 # 1950.01 - 2013.02 (ignore last year)

# density domain
domrho = [21., 26., 28.]  # min/mid/max
delrho = [.5, .2]
#
# -------------------------------------------------------------------------------

# -- Define variable properties

var = varname['var']
minmax = varname['minmax']
clevsm = varname['clevsm']
legVar = varname['legVar']
unit = varname['unit']

# -- Open netcdf files
nc1d = open_ncfile(indir + '/' + file1d)
nc2d = open_ncfile(indir + '/' + file2d)

# -- Read variables
# Restrict variables to bowl
if restrictBowl:
    varb = 'Bowl'
else:
    varb=''
tvara = nc2d.variables[var + varb][:, 1, :, :].squeeze()
tvarp = nc2d.variables[var + varb][:, 2, :, :].squeeze()
tvari = nc2d.variables[var + varb][:, 3, :, :].squeeze()
lev = nc2d.variables['lev'][:]
lat = nc2d.variables['latitude'][:]
# Read model agreement variables
if modelAgree:
    tvaraa = nc2d.variables[var + 'Agree'][:, 1, :, :].squeeze()
    tvarap = nc2d.variables[var + 'Agree'][:, 2, :, :].squeeze()
    tvarai = nc2d.variables[var + 'Agree'][:, 3, :, :].squeeze()
else:
    tvaraa = 1.; tvarap=1. ; tvarai = 1.
# Read lightest density of persistent ocean (ptopsigma)
ptopsiga = nc1d.variables['ptopsigma'][:, 1, :].squeeze()
ptopsigp = nc1d.variables['ptopsigma'][:, 2, :].squeeze()
ptopsigi = nc1d.variables['ptopsigma'][:, 3, :].squeeze()

# -- Build plot variables
# difference
vara = np.ma.average(tvara[y11:y12], axis=0) - np.ma.average(tvara[y21:y22], axis=0)
varp = np.ma.average(tvarp[y11:y12], axis=0) - np.ma.average(tvarp[y21:y22], axis=0)
vari = np.ma.average(tvari[y11:y12], axis=0) - np.ma.average(tvari[y21:y22], axis=0)
# mean
varam = np.ma.average(tvara, axis=0)
varpm = np.ma.average(tvarp, axis=0)
varim = np.ma.average(tvari, axis=0)
# Average model agreement over final period
if modelAgree:
    varaa = np.ma.average(tvaraa[y11:y12], axis=0)
    varap = np.ma.average(tvarap[y11:y12], axis=0)
    varai = np.ma.average(tvarai[y11:y12], axis=0)
else:
    varaa = 1.; varap=1. ; varai = 1.
# Periods ptopsigma
ptopsigyr1a = np.ma.average(ptopsiga[y11:y12], axis=0)
ptopsigyr1p = np.ma.average(ptopsigp[y11:y12], axis=0)
ptopsigyr1i = np.ma.average(ptopsigi[y11:y12], axis=0)
ptopsigyr2a = np.ma.average(ptopsiga[y21:y22], axis=0)
ptopsigyr2p = np.ma.average(ptopsigp[y21:y22], axis=0)
ptopsigyr2i = np.ma.average(ptopsigi[y21:y22], axis=0)

# -- Create variable bundles
varAtl = {'name': 'Atlantic', 'diffBowl': vara, 'meanBowl': varam, 'agree': varaa}
varPac = {'name': 'Pacific', 'diffBowl': varp, 'meanBowl': varpm, 'agree': varap}
varInd = {'name': 'Indian', 'diffBowl': vari, 'meanBowl': varim, 'agree': varai}
vartsiga = {'yr1': ptopsigyr1a, 'yr2': ptopsigyr2a}
vartsigp = {'yr1': ptopsigyr1p, 'yr2': ptopsigyr2p}
vartsigi = {'yr1': ptopsigyr1i, 'yr2': ptopsigyr2i}

# -- Create figure and axes instances
fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(17, 5))

# -- color map
# fileGMT='palet_diff_6'
# cdict = gmtColormap(fileGMT,GMTPath = '/Users/ericg/POST_IT/config/palettes')
# cmap = LinearSegmentedColormap('campgmt',cdict)
cmap = plt.get_cmap('bwr')  # red/white/blue difference map

#
# -------- Make plot ----------------
#
cnplot = zon_2dom(plt, axes[0, 0], axes[1, 0], lat, lev, varAtl, vartsiga, unit, minmax, clevsm, cmap, domrho, agreelev,
                  modelAgree, 'F', labBowl)

cnplot = zon_2dom(plt, axes[0, 1], axes[1, 1], lat, lev, varPac, vartsigp, unit, minmax, clevsm, cmap, domrho, agreelev,
                  modelAgree, 'T', labBowl)

cnplot = zon_2dom(plt, axes[0, 2], axes[1, 2], lat, lev, varInd, vartsigi, unit, minmax, clevsm, cmap, domrho, agreelev,
                  modelAgree, 'R', labBowl)

plt.subplots_adjust(hspace=.00001, wspace=0.05, left=0.04, right=0.86)

# -- Add colorbar
cbar = fig.colorbar(cnplot[0], ax=axes.ravel().tolist(), fraction=0.015, shrink=2.0, pad=0.05)
cbar.set_label(unit)

# add Title text
ttxt = fig.suptitle(legVar + ' for ' + work+'/ '+dow, fontsize=14, fontweight='bold')

# -- Output  # TODO read as argument

if outfmt == 'view':
    plt.show()
else:
    print 'Save',plotName+'.pdf'
    plt.savefig(plotName+'.pdf', bbox_inches='tight')

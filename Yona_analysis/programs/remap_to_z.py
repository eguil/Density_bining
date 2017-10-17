#!/bin/env python
# -*- coding: utf-8 -*-

import os
from netCDF4 import Dataset as open_ncfile
import matplotlib.pyplot as plt
from maps_matplot_lib import defVarmme, remapToZ, zon_2Dz, custom_div_cmap
import numpy as np


# -------------------------------------------------------------------------------
#                                Define work
# -------------------------------------------------------------------------------

# -- Choose what to compute
name = 'mme_hist_histNat'
# name = 'mme_1pctCO2vsPiC'

# -- Choose where to stop for 1%CO2 simulations : 2*CO2 (70 years) or 4*CO2 (140 years) or 1.4*CO2 (34 years0
focus_1pctCO2 = '2*CO2'  # 1.4 or 2*CO2 or 4*CO2

# output format
outfmt = 'view'
# outfmt = 'save'

valmask = 1.e20
basinN = 4

# -- Choose work files

if name == 'mme_hist_histNat':
    indirh = '/data/ericglod/Density_binning/Prod_density_april15/historical/'
    # fileh_2d = 'cmip5.multimodel_Nat.historical.ensm.an.ocn.Omon.density_zon2D.nc'
    fileh_2d = 'cmip5.CCSM4.historical.r4i1p1.an.ocn.Omon.density.ver-v20121128_zon2D.nc'
    # fileh_1d = 'cmip5.multimodel_Nat.historical.ensm.an.ocn.Omon.density_zon1D.nc'
    fileh_1d = 'cmip5.CCSM4.historical.r4i1p1.an.ocn.Omon.density.ver-v20121128_zon1D.nc'
    datah_2d = indirh + fileh_2d; datah_1d = indirh + fileh_1d
    indirhn = '/data/ericglod/Density_binning/Prod_density_april15/historicalNat/'
    # filehn_2d = 'cmip5.multimodel_All.historicalNat.ensm.an.ocn.Omon.density_zon2D.nc'
    filehn_2d = 'cmip5.CCSM4.historicalNat.r4i1p1.an.ocn.Omon.density.ver-v20120614_zon2D.nc'
    # filehn_1d = 'cmip5.multimodel_All.historicalNat.ensm.an.ocn.Omon.density_zon1D.nc'
    filehn_1d = 'cmip5.CCSM4.historicalNat.r4i1p1.an.ocn.Omon.density.ver-v20120614_zon1D.nc'
    datahn_2d = indirhn + filehn_2d; datahn_1d = indirhn + filehn_1d
    fh2d = open_ncfile(datah_2d,'r')
    fh1d = open_ncfile(datah_1d,'r')
    fhn2d = open_ncfile(datahn_2d,'r')
    fhn1d = open_ncfile(datahn_1d,'r')

if name == 'mme_1pctCO2vsPiC':
    indir_1pctCO2 = '/data/ericglod/Density_binning/Prod_density_april15/mme_1pctCO2/'
    file_2d = 'cmip5.multimodel_piCtl.1pctCO2.ensm.an.ocn.Omon.density_zon2D.nc'
    # file_2d = 'cmip5.GFDL-ESM2G.1pctCO2.ensm.an.ocn.Omon.density.ver-v20120820_zon2D.nc'
    file_1d = 'cmip5.multimodel_piCtl.1pctCO2.ensm.an.ocn.Omon.density_zon1D.nc'
    # file_1d = 'cmip5.GFDL-ESM2G.1pctCO2.ensm.an.ocn.Omon.density.ver-v20120820_zon1D.nc'
    data_2d = indir_1pctCO2 + file_2d
    data_1d = indir_1pctCO2 + file_1d
    fh2d = open_ncfile(data_2d,'r')
    fh1d = open_ncfile(data_1d,'r')

    indir_piC = '/data/ericglod/Density_binning/Prod_density_april15/mme_piControl/'
    file_2d = 'cmip5.multimodel_1pct.piControl.ensm.an.ocn.Omon.density_zon2D.nc'
    # file_2d = 'cmip5.GFDL-ESM2G.piControl.ensm.an.ocn.Omon.density.ver-v20110601_zon2D.nc'
    file_1d = 'cmip5.multimodel_1pct.piControl.ensm.an.ocn.Omon.density_zon1D.nc'
    # file_1d = 'cmip5.GFDL-ESM2G.piControl.ensm.an.ocn.Omon.density.ver-v20110601_zon1D.nc'
    data_2d = indir_piC + file_2d
    data_1d = indir_piC + file_1d
    fhn2d = open_ncfile(data_2d,'r')
    fhn1d = open_ncfile(data_1d,'r')


# -------------------------------------------------------------------------------
#                                Build variables
# -------------------------------------------------------------------------------

# == Define variables ==
varname = defVarmme('salinity'); v = 'S'
# varname = defVarmme('temp'); v = 'T'
# varname = defVarmme('volume'); v = 'V'

minmax = varname['minmax_zonal']
clevsm = varname['clevsm_zonal']
legVar = varname['legVar']
unit = varname['unit']

var = varname['var_zonal_w/bowl']

# == Read variables ==
lat = fh2d.variables['latitude'][:]

if name == 'mme_hist_histNat':
    field2r = fh2d.variables[var][-10:,:,:,:] # historical
    field1r = fhn2d.variables[var][-10:,:,:] # historicalNat
    depthr = fh2d.variables['isondepth'][-10:,:,:,:]
    volumr = fh2d.variables['isonvol'][-10:,:,:,:]
    bowl2z = fh1d.variables['ptopdepth'][-5:,:,:]
    bowl1z = fhn1d.variables['ptopdepth'][-5:,:,:]
    labBowl = ['histNat', 'hist']

if name == 'mme_1pctCO2vsPiC':
    if focus_1pctCO2 == '4*CO2':
        y1 = 134
        y2 = 140
    if focus_1pctCO2 == '2*CO2':
        y1 = 69
        y2 = 74
    if focus_1pctCO2 == '1.4*CO2':
        y1 = 33
        y2 = 38
    field2r = fh2d.variables[var][y1:y2,:,:,:] # 1pctCO2
    field1r = fhn2d.variables[var][-20:,:,:,:] # PiControl
    depthr = fh2d.variables['isondepth'][y1:y2,:,:,:]
    volumr = fh2d.variables['isonvol'][y1:y2,:,:,:]
    bowl2z = fh1d.variables['ptopdepth'][y1:y2,:,:]
    bowl1z = fhn1d.variables['ptopdepth'][-10:,:,:]
    labBowl = ['PiControl', focus_1pctCO2]

# if v != 'V':
#     field2r[np.ma.nonzero(field2r>60)] = np.ma.masked
#     field1r[np.ma.nonzero(field1r>60)] = np.ma.masked

# == Compute signal hist - histNat or 1pctCO2 - PiControl ==
vardiffr = np.ma.average(field2r, axis=0) - np.ma.average(field1r, axis=0)
# == Average other variables ==
depthr = np.ma.average(depthr, axis=0)
volumr = np.ma.average(volumr, axis=0)
bowl2z = np.ma.average(bowl2z, axis=0)
bowl1z = np.ma.average(bowl1z, axis=0)

# # Mask data
# depthr[np.ma.nonzero(depthr>valmask/2)] = np.ma.masked
# volumr[np.ma.nonzero(volumr>valmask/2)] = np.ma.masked
# vardiffr[np.ma.nonzero(np.abs(vardiffr)>valmask/2)] = np.ma.masked

# == Bathymetry ==
# Read masks
fmask = open_ncfile('/home/ysilvy/Density_bining/Yona_analysis/data/170224_WOD13_masks.nc','r')
basinmask = fmask.variables['basinmask3'][:] # (latitude, longitude)
depthmask = fmask.variables['depthmask'][:] # (latitude, longitude)
# Create basin masks
mask_a = basinmask != 1
mask_p = basinmask != 2
mask_i = basinmask != 3
# Take max value along lon for each latitude of each basin --> gives us the bathymetry
depthmask_a = np.ma.array(depthmask, mask=mask_a)
bathy_a = np.ma.max(depthmask_a, axis=1)
depthmask_p = np.ma.array(depthmask, mask=mask_p)
bathy_p = np.ma.max(depthmask_p, axis=1)
depthmask_i = np.ma.array(depthmask, mask=mask_i)
bathy_i = np.ma.max(depthmask_i, axis=1)

bathy = np.ma.masked_all((basinN,len(lat)))
bathy[1,:] = bathy_a
bathy[2,:] = bathy_p
bathy[3,:] = bathy_i


# == Target grid for remapping ==

# targetz = [0.,5.,15.,25.,35.,45.,55.,65.,75.,85.,95.,105.,116.,128.,142.,158.,181.,216.,272.,364.,511.,732.,
# 1033.,1405.,1830.,2289.,2768.,3257.,3752.,4350.,4749.,5250.]

# WOA09 grid
# targetz = [0, 10, 20, 30, 50, 75, 100, 125, 150, 200, 250, 300, 400, 500, 600,
# 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1750, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 5500]

# WOA13 grid
targetz = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80,
   85, 90, 95, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 375,
   400, 425, 450, 475, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950,
   1000, 1050, 1100, 1150, 1200, 1250, 1300, 1350, 1400, 1450, 1500, 1550,
   1600, 1650, 1700, 1750, 1800, 1850, 1900, 1950, 2000, 2100, 2200, 2300,
   2400, 2500, 2600, 2700, 2800, 2900, 3000, 3100, 3200, 3300, 3400, 3500,
   3600, 3700, 3800, 3900, 4000, 4100, 4200, 4300, 4400, 4500, 4600, 4700,
   4800, 4900, 5000, 5100, 5200, 5300, 5400, 5500]


# == Remap ==
fieldz = remapToZ(vardiffr.data, depthr.data, volumr.data, targetz, bowl1z, v, bathy)


# -- Make variable bundles for each basin
varAtl = {'name': 'Atlantic', 'var_change': fieldz[1,:,:], 'bowl1': bowl1z[1,:], 'bowl2': bowl2z[1,:], 'labBowl': labBowl}
varPac = {'name': 'Pacific', 'var_change': fieldz[2,:,:], 'bowl1': bowl1z[2,:], 'bowl2': bowl2z[2,:], 'labBowl': labBowl}
varInd = {'name': 'Indian', 'var_change': fieldz[3,:,:], 'bowl1': bowl1z[3,:], 'bowl2': bowl2z[3,:], 'labBowl': labBowl}


# -------------------------------------------------------------------------------
#                                Plot
# -------------------------------------------------------------------------------

domzed = [0,500,2000]

# -- Create figure and axes instances
fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(17, 5))

# -- color map
cmap = custom_div_cmap()

# -- levels
levels = np.linspace(minmax[0],minmax[1],minmax[2]) # Change
#levels = np.arange(33.5,35.5,0.2) # Mean salinity
#levels = np.arange(-2,30,2) # Mean temperature
#levels = np.arange(0,701,50) # Mean volume

# -- Contourf of signal
cnplot = zon_2Dz(plt, axes[0,0], axes[1,0], 'left', lat, targetz, varAtl,
                 clevsm, cmap, levels, domzed)
cnplot = zon_2Dz(plt, axes[0,1], axes[1,1], 'mid', lat, targetz, varPac,
                 clevsm, cmap, levels, domzed)
cnplot = zon_2Dz(plt, axes[0,2], axes[1,2], 'right', lat, targetz, varInd,
                 clevsm, cmap, levels, domzed)


plt.subplots_adjust(hspace=.0001, wspace=0.05, left=0.04, right=0.86)

# -- Add colorbar
cb = plt.colorbar(cnplot, ax=axes.ravel().tolist(), ticks=levels[::3], fraction=0.015, shrink=2.0, pad=0.05)
cb.set_label('%s (%s)' % (legVar, unit), fontweight='bold')

# -- Add Title text
if name == 'mme_hist_histNat':
    plotTitle = 'Remapping %s changes %s from rho to z' %(legVar,name)
if name == 'mme_1pctCO2vsPiC':
    plotTitle  = 'Remapping %s changes %s (%s) from rho to z' %(legVar,name,focus_1pctCO2)
fig.suptitle(plotTitle, fontsize=14, fontweight='bold')

plt.figtext(.5,.01,'Computed by : remap_to_z.py',fontsize=9,ha='center')

# -- Output  # TODO read as argument
if outfmt == 'view':
    plt.show()



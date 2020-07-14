#!/bin/env python
# -*- coding: utf-8 -*-

import os
from netCDF4 import Dataset as open_ncfile
import matplotlib.pyplot as plt
from maps_matplot_lib import defVarmme, defVarDurack, remapToZ
import numpy as np

# -------------------------------------------------------------------------------
#                                Define work
# -------------------------------------------------------------------------------

# -- Choose what to compute
# name = 'mme_hist'
name = 'mme_1pctCO2'
# name = 'mme_rcp85'

# -- Choose where to stop for 1%CO2 simulations : 2*CO2 (70 years) or 4*CO2 (140 years) or 1.4*CO2 (34 years)
focus_1pctCO2 = '2*CO2'  # 1.4 or 2*CO2 or 4*CO2

# -- Choose work files

if name == 'mme_hist':
    indirh = '/data/ericglod/Density_binning/Prod_density_april15/mme_hist/'
    fileh_2d = 'cmip5.multimodel_Nat.historical.ensm.an.ocn.Omon.density_zon2D.nc'
    fileh_1d = 'cmip5.multimodel_Nat.historical.ensm.an.ocn.Omon.density_zon1D.nc'
    datah_2d = indirh + fileh_2d; datah_1d = indirh + fileh_1d
    fh2d = open_ncfile(datah_2d,'r')
    fh1d = open_ncfile(datah_1d,'r')

if name == 'mme_1pctCO2':
    indir_1pctCO2 = '/data/ericglod/Density_binning/Prod_density_april15/mme_1pctCO2/'
    file_2d = 'cmip5.multimodel_piCtl.1pctCO2.ensm.an.ocn.Omon.density_zon2D.nc'
    file_1d = 'cmip5.multimodel_piCtl.1pctCO2.ensm.an.ocn.Omon.density_zon1D.nc'
    data_2d = indir_1pctCO2 + file_2d
    data_1d = indir_1pctCO2 + file_1d
    fh2d = open_ncfile(data_2d,'r')
    fh1d = open_ncfile(data_1d,'r')

if name == 'mme_rcp85':
    indir_rcp85 = '/data/ericglod/Density_binning/Prod_density_april15/mme_rcp85/'
    filercp85_2d = 'cmip5.multimodel_Nat.rcp85.ensm.an.ocn.Omon.density_zon2D.nc'
    filercp85_1d = 'cmip5.multimodel_Nat.rcp85.ensm.an.ocn.Omon.density_zon1D.nc'
    fh2d = open_ncfile(indir_rcp85 + filercp85_2d, 'r')
    fh1d = open_ncfile(indir_rcp85 + filercp85_1d, 'r')

# -------------------------------------------------------------------------------
#                                Read variables
# -------------------------------------------------------------------------------

# Read dimension variables
lat = fh2d.variables['latitude'][:]
lev = fh2d.variables['lev'][:]
basin = fh2d.variables['basin'][:]

# Repeat density into (basin,density,latitude)
lat2d,density2d = np.meshgrid(lat,lev)
density3d = np.repeat(density2d[np.newaxis,:,:],4,axis=0)

# Read isondepth, isonvol, ptopdepth, last 5 years
if name != 'mme_1pctCO2':
    depthr = fh2d.variables['isondepth'][-5:,:,:,:]
    volumr = fh2d.variables['isonvol'][-5:,:,:,:]
    bowlz = fh1d.variables['ptopdepth'][-5:,:,:]

else:
    if focus_1pctCO2 == '4*CO2':
        y1 = 134
        y2 = 140
    if focus_1pctCO2 == '2*CO2':
        y1 = 69
        y2 = 74
    if focus_1pctCO2 == '1.4*CO2':
        y1 = 33
        y2 = 38
    depthr = fh2d.variables['isondepth'][y1:y2,:,:,:]
    volumr = fh2d.variables['isonvol'][y1:y2,:,:,:]
    bowlz = fh1d.variables['ptopdepth'][y1:y2,:,:]

# Average
depthr = np.ma.average(depthr, axis=0)
volumr = np.ma.average(volumr, axis=0)
bowlz = np.ma.average(bowlz, axis=0)

# -------------------------------------------------------------------------------
#                                Remapping
# -------------------------------------------------------------------------------

# Target grid for remapping

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
density_z = remapToZ(density3d, depthr, volumr, targetz, bowlz, bathy=None)
print('density remapped')

# -------------------------------------------------------------------------------
#                                Save in file
# -------------------------------------------------------------------------------

if name == 'mme_hist':
    fileName = 'cmip5.multimodel_Nat.historical.remaptoz_density.zon2D.nc'
if name == 'mme_rcp85':
    fileName = 'cmip5.multimodel_Nat.rcp85.remaptoz_density.zon2D.nc'
if name == 'mme_1pctCO2':
    fileName = 'cmip5.multimodel_piCtl.1pctCO2.'+focus_1pctCO2+'.remaptoz_density.zon2D.nc'

dir = '/home/ysilvy/Density_bining/Yona_analysis/data/remaptoz/density/'

fout = open_ncfile(dir+fileName,'w', format='NETCDF4')

fout.description = 'Density field remapped back to pseudo-depth coordinate, using isondepth and weighting' \
                   'by isonvol (averaged over last 5 years).'

# dimensions
fout.createDimension('basin',4)
fout.createDimension('pseudo_depth',len(targetz))
fout.createDimension('latitude',len(lat))

# variables
basin = fout.createVariable('basin', 'i4', ('basin',))
pseudo_depth = fout.createVariable('pseudo_depth', 'i4', ('pseudo_depth',))
latitude = fout.createVariable('latitude', 'f4', ('latitude',))
density = fout.createVariable('density', 'f4', ('basin','pseudo_depth','latitude'))

# data
basin[:] = fh2d.variables['basin'][:]
latitude[:] = fh2d.variables['latitude'][:]
pseudo_depth[:] = targetz
density[:,:,:] = density_z

#units
basin.units = 'basin index : 0=Global, 1=Atlantic, 2=Pacific, 3=Indian'
density.units = 'kg.m-3'

fout.close()
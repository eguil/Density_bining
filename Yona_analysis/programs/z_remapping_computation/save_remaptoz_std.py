#!/bin/env python
# -*- coding: utf-8 -*-

"""
Python matplotlib
Compute pseudo-depth/latitude standard deviation of salinity for each basin, each model, and save in output file
Choose between :
- HistoricalNat (ensemble mean, max std)
- PiControl

"""

import numpy as np
from netCDF4 import Dataset as open_ncfile
from maps_matplot_lib import defVarmme, remapToZ
import glob, os

# ===== Define Work ======

# work  = 'histNat'
work = 'piC'

indir_histNat = '/data/ericglod/Density_binning/Prod_density_april15/mme_histNat/'
indir_piC = '/data/ericglod/Density_binning/Prod_density_april15/mme_piControl/'
fig_dir = '/home/ysilvy/Density_bining/Yona_analysis/figures/models/zonal_remaptoz/'

varname = defVarmme('salinity'); v = 'S'
# varname = defVarmme('temp'); v = 'T'

var = varname['var_zonal_w/bowl']
legVar = varname['legVar']
unit = varname['unit']

# Density domain
rhomin = 21
rhomid = 26
rhomax = 28
domrho = [rhomin, rhomid, rhomax]

if work == 'histNat':
    indir = indir_histNat
    longName = 'historicalNat'
    dir = 'histNat_std/'
else:
    indir = indir_piC
    longName = 'piControl'
    dir = 'piControl_std/'

# ==== Read std for each model in list ====

# -- Read files in directory
listfiles = sorted(glob.glob(indir + '/*2D.nc'))
nmodels = len(listfiles)

# -- Loop on models
# for i in range(nmodels):
i=0
file = os.path.basename(listfiles[i])
f = open_ncfile(indir+file,'r')
name = file.split('.')[1] # Read model name

print('Reading '+work+' for '+name)

if work == 'histNat':
    if name != 'multimodel_Nat':
        # Read varstd of histNat (max std of all runs for each model)
        stdvar = f.variables[var+'Std'][:,:,:]
    else:
        # Read var of histNat for the multimodel
        varmme = f.variables[var][:]
        stdvar = np.ma.std(varmme[:,:,:,:],axis=0)

else:
    # Read PiControl over 240 years + compute std of PiControl
    varpiC = f.variables[var][-240:,:,:,:]
    stdvar = np.ma.std(varpiC[:,:,:,:], axis=0)

density = f.variables['lev'][:]
lat = f.variables['latitude'][:]
depthr = np.ma.average(f.variables['isondepth'][:,:,:,:],axis=0)
volumr = np.ma.average(f.variables['isonvol'][:,:,:,:],axis=0)

# Read bowl and take the average
file2 = glob.glob(indir+'/*'+name+'*1D.nc')[0]
f2 = open_ncfile(file2,'r')
if work == 'histNat':
    bowl = f2.variables['ptopsigma'][:]
else:
    bowl = f2.variables['ptopsigma'][-240:,:,:]
bowl = np.ma.average(bowl,axis=0)


# ==== Remap std ====

# -- Target grid for remapping
# WOA13 grid
targetz = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80,
   85, 90, 95, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 375,
   400, 425, 450, 475, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950,
   1000, 1050, 1100, 1150, 1200, 1250, 1300, 1350, 1400, 1450, 1500, 1550,
   1600, 1650, 1700, 1750, 1800, 1850, 1900, 1950, 2000, 2100, 2200, 2300,
   2400, 2500, 2600, 2700, 2800, 2900, 3000, 3100, 3200, 3300, 3400, 3500,
   3600, 3700, 3800, 3900, 4000, 4100, 4200, 4300, 4400, 4500, 4600, 4700,
   4800, 4900, 5000, 5100, 5200, 5300, 5400, 5500]

# -- Remap
stdvar_z = remapToZ(stdvar, depthr, volumr, targetz, bowlz=bowl, bathy=None)

# ==== Save std in output file ====

version = file.split('.')[-2]
if name ==  'multimodel_Nat' or name == 'multimodel_1pctCO2':
    version = version.split('_')[-1]
fileOut = 'cmip5.'+name+'.'+longName+'.std.ensm.an.ocn.Omon.remaptoz.'+version+'.nc'

if work == 'histNat':
    dir = '/home/ysilvy/Density_bining/Yona_analysis/data/remaptoz/histNat_std/'
else:
    dir = '/home/ysilvy/Density_bining/Yona_analysis/data/remaptoz/piControl_std/'

fout = open_ncfile(dir+fileOut,'w', format='NETCDF4')

fout.description = 'Interannual standard deviation in lat/rho remapped back to pseudo-depth coordinate.'

# dimensions
fout.CreateDimension('basin',4)
fout.createDimension('pseudo_depth',len(targetz))
fout.createDimension('latitude',len(lat))

# variables
basin = fout.createVariable('basin', 'i4', ('basin',))
pseudo_depth = fout.createVariable('pseudo_depth', 'i4', ('pseudo_depth',))
latitude = fout.createVariable('latitude', 'f4', ('latitude',))
varStd = fout.createVariable(var+'Std', 'f4', ('basin','pseudo_depth','latitude'))

# data
basin[:] = f.variables['basin'][:]
latitude[:] = lat
pseudo_depth[:] = targetz
varStd[:,:,:] = stdvar_z

#units
basin.units = 'basin index : 0=Global, 1=Atlantic, 2=Pacific, 3=Indian'
varStd.units = unit
if work == 'histNat':
    varStd.long_name = 'Interannual standard deviation of '+ legVar
else:
    varStd.long_name = 'Interannual standard deviation of '+ legVar + ' for the last 240 years of the piControl'

fout.close()
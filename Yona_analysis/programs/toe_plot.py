#!/bin/env python
# -*- coding: utf-8 -*-

"""
Python matplotlib
Plot signal (hist-histNat) and noise (2*std of histNat)


"""

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset as open_ncfile
from maps_matplot_lib import zonal_2D, defVarmme, averageDom
from modelsDef import defModels, defModelsCO2piC
from libToE import ToEdomainhistvshistNat

# ----- Workspace ------

indirh = '/data/ericglod/Density_binning/Prod_density_april15/mme_hist/'
fileh = 'cmip5.multimodel_Nat.historical.ensm.an.ocn.Omon.density_zon2D.nc'
indirhn = '/data/ericglod/Density_binning/Prod_density_april15/mme_histNat/'
filehn = 'cmip5.multimodel_All.historicalNat.ensm.an.ocn.Omon.density_zon2D.nc'
fh = open_ncfile(indirh + fileh,'r')
fhn = open_ncfile(indirhn + filehn,'r')


# ----- Work/Variables ------

domains = ['Southern ST', 'SO', 'Northern ST', 'North Atlantic', 'North Pacific']

multStd = 2. # detect ToE at multStd std dev of histNat

varname = defVarmme('salinity'); v = 'S'
#varname = defVarmme('temp'); v = 'T'
#varname= defVarmme('depth'); v = 'Z'

var = varname['var_zonal']
legVar = varname['legVar']
unit = varname['unit']

# Time
iniyear = 1860
tmin = 0
tmax = 150

# Domains
SouthST = ToEdomainhistvshistNat('MME', domains[0])[0]
SO = ToEdomainhistvshistNat('MME', domains[1])[0]
NorthST = ToEdomainhistvshistNat('MME', domains[2])[0]
NA = ToEdomainhistvshistNat('MME', domains[3])[0]
NP = ToEdomainhistvshistNat('MME', domains[4])[0]


# Read variables
lat = fh.variables['latitude'][:]
density = fh.variables['lev'][:]
time = fh.variables['time'][:]

varh = fh.variables[var][:]
varhn = fhn.variables[var][:]

# Take difference hist - histNat
var_diff = varh - varhn
# Compute std of histNat
varstd = np.ma.std(varhn, axis=0)


# Average in domains
if SouthST['Atlantic'] != None :
    varSouthST = averageDom(var_diff[:,1,:,:], 3, SouthST['Atlantic'], lat, density)
    noiseSouthST = averageDom(2*varstd[1,:,:], 2, SouthST['Atlantic'], lat, density)
if SO['Atlantic'] != None :
    varSO = averageDom(var_diff[:,1,:,:], 3, SO['Atlantic'], lat, density)
    noiseSO = averageDom(2*varstd[1,:,:], 2, SO['Atlantic'], lat, density)
if NA['Atlantic'] != None :
    varNA = averageDom(var_diff[:,1,:,:], 3, NA['Atlantic'], lat, density)
    noiseNA = averageDom(2*varstd[1,:,:], 2, NA['Atlantic'], lat, density)


# ----- Plot ------

fig, ax = plt.subplots(figsize=(10,5))

ax.set_xlim([tmin+iniyear,tmax+iniyear])
# North Atlantic
ax.plot(time+iniyear, varNA, linewidth=2, color='blue', label=domains[3])
ax.plot([tmin+iniyear, tmax+iniyear], [noiseNA,noiseNA], ls='--', linewidth=2, color='blue')
ax.plot([tmin+iniyear, tmax+iniyear], [-noiseNA,-noiseNA], ls='--', linewidth=2, color='blue')
# Southern ST
ax.plot(time+iniyear, varSouthST, linewidth=2, color='green', label=domains[0])
ax.plot([tmin+iniyear, tmax+iniyear], [noiseSouthST,noiseSouthST], ls='--', linewidth=2, color='green')
ax.plot([tmin+iniyear, tmax+iniyear], [-noiseSouthST,-noiseSouthST], ls='--', linewidth=2, color='green')
# Southern Ocean
ax.plot(time+iniyear, varSO, linewidth=2, color='red', label=domains[1])
ax.plot([tmin+iniyear, tmax+iniyear], [noiseSO,noiseSO], ls='--', linewidth=2, color='red')
ax.plot([tmin+iniyear, tmax+iniyear], [-noiseSO,-noiseSO], ls='--', linewidth=2, color='red')


ax.axhline(y=0, color='black')
ax.set_xlabel('Years')
ax.set_ylabel('PSU')

ax.legend(loc='upper left', title='Atlantic', fontsize=11)

plt.title('Time of Emergence [>2std] mme_hist vs. mme_histNat (15 models) for %s' %(legVar,),
          va='bottom', fontweight='bold')

plotName = 'ToE_mme_histNat_'+legVar

plt.show()
#plt.savefig('/home/ysilvy/Density_bining/Yona_analysis/figures/models/ToE/histvshistNat/average_signal/'+plotName+'.png', bbox_inches='tight')

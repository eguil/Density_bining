#!/bin/env python
# -*- coding: utf-8 -*-
"""
Python matplotlib
Make time emergence plots from
 1) <var>agree variable in Historical runs
 2) Hist vs. histNat runs
 3) using piControl stdev (TODO)

(c) Eric Guilyardi March 2016

TODO: - add arguments for variable and output type

"""
from netCDF4 import Dataset as open_ncfile
import matplotlib.pyplot as plt
from densit_matplot_lib import defVar, averageDom
import numpy as np

# -------------------------------------------------------------------------------
#                                Define work
# -------------------------------------------------------------------------------

inDir = '/Users/ericg/Projets/Density_bining/'
workh = 'Prod_density_april15/mme_hist/'
workhn = 'Prod_density_april15/mme_histNat/'
inDirh = inDir + workh
inDirhn = inDir + workhn

file2dh  = 'cmip5.multimodel_Nat.historical.ensm.an.ocn.Omon.density_zon2D.nc'
file2dhn = 'cmip5.multimodel_All.historicalNat.ensm.an.ocn.Omon.density_zon2D.nc'


# Choose method
#ToE_Method = 'agreehist'
ToE_Method = 'usehistNat'
#ToE_Method = 'usepiControl'

# Define variable  TODO: read as argument
varname = defVar('salinity')
varname = defVar('temp')
varname = defVar('depth')
#varname = defVar('volume')
# varname = defVar('persist')

# Define plot name and title
plotName = 'cmip5_ToE_mme_hist_Nat_16models_' + varname['var']
figTitle = 'Time of emergence mme_hist (16 models) vs. HistNat ['+ varname['var']+']'

# density/latitude domains (where agreement level was noted > agreelev on ys plots)
# Atl
ToEA1 = {'domain': [-40., -30, 25.5, 26.7],  # [minlat, maxlat, minrho, maxrho]
         'name': 'Southern ST', 'color': 'blue'}
ToEA2 = {'domain': [20., 40., 26., 27.2], 'name': 'Northern ST', 'color': 'green'}
ToEA3 = {'domain': [-55, -40., 26.8, 27.8], 'name': 'Southern O', 'color': 'red'}
ToEA4 = {'domain': [5, 25, 23, 25.5], 'name': 'Tropics N', 'color': 'purple'}
ToEA5 = {'domain': [65, 90, 26.5, 28.], 'name': 'Arctic', 'color': 'black'}
# Pac
ToEP1 = {'domain': [-40, 0, 24.5, 26.7], 'name': 'Southern ST', 'color': 'blue'}
ToEP2 = {'domain': [-60, -45, 27.3, 27.7], 'name': 'Southern O', 'color': 'red'}
ToEP3 = {'domain': [45, 60, 26., 26.8], 'name': 'North Pac', 'color': 'green'}
ToEP4 = {'domain': [-30, -10, 25., 26.], 'name': 'Southern STz', 'color': 'purple'}
# Ind
ToEI1 = {'domain': [-40, -15, 25.6, 26.8], 'name': 'Southern ST', 'color': 'blue'}
ToEI2 = {'domain': [-60, -40, 27.2, 27.8], 'name': 'Southern O', 'color': 'red'}

# Plot domain (time)
iniyear = 1861
tmin0 = 0
tmax0 = 149

# Plot domain (vertical)
if ToE_Method == 'agreehist':
    # Model agreement level
    agreelev = 0.6
#
#
# -------------------------------------------------------------------------------

# -- Define variable properties

var = varname['var']
minmax = varname['1dminmax']
clevsm = varname['clevsm']
legVar = varname['legVar']
unit = varname['unit']

minmax = varname['1dminmax']
if ToE_Method == 'agreehist':
    minmax=[-1.,1.]
#
# -- Open netcdf files
nc2dh = open_ncfile(inDirh + '/' + file2dh)
if ToE_Method == 'usehistNat':
    nc2dhn = open_ncfile(inDirhn + '/' + file2dhn)

# -- Read variables
if ToE_Method == 'agreehist':
    # Read model agreement variables
    tvaraa = nc2dh.variables[var + 'Agree'][:, 1, :, :].squeeze()
    tvarap = nc2dh.variables[var + 'Agree'][:, 2, :, :].squeeze()
    tvarai = nc2dh.variables[var + 'Agree'][:, 3, :, :].squeeze()
if ToE_Method == 'usehistNat':
    print 'use histNat'
    # Read var in bowl for hist and histnat
    tvaraah = nc2dh.variables[var + 'Bowl'][:, 1, :, :].squeeze()
    tvaraph = nc2dh.variables[var + 'Bowl'][:, 2, :, :].squeeze()
    tvaraih = nc2dh.variables[var + 'Bowl'][:, 3, :, :].squeeze()
    tvaraahn = nc2dhn.variables[var + 'Bowl'][:, 1, :, :].squeeze()
    tvaraphn = nc2dhn.variables[var + 'Bowl'][:, 2, :, :].squeeze()
    tvaraihn = nc2dhn.variables[var + 'Bowl'][:, 3, :, :].squeeze()
    # Take difference Hist minus histNat
    tvaraa = tvaraah - tvaraahn
    tvarap = tvaraph - tvaraphn
    tvarai = tvaraih - tvaraihn
    # Compute stddev of histNat
    tstdahn = np.ma.std(tvaraahn, axis=0)
    tstdphn = np.ma.std(tvaraphn, axis=0)
    tstdihn = np.ma.std(tvaraihn, axis=0)

levr = nc2dh.variables['lev'][:]
lats = nc2dh.variables['latitude'][:]
time = nc2dh.variables['time'][:]

# -- Build plot variables (average in domain)
# Atl
varToEA1 = averageDom(tvaraa, 3, ToEA1['domain'], lats, levr)
varToEA2 = averageDom(tvaraa, 3, ToEA2['domain'], lats, levr)
varToEA3 = averageDom(tvaraa, 3, ToEA3['domain'], lats, levr)
varToEA4 = averageDom(tvaraa, 3, ToEA4['domain'], lats, levr)
varToEA5 = averageDom(tvaraa, 3, ToEA5['domain'], lats, levr)
# Pac
varToEP1 = averageDom(tvarap, 3, ToEP1['domain'], lats, levr)
varToEP2 = averageDom(tvarap, 3, ToEP2['domain'], lats, levr)
varToEP3 = averageDom(tvarap, 3, ToEP3['domain'], lats, levr)
varToEP4 = averageDom(tvarap, 3, ToEP4['domain'], lats, levr)
# Ind
varToEI1 = averageDom(tvarai, 3, ToEI1['domain'], lats, levr)
varToEI2 = averageDom(tvarai, 3, ToEI2['domain'], lats, levr)

# -- Define detection level
if ToE_Method == 'agreehist':
    # set to agreelev
    varToEA1d = agreelev
    varToEA2d = agreelev
    varToEA3d = agreelev
    varToEA4d = agreelev
    varToEA5d = agreelev
    varToEP1d = agreelev
    varToEP2d = agreelev
    varToEP3d = agreelev
    varToEP4d = agreelev
    varToEI1d = agreelev
    varToEI2d = agreelev
if ToE_Method == 'usehistNat':
    # set to stddev of region
    varToEA1d = averageDom(tstdahn, 2, ToEA1['domain'], lats, levr)
    varToEA2d = averageDom(tstdahn, 2, ToEA2['domain'], lats, levr)
    varToEA3d = averageDom(tstdahn, 2, ToEA3['domain'], lats, levr)
    varToEA4d = averageDom(tstdahn, 2, ToEA4['domain'], lats, levr)
    varToEA5d = averageDom(tstdahn, 2, ToEA5['domain'], lats, levr)
    varToEP1d = averageDom(tstdphn, 2, ToEP1['domain'], lats, levr)
    varToEP2d = averageDom(tstdphn, 2, ToEP2['domain'], lats, levr)
    varToEP3d = averageDom(tstdphn, 2, ToEP3['domain'], lats, levr)
    varToEP4d = averageDom(tstdphn, 2, ToEP4['domain'], lats, levr)
    varToEI1d = averageDom(tstdihn, 2, ToEI1['domain'], lats, levr)
    varToEI2d = averageDom(tstdihn, 2, ToEI2['domain'], lats, levr)

#
# -------- Make plot ----------------
#
tmin = tmin0 + iniyear
tmax = tmax0 + iniyear

# -- Create figure and axes instances
fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(10, 10))

# Atlantic

axes[0].axis([tmin, tmax, minmax[0], minmax[1]])

xaxisa0 = axes[0].plot([tmin, tmax], [0., 0.], linestyle='-', color='black', linewidth=1)

print time.shape, varToEA1.shape
pltA1 = axes[0].plot(time + iniyear, varToEA1, linestyle='-', color=ToEA1['color'], linewidth=2, label=ToEA1['name'])
xaxisap = axes[0].plot([tmin, tmax], [varToEA1d, varToEA1d], linestyle='--', color=ToEA1['color'], linewidth=1)
xaxisam = axes[0].plot([tmin, tmax], [-varToEA1d, -varToEA1d], linestyle='--', color=ToEA1['color'], linewidth=1)

pltA2 = axes[0].plot(time + iniyear, varToEA2, linestyle='-', color=ToEA2['color'], linewidth=2, label=ToEA2['name'])
xaxisap = axes[0].plot([tmin, tmax], [varToEA2d, varToEA2d], linestyle='--', color=ToEA2['color'], linewidth=1)
xaxisam = axes[0].plot([tmin, tmax], [-varToEA2d, -varToEA2d], linestyle='--', color=ToEA2['color'], linewidth=1)
pltA3 = axes[0].plot(time + iniyear, varToEA3, linestyle='-', color=ToEA3['color'], linewidth=2, label=ToEA3['name'])
xaxisap = axes[0].plot([tmin, tmax], [varToEA3d, varToEA3d], linestyle='--', color=ToEA3['color'], linewidth=1)
xaxisam = axes[0].plot([tmin, tmax], [-varToEA3d, -varToEA3d], linestyle='--', color=ToEA3['color'], linewidth=1)
pltA4 = axes[0].plot(time + iniyear, varToEA4, linestyle='-', color=ToEA4['color'], linewidth=2, label=ToEA4['name'])
xaxisap = axes[0].plot([tmin, tmax], [varToEA4d, varToEA4d], linestyle='--', color=ToEA4['color'], linewidth=1)
xaxisam = axes[0].plot([tmin, tmax], [-varToEA4d, -varToEA4d], linestyle='--', color=ToEA4['color'], linewidth=1)
pltA5 = axes[0].plot(time + iniyear, varToEA5, linestyle='-', color=ToEA5['color'], linewidth=2, label=ToEA5['name'])
xaxisap = axes[0].plot([tmin, tmax], [varToEA5d, varToEA5d], linestyle='--', color=ToEA5['color'], linewidth=1)
xaxisam = axes[0].plot([tmin, tmax], [-varToEA5d, -varToEA5d], linestyle='--', color=ToEA5['color'], linewidth=1)
axes[0].legend(loc='upper left', title='Atlantic', prop={'size': 10})

# Pac

axes[1].axis([tmin, tmax, minmax[0], minmax[1]])
xaxisp0 = axes[1].plot([tmin, tmax], [0., 0.], linestyle='-', color='black', linewidth=1)
#xaxispp = axes[1].plot([tmin, tmax], [agreelev, agreelev], linestyle='--', color='red', linewidth=1)
#xaxispm = axes[1].plot([tmin, tmax], [-agreelev, -agreelev], linestyle='--', color='blue', linewidth=1)

pltP1 = axes[1].plot(time + iniyear, varToEP1, linestyle='-', color=ToEP1['color'], linewidth=2, label=ToEP1['name'])
xaxisap = axes[1].plot([tmin, tmax], [varToEP1d, varToEP1d], linestyle='--', color=ToEP1['color'], linewidth=1)
xaxisam = axes[1].plot([tmin, tmax], [-varToEP1d, -varToEP1d], linestyle='--', color=ToEP1['color'], linewidth=1)
pltP2 = axes[1].plot(time + iniyear, varToEP2, linestyle='-', color=ToEP2['color'], linewidth=2, label=ToEP2['name'])
xaxisap = axes[1].plot([tmin, tmax], [ varToEP2d,  varToEP2d], linestyle='--', color=ToEP2['color'], linewidth=1)
xaxisam = axes[1].plot([tmin, tmax], [-varToEP2d, -varToEP2d], linestyle='--', color=ToEP2['color'], linewidth=1)
pltP3 = axes[1].plot(time + iniyear, varToEP3, linestyle='-', color=ToEP3['color'], linewidth=2, label=ToEP3['name'])
xaxisap = axes[1].plot([tmin, tmax], [ varToEP3d,  varToEP3d], linestyle='--', color=ToEP3['color'], linewidth=1)
xaxisam = axes[1].plot([tmin, tmax], [-varToEP3d, -varToEP3d], linestyle='--', color=ToEP3['color'], linewidth=1)
pltP4 = axes[1].plot(time + iniyear, varToEP4, linestyle='-', color=ToEP4['color'], linewidth=2, label=ToEP4['name'])
xaxisap = axes[1].plot([tmin, tmax], [ varToEP4d,  varToEP4d], linestyle='--', color=ToEP4['color'], linewidth=1)
xaxisam = axes[1].plot([tmin, tmax], [-varToEP4d, -varToEP4d], linestyle='--', color=ToEP4['color'], linewidth=1)
axes[1].legend(loc='upper left', title='Pacific', prop={'size': 10})

# Ind

axes[2].axis([tmin, tmax, minmax[0], minmax[1]])
xaxisi0 = axes[2].plot([tmin, tmax], [0., 0.], linestyle='-', color='black', linewidth=1)
#xaxisip = axes[2].plot([tmin, tmax], [agreelev, agreelev], linestyle='--', color='red', linewidth=1)
#xaxisim = axes[2].plot([tmin, tmax], [-agreelev, -agreelev], linestyle='--', color='blue', linewidth=1)

pltI1 = axes[2].plot(time + iniyear, varToEI1, linestyle='-', color=ToEI1['color'], linewidth=2, label=ToEI1['name'])
xaxisap = axes[2].plot([tmin, tmax], [varToEI1d, varToEI1d], linestyle='--', color=ToEI1['color'], linewidth=1)
xaxisam = axes[2].plot([tmin, tmax], [-varToEI1d, -varToEI1d], linestyle='--', color=ToEI1['color'], linewidth=1)
pltI2 = axes[2].plot(time + iniyear, varToEI2, linestyle='-', color=ToEI2['color'], linewidth=2, label=ToEI2['name'])
xaxisap = axes[2].plot([tmin, tmax], [varToEI2d, varToEI2d], linestyle='--', color=ToEI2['color'], linewidth=1)
xaxisam = axes[2].plot([tmin, tmax], [-varToEI2d, -varToEI2d], linestyle='--', color=ToEI2['color'], linewidth=1)

# -- Add legend for bowl position
axes[2].legend(loc='upper left', title='Indian', prop={'size': 10})

# -- add plot title
plt.suptitle(figTitle, fontsize=14, fontweight='bold')

# -- Output  # TODO read as argument

#plt.show()
plt.savefig(plotName + '.pdf', bbox_inches='tight')

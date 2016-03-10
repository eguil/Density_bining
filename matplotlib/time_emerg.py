#!/bin/env python
# -*- coding: utf-8 -*-
"""
Python matplotlib
Make time emergence plots from <var>agree variable in Historicla runs

(c) Eric Guilyardi March 2016

TODO: - add arguments for variable and output type

"""
import numpy as np
from   netCDF4 import Dataset as open_ncfile
import matplotlib.pyplot as plt
from densitlib import defVar, averageDom

# -------------------------------------------------------------------------------
#                                Define work
# -------------------------------------------------------------------------------

indir = '/Users/ericg/Projets/Density_bining/'
work = 'Prod_density_april15/mme_hist/'
indir = indir + work
file2d = 'cmip5.multimodel.historical.ensm.an.ocn.Omon.density_zon2D.nc'


# Model agreement level
agreelev = 0.6

# Define variable  TODO: read as argument
varname = defVar('salinity')
# varname = defVar('temp')
# varname = defVar('depth')
# varname = defVar('volume')
# varname = defVar('persist')

# Define plot name
plotName = 'cmip5_ToE_mme_hist_r1i1p1_' + varname['var']

# density/latitude domains (where agreement level was noted > agreelev on ys plots)
# Atl
ToEA1 = {'domain': [-40., -20, 26., 26.7],  # [minlat, maxlat, minrho, maxrho]
         'name': 'Southern ST', 'color': 'blue'}
ToEA2 = {'domain': [20., 40., 26., 27.2], 'name': 'Northern ST', 'color': 'green'}
ToEA3 = {'domain': [-55, -40., 26.8, 27.8], 'name': 'Southern O', 'color': 'red'}
ToEA4 = {'domain': [5, 25, 23.5, 25.5], 'name': 'Tropics N', 'color': 'purple'}
ToEA5 = {'domain': [65, 90, 26.4, 28.], 'name': 'Arctic', 'color': 'black'}
# Pac
ToEP1 = {'domain': [-45, 0, 24.5, 26.8], 'name': 'Southern ST', 'color': 'blue'}
ToEP2 = {'domain': [-60, -45, 27.2, 27.8], 'name': 'Southern O', 'color': 'red'}
ToEP3 = {'domain': [45, 60, 26.4, 27.], 'name': 'North Pac', 'color': 'green'}
ToEP4 = {'domain': [-30, -20, 25., 26.1], 'name': 'Southern STz', 'color': 'blue'}
# Ind
ToEI1 = {'domain': [-40, -20, 25.8, 26.6], 'name': 'Southern ST', 'color': 'blue'}
ToEI2 = {'domain': [-55, -45, 27., 27.8], 'name': 'Southern O', 'color': 'red'}

# Plot domain
iniyear = 1860
tmin0 = 0
tmax0 = 150
agreemin = -1.
agreemax = 1.
#
figTitle = "Time of emergence mme_hist for agreelev = ",agreelev
#
# -------------------------------------------------------------------------------

# -- Define variable properties

var = varname['var']

# -- Open netcdf files
nc2d = open_ncfile(indir + '/' + file2d)
# -- Read variables
# Read model agreement variables
tvaraa = nc2d.variables[var + 'Agree'][:, 1, :, :].squeeze()
tvarap = nc2d.variables[var + 'Agree'][:, 2, :, :].squeeze()
tvarai = nc2d.variables[var + 'Agree'][:, 3, :, :].squeeze()
levr = nc2d.variables['lev'][:]
lats = nc2d.variables['latitude'][:]
time = nc2d.variables['time'][:]

# -- Build plot variables (average in domain)
# Atl
varToEA1 = averageDom(tvaraa, ToEA1['domain'], lats, levr)
varToEA2 = averageDom(tvaraa, ToEA2['domain'], lats, levr)
varToEA3 = averageDom(tvaraa, ToEA3['domain'], lats, levr)
varToEA4 = averageDom(tvaraa, ToEA4['domain'], lats, levr)
varToEA5 = averageDom(tvaraa, ToEA5['domain'], lats, levr)
# Pac
varToEP1 = averageDom(tvarap, ToEP1['domain'], lats, levr)
varToEP2 = averageDom(tvarap, ToEP2['domain'], lats, levr)
varToEP3 = averageDom(tvarap, ToEP3['domain'], lats, levr)
varToEP4 = averageDom(tvarap, ToEP4['domain'], lats, levr)
# Ind
varToEI1 = averageDom(tvarai, ToEI1['domain'], lats, levr)
varToEI2 = averageDom(tvarai, ToEI2['domain'], lats, levr)

#
# -------- Make plot ----------------
#
tmin = tmin0 + iniyear
tmax = tmax0 + iniyear

# -- Create figure and axes instances
fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(10, 10))

# Atlantic

axes[0].axis([tmin, tmax, agreemin, agreemax])

xaxisa0 = axes[0].plot([tmin, tmax], [0., 0.], linestyle='--', color='black', linewidth=1)
xaxisap = axes[0].plot([tmin, tmax], [agreelev, agreelev], linestyle='--', color='red', linewidth=1)
xaxisam = axes[0].plot([tmin, tmax], [-agreelev, -agreelev], linestyle='--', color='blue', linewidth=1)

pltA1 = axes[0].plot(time + iniyear, varToEA1, linestyle='-', color=ToEA1['color'], linewidth=2, label=ToEA1['name'])
pltA2 = axes[0].plot(time + iniyear, varToEA2, linestyle='-', color=ToEA2['color'], linewidth=2, label=ToEA2['name'])
pltA3 = axes[0].plot(time + iniyear, varToEA3, linestyle='-', color=ToEA3['color'], linewidth=2, label=ToEA3['name'])
pltA4 = axes[0].plot(time + iniyear, varToEA4, linestyle='-', color=ToEA4['color'], linewidth=2, label=ToEA4['name'])
pltA5 = axes[0].plot(time + iniyear, varToEA5, linestyle='-', color=ToEA5['color'], linewidth=2, label=ToEA5['name'])
axes[0].legend(loc='upper left', title='Atlantic', prop={'size': 10})

# Pac

axes[1].axis([tmin, tmax, agreemin, agreemax])
xaxisp0 = axes[1].plot([tmin, tmax], [0., 0.], linestyle='--', color='black', linewidth=1)
xaxispp = axes[1].plot([tmin, tmax], [agreelev, agreelev], linestyle='--', color='red', linewidth=1)
xaxispm = axes[1].plot([tmin, tmax], [-agreelev, -agreelev], linestyle='--', color='blue', linewidth=1)

pltP1 = axes[1].plot(time + iniyear, varToEP1, linestyle='-', color=ToEP1['color'], linewidth=2, label=ToEP1['name'])
pltP2 = axes[1].plot(time + iniyear, varToEP2, linestyle='-', color=ToEP2['color'], linewidth=2, label=ToEP2['name'])
pltP3 = axes[1].plot(time + iniyear, varToEP3, linestyle='-', color=ToEP3['color'], linewidth=2, label=ToEP3['name'])
pltP4 = axes[1].plot(time + iniyear, varToEP4, linestyle='--', color=ToEP4['color'], linewidth=2, label=ToEP4['name'])
axes[1].legend(loc='upper left', title='Pacific', prop={'size': 10})

# Ind

axes[2].axis([tmin, tmax, agreemin, agreemax])
xaxisi0 = axes[2].plot([tmin, tmax], [0., 0.], linestyle='--', color='black', linewidth=1)
xaxisip = axes[2].plot([tmin, tmax], [agreelev, agreelev], linestyle='--', color='red', linewidth=1)
xaxisim = axes[2].plot([tmin, tmax], [-agreelev, -agreelev], linestyle='--', color='blue', linewidth=1)

pltI1 = axes[2].plot(time + iniyear, varToEI1, linestyle='-', color=ToEI1['color'], linewidth=2, label=ToEI1['name'])
pltI2 = axes[2].plot(time + iniyear, varToEI2, linestyle='-', color=ToEI2['color'], linewidth=2, label=ToEI2['name'])

# -- Add legend for bowl position
axes[2].legend(loc='upper left', title='Indian', prop={'size': 10})

# -- add plot title
plt.suptitle(figTitle, fontsize=14, fontweight='bold')

# -- Output  # TODO read as argument

#plt.show()
plt.savefig(plotName+'.pdf', bbox_inches='tight')

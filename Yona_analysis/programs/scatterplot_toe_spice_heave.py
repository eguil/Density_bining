#!/bin/env python
# -*- coding: utf-8 -*-

"""
Python matplotlib
Make scatter plots of ToE Salinity (spice) / Depth (heave) hist+RCP85 vs. histNat in all 5 domains
"""

import numpy as np
from numpy.polynomial.polynomial import polyfit
import matplotlib.pyplot as plt
from netCDF4 import Dataset as open_ncfile
from maps_matplot_lib import defVarmme
from modelsDef import defModels, modelcolors
from matplotlib.ticker import MultipleLocator, AutoMinorLocator
import glob
import os
import datetime

# ----- Work -----

# Directory
indir_rcphn = '/home/ysilvy/Density_bining/Yona_analysis/data/toe_rcp85_histNat_average_signal/'
indir_rcppiC = '/home/ysilvy/Density_bining/Yona_analysis/data/toe_rcp85_PiControl_average_signal/'
indir_noise = '/home/ysilvy/Density_bining/Yona_analysis/data/noise_estimate/RCP85vshistNat_domains/'

models = defModels()

domains = ['Southern ST', 'SO', 'Northern ST', 'North Atlantic', 'North Pacific']

varname2 = defVarmme('salinity'); v = 'S'
varname1 = defVarmme('depth'); v = 'Z'

method = 'average_signal' # Average signal and noise in the box, then compute ToE

# === INPUTS ===

method_noise_rcp = 'average_histNat' # Average histNat in the specified domains then determine the std of this averaged value

use_piC = False # signal = hist+RCP - histNat, noise = std(histNat)
# use_piC = True # signal = hist+RCP - PiControl, noise = std(PiControl)

# output format
# outfmt = 'view'
outfmt = 'save'

# ===========

# ----- Variables ------
var1 = varname1['var_zonal_w/bowl']
legVar1 = varname1['legVar']
unit1 = varname1['unit']

var2 = varname2['var_zonal_w/bowl']
legVar2 = varname2['legVar']
unit2 = varname2['unit']

# ----- Read ToE for each model ------

# == Historical + RCP8.5 vs. historicalNat or vs. PiControl ==

nruns = 0 # Initialize total number of runs
nrunmax = 100

# -- Initialize varToE 1 and 2 containing ToE of all runs
varToEA1 = np.ma.masked_all((nrunmax, len(domains)))
varToEP1 = np.ma.masked_all((nrunmax, len(domains)))
varToEI1 = np.ma.masked_all((nrunmax, len(domains)))
varToEA2 = np.ma.masked_all((nrunmax, len(domains)))
varToEP2 = np.ma.masked_all((nrunmax, len(domains)))
varToEI2 = np.ma.masked_all((nrunmax, len(domains)))

# -- Loop over models

listfiles1 = sorted(glob.glob(indir_rcphn + method_noise_rcp + '/*Depth*.nc'))
listfiles2 = sorted(glob.glob(indir_rcphn + method_noise_rcp + '/*.toe*.nc'))

nmodels = len(listfiles2)
nMembers = np.ma.zeros(len(models)) # Initialize array for keeping number of members per model
names = ['']*nmodels # Initialize array for keeping names of models

for i in range(nmodels):

    # Depth
    file_toe1 = listfiles1[i]
    ftoe1 = open_ncfile(file_toe1, 'r')
    
    # Spice
    file_toe2 = listfiles2[i]
    ftoe2 = open_ncfile(file_toe2, 'r')
    
    name = os.path.basename(file_toe1).split('.')[1]
    name_check = os.path.basename(file_toe2).split('.')[1]
    
    if name != name_check:
        print('ERROR: NOT READING SAME MODEL FOR HEAVE AND SPICE')

    # Read ToE (members, basin, domain)
    toe1read = ftoe1.variables[var1 + 'ToE2'][:]
    toe2read = ftoe2.variables[var2 + 'ToE2'][:]
    nMembers[i] = toe2read.shape[0]
    print('- Reading ToE of %s with %d members'%(name,nMembers[i]))
    nruns1 = nruns + nMembers[i]

    # Save ToE
    varToEA1[nruns:nruns1,:] = toe1read[:,1,:]
    varToEP1[nruns:nruns1,:] = toe1read[:,2,:]
    varToEI1[nruns:nruns1,:] = toe1read[:,3,:]
    varToEA2[nruns:nruns1,:] = toe2read[:,1,:]
    varToEP2[nruns:nruns1,:] = toe2read[:,2,:]
    varToEI2[nruns:nruns1,:] = toe2read[:,3,:]

    nruns = nruns1

    names[i] = name

print('Total number of runs:', nruns)
varToEA1 = varToEA1[0:nruns,:]
varToEP1 = varToEP1[0:nruns,:]
varToEI1 = varToEI1[0:nruns,:]
varToEA2 = varToEA2[0:nruns,:]
varToEP2 = varToEP2[0:nruns,:]
varToEI2 = varToEI2[0:nruns,:]

nruns = int(nruns)
cumMembers = np.cumsum(nMembers)

# ----- Turn masked data into nans -----

varToEA1[np.ma.getmask(varToEA1)] = np.nan
varToEP1[np.ma.getmask(varToEP1)] = np.nan
varToEI1[np.ma.getmask(varToEI1)] = np.nan
varToEA2[np.ma.getmask(varToEA2)] = np.nan
varToEP2[np.ma.getmask(varToEP2)] = np.nan
varToEI2[np.ma.getmask(varToEI2)] = np.nan

# ----- Plot ------

maskdata = np.nan

dataToE1 = np.array([varToEA1[:,0], varToEP1[:,0], varToEI1[:,0], varToEP1[:,2], varToEI1[:,2],
          varToEA1[:,1], varToEP1[:,1], varToEI1[:,1], varToEA1[:,3], varToEP1[:,4]])

dataToE2 = np.array([varToEA2[:,0], varToEP2[:,0], varToEI2[:,0], varToEP2[:,2], varToEI2[:,2],
          varToEA2[:,1], varToEP2[:,1], varToEI2[:,1], varToEA2[:,3], varToEP2[:,4]])


# ===== One plot per domain ======

# Create figure
fig, axes = plt.subplots(nrows=2,ncols=5,sharex=True,sharey=True,figsize=(14.5,6))

plottitles = ['Southern ST, Atlantic','Southern ST, Pacific', 'Southern ST, Indian',
              'Northern ST, Pacific', 'Northern ST, Indian',
              'South. Ocean, Atlantic','South. Ocean, Pacific', 'South. Ocean, Indian',
              'North Atlantic',
              'North Pacific']

ax1D = axes.ravel().tolist()
l=['']*nmodels

for idomain in range(len(ax1D)):

    ax = ax1D[idomain]
    l[0] = ax.scatter(dataToE1[idomain,0:cumMembers[0]],dataToE2[idomain,0:cumMembers[0]],color=modelcolors(names[0])['color'],
           s=20,facecolors='none',marker=modelcolors(names[0])['marker'], label=names[0])

    for i in range(1,nmodels):
        col = modelcolors(names[i])['color']
        mk = modelcolors(names[i])['marker']
        l[i] = ax.scatter(dataToE1[idomain,cumMembers[i-1]:cumMembers[i]],dataToE2[idomain,cumMembers[i-1]:cumMembers[i]],
                   s=20, color=col, facecolors='none', marker=mk, label=names[i])

    ax.plot([1860,2000,2100],[1860,2000,2100], color='black', linestyle='--')

    ax.set_title(plottitles[idomain],fontsize=11,fontweight='bold')

    for tick in ax.get_xticklabels():
        tick.set_rotation(35)

ax1D[7].set_xlabel('Time of Emergence ['+legVar1+']', fontweight='bold',fontsize=13,va='baseline')

title = 'Hist+RCP8.5 vs. histNat ('+str(nruns)+' runs)'
noise = 'histNat'
plotName = 'Scatterplot_ToE_spice_heave_RCP85vshistNat'
  

ax.set_ylim([1860,2105])
ax.set_xlim([1860,2105])
majorLocator = MultipleLocator(40)
minorLocator = AutoMinorLocator(2)
ax.yaxis.set_major_locator(majorLocator)
ax.yaxis.set_minor_locator(minorLocator)
ax.xaxis.set_major_locator(majorLocator)
ax.xaxis.set_minor_locator(minorLocator)

plt.suptitle('ToE['+legVar2+']/ToE['+legVar1+'] for '+title, fontweight='bold', fontsize=14)

# plt.figtext(0.3,0.033,'Time of Emergence [>1std]', fontweight='bold',fontsize=13)
plt.figtext(0.007,0.7,'Time of Emergence ['+legVar2+']', fontweight='bold', rotation='vertical',fontsize=13)

# Put a legend to the right of the current axis
lgd = fig.legend(l,names,loc='center right',scatterpoints=1,frameon=False, markerscale=2)

plt.subplots_adjust(left=0.05,right=0.82,hspace=0.13,bottom=0.12,wspace=0.09)

# Date
now = datetime.datetime.now()
date = now.strftime("%Y-%m-%d")

# Text at the bottom of the figure
plt.figtext(.8,.01,'Computed by : scatterplot_toe_spice_heave.py, '+date, fontsize=8, ha='center')
plt.figtext(.1,.01,'Noise: %s' %(method_noise_rcp), fontsize=8, ha='center')

if outfmt == 'view':
    plt.show()
else:
    dir = '/home/ysilvy/figures/models/ToE/scatterplots/'
    plt.savefig(dir+plotName+'.pdf',bbox_extra_artists=(lgd,))

#!/bin/env python
# -*- coding: utf-8 -*-

"""
Python matplotlib
Make scatter plots of ToE1 [>1std] / ToE2[>2std] hist+RCP85 vs. histNat (or Picontrol) in all 5 domains
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

varname = defVarmme('salinity'); v = 'S'

method = 'average_signal' # Average signal and noise in the box, then compute ToE

# === INPUTS ===

# method_noise_rcphn = 'average_std' # Average the standard deviation of histNat in the specified domains
method_noise_rcp = 'average_histNat' # Average histNat in the specified domains then determine the std of this averaged value

use_piC = False # signal = hist+RCP - histNat, noise = std(histNat)
# use_piC = True # signal = hist+RCP - PiControl, noise = std(PiControl)

if use_piC == True and method_noise_rcp == 'average_histNat':
    method_noise_rcp = 'average_piC' # Just change the notation for coherence and for data path

# runs_rcp = 'same' # Same runs (30 runs) for hist+RCP8.5 vs. histNat as for hist+RCP8.5 vs. PiControl (use only with histNat method)
runs_rcp = 'all' # All runs (35)

# output format
# outfmt = 'view'
outfmt = 'save'

# ===========

# ----- Variables ------
var = varname['var_zonal_w/bowl']
legVar = varname['legVar']
unit = varname['unit']

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
if use_piC == True:
    indir = indir_rcppiC
else:
    indir = indir_rcphn
listfiles = glob.glob(indir + method_noise_rcp + '/*.nc')
listfiles = sorted(listfiles)
nmodels = len(listfiles)
nMembers = np.ma.zeros(len(models)) # Initialize array for keeping number of members per model
names = ['']*nmodels # Initialize array for keeping names of models

for i in range(nmodels):

    file_toe = listfiles[i]
    ftoe = open_ncfile(file_toe, 'r')
    name = os.path.basename(file_toe).split('.')[1]

    # If use same runs in vs. histNat as in vs. PiControl, take out deficient models
    if (runs_rcp == 'all') or (runs_rcp =='same' and name != 'GISS-E2-R' and name != 'FGOALS-g2' and name != 'MIROC-ESM'):

        # Read ToE (members, basin, domain)
        toe1read = ftoe.variables[var + 'ToE1'][:]
        toe2read = ftoe.variables[var + 'ToE2'][:]
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

if use_piC == False and runs_rcp == 'same':
        nmodels=nmodels-3
        names = filter(None,names)
        nMembers = nMembers[np.ma.nonzero(nMembers)]

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
fig, axes = plt.subplots(nrows=2,ncols=5,sharex=True,sharey=True,figsize=(14.5,5.5))

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

    ax.set_title(plottitles[idomain],fontsize=12,fontweight='bold')

    for tick in ax.get_xticklabels():
        tick.set_rotation(45)

ax1D[7].set_xlabel('Time of Emergence [>1std]', fontweight='bold',fontsize=13,va='baseline')

if use_piC:
    title = 'Hist+RCP8.5 vs. PiControl ('+str(nruns)+' runs)'
    plotName = 'Scatterplot_ToE1_ToE2_RCP85vsPiControl'
    noise = 'PiControl'
else:
    title = 'Hist+RCP8.5 vs. histNat ('+str(nruns)+' runs)'
    noise = 'histNat'
    if runs_rcp == 'all':
        plotName = 'Scatterplot_ToE1_ToE2_RCP85vshistNat'
    else:
        plotName = 'Scatterplot_ToE1_ToE2_RCP85vshistNat_samerunsvsPiC'


ax.set_ylim([1860,2105])
ax.set_xlim([1860,2105])
majorLocator = MultipleLocator(40)
minorLocator = AutoMinorLocator(2)
ax.yaxis.set_major_locator(majorLocator)
ax.yaxis.set_minor_locator(minorLocator)
ax.xaxis.set_major_locator(majorLocator)
ax.xaxis.set_minor_locator(minorLocator)

plt.suptitle('ToE[>2std]/ToE[>1std] for '+title+ ' in different regions', fontweight='bold', fontsize=14)

# plt.figtext(0.3,0.033,'Time of Emergence [>1std]', fontweight='bold',fontsize=13)
plt.figtext(0.007,0.7,'Time of Emergence [>2std]', fontweight='bold', rotation='vertical',fontsize=13)

# Put a legend to the right of the current axis
lgd = fig.legend(l,names,loc='center right',scatterpoints=1,frameon=False, markerscale=2)

plt.subplots_adjust(left=0.05,right=0.82,hspace=0.13,bottom=0.12,wspace=0.09)

# Date
now = datetime.datetime.now()
date = now.strftime("%Y-%m-%d")

# Text at the bottom of the figure
plt.figtext(.8,.01,'Computed by : scatterplot_toe1_toe2.py, '+date, fontsize=8, ha='center')
plt.figtext(.1,.01,'Noise: %s' %(method_noise_rcp), fontsize=8, ha='center')
if use_piC :
    plt.figtext(.25,.01,'PiControl : mean(last_240_years)',fontsize=8,ha='center')


if outfmt == 'view':
    plt.show()
else:
    dir = '/home/ysilvy/Density_bining/Yona_analysis/figures/models/ToE/scatterplots/ToE1_ToE2/'
    plt.savefig(dir+plotName+'.png',bbox_extra_artists=(lgd,))

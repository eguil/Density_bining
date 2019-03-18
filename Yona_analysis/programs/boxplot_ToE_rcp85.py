#!/bin/env python
# -*- coding: utf-8 -*-

"""
Python matplotlib
Make whisker plots of ToE hist+RCP85 vs. histNat (or Picontrol) in all 5 domains
>1std and >2std
"""

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset as open_ncfile
from maps_matplot_lib import defVarmme
from modelsDef import defModels
from matplotlib.ticker import AutoMinorLocator, MultipleLocator
import glob
import os
import datetime

# ----- Work -----

# Directory
indir_rcphn = '/home/ysilvy/Density_bining/Yona_analysis/data/toe_rcp85_histNat_average_signal/'
indir_rcppiC = '/home/ysilvy/Density_bining/Yona_analysis/data/toe_rcp85_PiControl_average_signal/'

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

# runs_rcp = 'same' # Same runs (30 runs) for hist+RCP8.5 vs. histNat as for hist+RCP8.5 vs. PiControl
runs_rcp = 'all' # All runs (35)

# output format
outfmt = 'view'
# outfmt = 'save'

# ===========

# ----- Variables ------
var = varname['var_zonal_w/bowl']
legVar = varname['legVar']
unit = varname['unit']


# ----- Read ToE for each model ------

# == Historical + RCP8.5 vs. historicalNat or vs. PiControl ==

nruns = 0 # Initialize total number of runs
nrunmax = 100
nMembers = np.ma.empty(len(models)) # Initialize array for keeping number of members per model

# -- Initialize varToE 1 and 2 containing ToE of all runs
varToEA_1 = np.ma.masked_all((nrunmax, len(domains)))
varToEP_1 = np.ma.masked_all((nrunmax, len(domains)))
varToEI_1 = np.ma.masked_all((nrunmax, len(domains)))
varToEA_2 = np.ma.masked_all((nrunmax, len(domains)))
varToEP_2 = np.ma.masked_all((nrunmax, len(domains)))
varToEI_2 = np.ma.masked_all((nrunmax, len(domains)))

# -- Loop over models
if use_piC == True:
    indir = indir_rcppiC
else:
    indir = indir_rcphn
listfiles = glob.glob(indir + method_noise_rcp + '/*.nc')
nmodels = len(listfiles)

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
        varToEA_1[nruns:nruns1,:] = toe1read[:,1,:]
        varToEP_1[nruns:nruns1,:] = toe1read[:,2,:]
        varToEI_1[nruns:nruns1,:] = toe1read[:,3,:]
        varToEA_2[nruns:nruns1,:] = toe2read[:,1,:]
        varToEP_2[nruns:nruns1,:] = toe2read[:,2,:]
        varToEI_2[nruns:nruns1,:] = toe2read[:,3,:]

        nruns = nruns1


print('Total number of runs:', nruns)
varToEA_1 = varToEA_1[0:nruns,:]
varToEP_1 = varToEP_1[0:nruns,:]
varToEI_1 = varToEI_1[0:nruns,:]
varToEA_2 = varToEA_2[0:nruns,:]
varToEP_2 = varToEP_2[0:nruns,:]
varToEI_2 = varToEI_2[0:nruns,:]

nruns = int(nruns)

if runs_rcp == 'same':
        nmodels=nmodels-3

# ----- Turn masked data into nans -----

varToEA_1[np.ma.getmask(varToEA_1)] = np.nan
varToEP_1[np.ma.getmask(varToEP_1)] = np.nan
varToEI_1[np.ma.getmask(varToEI_1)] = np.nan
varToEA_2[np.ma.getmask(varToEA_2)] = np.nan
varToEP_2[np.ma.getmask(varToEP_2)] = np.nan
varToEI_2[np.ma.getmask(varToEI_2)] = np.nan

# ----- Plot ------

maskdata  = np.nan

# ToE 1 [>1std] hist+rcp8.5 vs. histNat (or vs. PiControl)
data1 = [varToEA_1[:,0], varToEP_1[:,0], varToEI_1[:,0], maskdata, varToEP_1[:,2], varToEI_1[:,2], maskdata,
          varToEA_1[:,1], varToEP_1[:,1], varToEI_1[:,1], maskdata, varToEA_1[:,3], maskdata, varToEP_1[:,4]]
data1 = data1[::-1]

# ToE 2 [>2std]
data2 = [varToEA_2[:,0], varToEP_2[:,0], varToEI_2[:,0], maskdata, varToEP_2[:,2], varToEI_2[:,2], maskdata,
          varToEA_2[:,1], varToEP_2[:,1], varToEI_2[:,1], maskdata, varToEA_2[:,3], maskdata, varToEP_2[:,4]]
data2 = data2[::-1]


labels = ['','','','','Indian','Pacific','Atlantic','','Indian','Pacific','','Indian','Pacific','Atlantic']
N = 15
ind = np.arange(1,N)
width = 0.25

fig, ax = plt.subplots(figsize=(10,12))

ax.axvline(x=2005, color='black', ls=':')

# ToE 2 Hist+rcp8.5 vs. HistNat (or vs. PiControl) boxes
boxes1 = ax.boxplot(data2, vert=0, positions=ind-width, widths=width, whis=0)
for box in boxes1['boxes']:
    box.set(color='#64000b', linewidth=2) #c90016
for whisker in boxes1['whiskers']:
    whisker.set(color='#64000b', linestyle='-', linewidth=1)
for cap in boxes1['caps']:
    cap.set(color='#64000b', linewidth=1)
for flier in boxes1['fliers']:
    flier.set(color='#64000b')
for median in boxes1['medians']:
    median.set(color='#64000b', linewidth=2)


ax.set_xlim([1860,2101])
ax.set_xlabel('Years', fontweight='bold')
ax.yaxis.set_tick_params(left='off', right='off', labelright='on', labelleft='off', pad=7)
xmajorLocator = MultipleLocator(20)
xminorLocator = AutoMinorLocator(2)
ax.xaxis.set_major_locator(xmajorLocator)
ax.xaxis.set_minor_locator(xminorLocator)
ax.xaxis.set_tick_params(which='major',width=2)

ax2 = ax.twiny()
# ToE 1 Hist+rcp8.5 vs. HistNat (or vs. PiControl) boxes
boxes2 = ax2.boxplot(data1, vert=0, positions=ind+width, widths=width, whis=0)
for box in boxes2['boxes']:
    box.set(color='#d94c5b', linewidth=2)
for whisker in boxes2['whiskers']:
    whisker.set(color='#d94c5b', linestyle='-', linewidth=1)
for cap in boxes2['caps']:
    cap.set(color='#d94c5b', linewidth=1)
for flier in boxes2['fliers']:
    flier.set(color='#d94c5b')
for median in boxes2['medians']:
    median.set(color='#d94c5b', linewidth=2)


ax2.set_xlim([1860,2101])
ax2.set_yticks(ind)
ax2.set_yticklabels(labels, fontweight='bold')
ax2.yaxis.set_tick_params(left='off', right='off')
ax2.set_ylim([0,15])
xmajorLocator2 = MultipleLocator(20)
xminorLocator2 = AutoMinorLocator(2)
ax2.xaxis.set_major_locator(xmajorLocator2)
ax2.xaxis.set_minor_locator(xminorLocator2)

plt.setp(ax.get_yticklabels(), visible=True)
plt.setp(ax.get_xticklabels(), fontweight='bold')

ax2.axhline(y=ind[1], color='black', ls='--')
ax2.axhline(y=ind[3], color='black', ls='--')
ax2.axhline(y=ind[7], color='black', ls='--')
ax2.axhline(y=ind[10], color='black', ls='--')

# Domain labels
ax2.text(1860-17,ind[0], 'North \n Pac', ha='center', va='center', fontweight='bold', fontsize=13)
ax2.text(1860-17,ind[2], 'North \n Atl', ha='center', va='center', fontweight='bold', fontsize=13)
ax2.text(1860-17,ind[5], 'Southern \n Ocean', ha='center', va='center', fontweight='bold', fontsize=13)
ax2.text(1860-17,ind[8]+width, 'Northern \n ST', ha='center', fontweight='bold', fontsize=13)
ax2.text(1860-17,ind[12], 'Southern \n ST', ha='center', va='center', fontweight='bold', fontsize=13)

if use_piC == True:
    title = 'Hist + RCP8.5 vs. PiControl'
    end_name = 'use_piC'
    end_noise = 'RCP8.5 vs. PiControl'
else:
    title = 'Hist + RCP8.5 vs. HistNat'
    end_name = 'use_histNat'
    end_noise = 'RCP8.5 vs. HistNat'
plotTitle = 'ToE distribution for '+legVar+ ' in different regions \n '+title+ ' ('+str(nruns)+' runs)'
ax.set_title(plotTitle, y=1.08, fontweight='bold', va='center')
ax2.text(0.4,1.045, 'ToE [>1std]', color='#d94c5b',
         va='center', ha='center',transform=ax2.transAxes, fontweight='bold')
ax2.text(0.6,1.045, 'ToE [>2std]', color='#64000b',
         va='center', ha='center',transform=ax2.transAxes, fontweight='bold')

# Date
now = datetime.datetime.now()
date = now.strftime("%Y-%m-%d")

# Text at the bottom of the figure
plt.figtext(.8,.01,'Computed by : boxplot_ToE_rcp85.py, '+date, fontsize=8, ha='center')
plt.figtext(.2,.01,'Noise: %s' %(method_noise_rcp), fontsize=8, ha='center')
if use_piC :
    plt.figtext(.5,.01,'PiControl : mean(last_240_years)',fontsize=9,ha='center')

if runs_rcp == 'all':
    plotName = 'ToE_boxplot_RCP85_' + method_noise_rcp + '_' + end_name
else:
    plotName = 'ToE_boxplot_RCP85_' + method_noise_rcp + '_' + end_name + '_samerunsvsPiC'


if outfmt == 'view':
    plt.show()
else:
    plt.savefig('/home/ysilvy/Density_bining/Yona_analysis/figures/models/ToE/boxplots/'+plotName+'.png')

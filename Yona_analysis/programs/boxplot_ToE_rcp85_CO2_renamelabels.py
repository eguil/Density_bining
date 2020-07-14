#!/bin/env python
# -*- coding: utf-8 -*-

"""
Python matplotlib
Make whisker plots of ToE hist+RCP85 vs. histNat (or Picontrol) and 1%CO2 vs. PiControl in all 5 domains
"""

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset as open_ncfile
from maps_matplot_lib import defVarmme
from modelsDef import defModels, defModelsCO2piC
from matplotlib.ticker import AutoMinorLocator, MultipleLocator
import glob
import os
import datetime

# ----- Work -----

# Directory
indir_CO2piC = '/home/ysilvy/Density_bining/Yona_analysis/data/toe_1pctCO2vsPiC_average_signal/'
indir_rcphn = '/home/ysilvy/Density_bining/Yona_analysis/data/toe_rcp85_histNat_average_signal/'
indir_rcppiC = '/home/ysilvy/Density_bining/Yona_analysis/data/toe_rcp85_PiControl_average_signal/'

models = defModels()
modelspiC = defModelsCO2piC()

domains = ['Southern ST', 'SO', 'Northern ST', 'North Atlantic', 'North Pacific']

varname = defVarmme('salinity'); v = 'S'

method = 'average_signal' # Average signal and noise in the box, then compute ToE

# === INPUTS ===

# method_noise_rcphn = 'average_std' # Average the standard deviation of histNat in the specified domains
# method_noise_piC = 'average_std' # Average the standard deviation of PiC in the specified domains
method_noise_rcp = 'average_histNat' # Average histNat in the specified domains then determine the std of this averaged value
# (for hist+RCP8.5 vs. histNat or for hist+RCP8.5 vs. PiControl)
method_noise_piC = 'average_piC' # Average PiC in the specified domains then determine the std of this averaged value
# (for 1%CO2 vs. PiControl)

use_piC = False # signal = hist+RCP - average(histNat), noise = std(histNat)
# use_piC = True # signal = hist+RCP - PiControl, noise = std(PiControl)

if use_piC == True and method_noise_rcp == 'average_histNat':
    method_noise_rcp = 'average_piC' # Just change the notation for coherence and for data path

multstd = 2 # >1st or >2std

# runs_rcp = 'same' # Same runs (30 runs) for hist+RCP8.5 vs. histNat as for hist+RCP8.5 vs. PiControl
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
nMembers = np.ma.empty(len(models)) # Initialize array for keeping number of members per model

# -- Initialize varToE containing ToE of all runs
varToEA = np.ma.masked_all((nrunmax, len(domains)))
varToEP = np.ma.masked_all((nrunmax, len(domains)))
varToEI = np.ma.masked_all((nrunmax, len(domains)))

# -- Loop over models
if use_piC == True:
    indir = indir_rcppiC
else:
    indir = indir_rcphn
listfiles = glob.glob(indir + method_noise_rcp + '/*'+legVar+'*.nc')
nmodels = len(listfiles)

for i in range(nmodels):

    file_toe = listfiles[i]
    ftoe = open_ncfile(file_toe, 'r')
    name = os.path.basename(file_toe).split('.')[1]

    # If use same runs in vs. histNat as in vs. PiControl, take out deficient models
    if (runs_rcp == 'all') or (runs_rcp =='same' and name != 'GISS-E2-R' and name != 'FGOALS-g2' and name != 'MIROC-ESM'):

        # Read ToE (members, basin, domain)
        if multstd == 1:
            toeread = ftoe.variables[var + 'ToE1'][:]
        else:
            toeread = ftoe.variables[var + 'ToE2'][:]
        nMembers[i] = toeread.shape[0]
        print('- Reading ToE of %s with %d members'%(name,nMembers[i]))
        nruns1 = int(nruns + nMembers[i])
        
        toeread_median = np.ma.median(toeread,axis=0)

        # Save ToE
#         varToEA[nruns:nruns1,:] = toeread[:,1,:]
#         varToEP[nruns:nruns1,:] = toeread[:,2,:]
#         varToEI[nruns:nruns1,:] = toeread[:,3,:]
        
        # Inter-member median
        varToEA[i,:] = toeread_median[1,:] #toeread[:,1,:]
        varToEP[i,:] = toeread_median[2,:] #toeread[:,2,:]
        varToEI[i,:] = toeread_median[3,:] #toeread[:,3,:]

        nruns = nruns1


print('Total number of runs:', nruns)
# varToEA = varToEA[0:nruns,:]
# varToEP = varToEP[0:nruns,:]
# varToEI = varToEI[0:nruns,:]
varToEA = varToEA[0:i+1,:]
varToEP = varToEP[0:i+1,:]
varToEI = varToEI[0:i+1,:]

nruns = int(nruns)
nMembers = nMembers[np.ma.nonzero(nMembers)]


if runs_rcp == 'same':
        nmodels=nmodels-3

# == 1%CO2 vs. Pi Control ==

# -- Initialize varToE containing ToE
varToEA_CO2 = np.ma.masked_all((len(modelspiC),len(domains)))
varToEP_CO2 = np.ma.masked_all((len(modelspiC),len(domains)))
varToEI_CO2 = np.ma.masked_all((len(modelspiC),len(domains)))

for i, model in enumerate(modelspiC):

    print('- Reading ' + model['name'])

    # Read file
    file_CO2piC = 'cmip5.' + model['name'] + '.'+legVar+'_toe_1pctCO2vsPiControl_method2_' + method_noise_piC + '.nc'
    fpiC = open_ncfile(indir_CO2piC + method_noise_piC + '/' + file_CO2piC, 'r')

    # Read ToE (basin, domain)
    if multstd == 1:
        toeread = fpiC.variables[var + 'ToE1'][:]
    else:
        toeread = fpiC.variables[var + 'ToE2'][:]

    # Save ToE
    varToEA_CO2[i,:] = toeread[1,:]
    varToEP_CO2[i,:] = toeread[2,:]
    varToEI_CO2[i,:] = toeread[3,:]


# ----- Turn masked data into nans -----

#varToEA[np.ma.getmask(varToEA)] = np.nan
#varToEP[np.ma.getmask(varToEP)] = np.nan
#varToEI[np.ma.getmask(varToEI)] = np.nan
#varToEA_CO2[np.ma.getmask(varToEA_CO2)] = np.nan
#varToEP_CO2[np.ma.getmask(varToEP_CO2)] = np.nan
#varToEI_CO2[np.ma.getmask(varToEI_CO2)] = np.nan

# ----- Plot ------

maskdata  = np.nan

# New domain labels
new_domains = ['SO subpolar', 'SH subtropics', 'NH subtropics', 'subpolar North Pacific']
# regroup previous "North Atlantic" with NH subtropics, remove NH Indian


# ToE hist+rcp8.5 vs. histNat (or vs. PiControl)
data1 = [varToEA[:,1], varToEP[:,1], varToEI[:,1], maskdata, varToEA[:,0], varToEP[:,0], varToEI[:,0], maskdata, varToEA[:,3], varToEP[:,2], maskdata, varToEP[:,4]]


# ToE 1%CO2 vs. PiControl
data2 = [varToEA_CO2[:,1], varToEP_CO2[:,1], varToEI_CO2[:,1], maskdata, varToEA_CO2[:,0], varToEP_CO2[:,0], varToEI_CO2[:,0], maskdata, varToEA_CO2[:,3], varToEP_CO2[:,2], maskdata, varToEP_CO2[:,4]]


labels = ['Atlantic','Pacific','Indian','','Atlantic','Pacific','Indian','','Atlantic','Pacific','','']
N = 13
ind = np.arange(1,N)
width = 0.25

fig, ax = plt.subplots(figsize=(10,12))

ax.axvline(x=2005, color='black', ls=':')

red_crosses = dict(markeredgecolor='#c90016', marker='+',linewidth=1)
# ToE Hist+rcp8.5 vs. HistNat (or vs. PiControl) boxes
boxes1 = ax.boxplot(data1, vert=0, positions=ind-width, widths=width, whis=0,flierprops=red_crosses)
for box in boxes1['boxes']:
    box.set(color='#c90016', linewidth=2)
for whisker in boxes1['whiskers']:
    whisker.set(color='#c90016', linestyle='-', linewidth=1) ##663366
for cap in boxes1['caps']:
    cap.set(color='#c90016', linewidth=1)
#for flier in boxes1['fliers']:
#    flier.set(color='#c90016')
for median in boxes1['medians']:
    median.set(color='#c90016', linewidth=2) #666699


ax.set_xlim([1860,2110])
ax.set_xlabel('Years', fontweight='bold', fontsize=14)
plotTitle = 'ToE [>'+str(multstd)+'std] distribution for '+legVar+ ' in different regions'
#ax.set_title(plotTitle, y=1.08, fontweight='bold', va='bottom')
ax.yaxis.set_tick_params(left=False, right=False, labelright=False, labelleft=True,pad=-10)
ax.set_yticklabels(labels, horizontalalignment = 'left')
xmajorLocator = MultipleLocator(20)
xminorLocator = AutoMinorLocator(2)
ax.xaxis.set_major_locator(xmajorLocator)
ax.xaxis.set_minor_locator(xminorLocator)
ax.xaxis.set_tick_params(which='major',width=2)

ax2 = ax.twiny()
blue_crosses = dict(markeredgecolor='#004f82', marker='+',linewidth=1)
# ToE 1%CO2 vs. PiControl
boxes2 = ax2.boxplot(data2, vert=0, positions=ind+width, widths=width, whis=0,flierprops=blue_crosses)
for box in boxes2['boxes']:
    box.set(color='#004f82', linewidth=2)
for whisker in boxes2['whiskers']:
    whisker.set(color='#0072bb', linestyle='-', linewidth=1)
for cap in boxes2['caps']:
    cap.set(color='#004f82', linewidth=1)
#for flier in boxes2['fliers']:
#    flier.set(color='#004f82')
for median in boxes2['medians']:
    median.set(color='#004f82', linewidth=2) #a1caf1


ax2.set_xlim([0,250])
ax2.set_yticks(ind)
ax2.set_yticklabels(labels, fontweight='bold')
ax2.yaxis.set_tick_params(left=False, right=False)
ax2.set_ylim([0,N])
xmajorLocator2 = MultipleLocator(20)
xminorLocator2 = AutoMinorLocator(2)
ax2.xaxis.set_major_locator(xmajorLocator2)
ax2.xaxis.set_minor_locator(xminorLocator2)

plt.setp(ax2.get_yticklabels(), fontsize=12,fontweight='bold')
plt.setp(ax.get_yticklabels(), fontsize=12,fontweight='bold')
plt.setp(ax.get_xticklabels(), fontweight='bold',fontsize=14, rotation=20)
plt.setp(ax2.get_xticklabels(), fontweight='bold',fontsize=14, rotation=20)

ax2.axhline(y=ind[3], color='black', ls='--')
ax2.axhline(y=ind[7], color='black', ls='--')
ax2.axhline(y=ind[10], color='black', ls='--')

degree_sign= u'\N{DEGREE SIGN}'

# Domain labels
ax2.text(-43,ind[1], 'SO \n subpolar', ha='center', va='center', fontweight='bold', fontsize=13)
ax2.text(275,ind[1], '~40-60'+degree_sign+'S \n ~27-28kg.m-3', ha='center', va='center', fontsize=12)
ax2.text(-43,ind[5], 'SH \n subtropics', ha='center', va='center', fontweight='bold', fontsize=13)
ax2.text(275,ind[5], '~20-40'+degree_sign+'S \n ~25-26.5kg.m-3', ha='center', va='center', fontsize=12)
ax2.text(-43,ind[8]+width, 'NH \n subtropics', ha='center', va='center', fontweight='bold', fontsize=13)
ax2.text(275,ind[8], '~20-40'+degree_sign+'N \n ~26-27kg.m-3', ha='center', va='center', fontsize=12)
ax2.text(275,ind[9], '~20-40'+degree_sign+'N \n ~25-26kg.m-3', ha='center', va='center', fontsize=12)
ax2.text(-43,ind[11], 'Subpolar \n North Pacific', ha='center', fontweight='bold', fontsize=13)
ax2.text(275,ind[11], '~40-60'+degree_sign+'N \n ~26-27kg.m-3', ha='center', va='center', fontsize=12)

plt.subplots_adjust(left=0.18,right=0.85,top=0.91,bottom=0.1)

# Legend
# legendlabel = 'Hist + RCP8.5 vs. HistNat ('+str(nruns)+' runs) \n 1%CO2 vs. PiControl ('+str(len(modelspiC))+' runs)'
if use_piC == True:
    title = 'Hist + RCP8.5 vs. PiControl'
    end_name = 'use_piC'
    end_noise = 'RCP8.5 vs. PiControl'
else :
    title = 'Hist + RCP8.5 vs. HistNat'
    end_name = 'use_histNat'
    end_noise = 'RCP8.5 vs. HistNat'
ax2.text(0.5,1.045, title, color='#c90016',
         va='center', ha='center',transform=ax2.transAxes, fontweight='bold')
ax2.text(0.5,1.065, '1%CO2 vs. PiControl', color='#004f82',
         va='center', ha='center',transform=ax2.transAxes, fontweight='bold')

# Date
now = datetime.datetime.now()
date = now.strftime("%Y-%m-%d")

# Text at the bottom of the figure
#plt.figtext(.8,.01,'Computed by : boxplot_ToE_rcp85_CO2_renamelabels.py, '+date, fontsize=8, ha='center')
#plt.figtext(.2,.01,'Noise: %s %s %d std' %(method_noise_rcp, method_noise_piC, multstd), fontsize=8, ha='center')
if use_piC :
    plt.figtext(.5,.01,'PiControl : mean(last_240_years)',fontsize=9,ha='center')

if runs_rcp == 'all':
    plotName = 'ToE_boxplot_RCP85_1pctCO2_' + method_noise_rcp + '_' + method_noise_piC + '_' + end_name + '_' + str(multstd) + 'std_renamelabels_paper'
else:
    plotName = 'ToE_boxplot_RCP85_1pctCO2_' + method_noise_rcp + '_' + method_noise_piC + '_' + end_name + '_' + str(multstd) + 'std_samerunsvsPiC_renamelabels'


if outfmt == 'view':
    plt.show()
else:
    plt.savefig(plotName+'.png',dpi=150)

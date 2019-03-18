#!/bin/env python
# -*- coding: utf-8 -*-

"""
Python matplotlib
Make whisker plots of ToE 1%CO2 vs. PiControl in all 5 domains
>1std and >2std
"""

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset as open_ncfile
from maps_matplot_lib import defVarmme
from modelsDef import defModelsCO2piC
from matplotlib.ticker import AutoMinorLocator, MultipleLocator
import glob
import os
import datetime

# ----- Work -----

# Directory
indir_CO2piC = '/home/ysilvy/Density_bining/Yona_analysis/data/toe_1pctCO2vsPiC_average_signal/'

modelspiC = defModelsCO2piC()

domains = ['Southern ST', 'SO', 'Northern ST', 'North Atlantic', 'North Pacific']

varname = defVarmme('salinity'); v = 'S'

method = 'average_signal' # Average signal and noise in the box, then compute ToE

# === INPUTS ===

# method_noise_piC = 'average_std' # Average the standard deviation of PiC in the specified domains
method_noise_piC = 'average_piC' # Average PiC in the specified domains then determine the std of this averaged value

# output format
# outfmt = 'view'
outfmt = 'save'

# ===========

# ----- Variables ------
var = varname['var_zonal_w/bowl']
legVar = varname['legVar']
unit = varname['unit']


# ----- Read ToE for each model ------

# == 1%CO2 vs. Pi Control ==

# -- Initialize varToE containing ToE
varToEA_1 = np.ma.masked_all((len(modelspiC),len(domains)))
varToEP_1 = np.ma.masked_all((len(modelspiC),len(domains)))
varToEI_1 = np.ma.masked_all((len(modelspiC),len(domains)))
varToEA_2 = np.ma.masked_all((len(modelspiC),len(domains)))
varToEP_2 = np.ma.masked_all((len(modelspiC),len(domains)))
varToEI_2 = np.ma.masked_all((len(modelspiC),len(domains)))

for i, model in enumerate(modelspiC):

    print('- Reading ' + model['name'])

    # Read file
    file_CO2piC = 'cmip5.' + model['name'] + '.toe_1pctCO2vsPiControl_method2_' + method_noise_piC + '.nc'
    fpiC = open_ncfile(indir_CO2piC + method_noise_piC + '/' + file_CO2piC, 'r')

    # Read ToE (basin, domain)
    toe1read = fpiC.variables[var + 'ToE1'][:]
    toe2read = fpiC.variables[var + 'ToE2'][:]

    # Save ToE
    varToEA_1[i,:] = toe1read[1,:]
    varToEP_1[i,:] = toe1read[2,:]
    varToEI_1[i,:] = toe1read[3,:]
    varToEA_2[i,:] = toe2read[1,:]
    varToEP_2[i,:] = toe2read[2,:]
    varToEI_2[i,:] = toe2read[3,:]

# ----- Turn masked data into nans -----

varToEA_1[np.ma.getmask(varToEA_1)] = np.nan
varToEP_1[np.ma.getmask(varToEP_1)] = np.nan
varToEI_1[np.ma.getmask(varToEI_1)] = np.nan
varToEA_2[np.ma.getmask(varToEA_2)] = np.nan
varToEP_2[np.ma.getmask(varToEP_2)] = np.nan
varToEI_2[np.ma.getmask(varToEI_2)] = np.nan

# ----- Plot ------

maskdata  = np.nan

# ToE 1 [>1std] 1%CO2 vs. PiControl
data1 = [varToEA_1[:,0], varToEP_1[:,0], varToEI_1[:,0], maskdata, varToEP_1[:,2], varToEI_1[:,2], maskdata,
          varToEA_1[:,1], varToEP_1[:,1], varToEI_1[:,1], maskdata, varToEA_1[:,3], maskdata, varToEP_1[:,4]]
data1 = data1[::-1]

# ToE 2 [>2std] 1%CO2 vs. PiControl
data2 = [varToEA_2[:,0], varToEP_2[:,0], varToEI_2[:,0], maskdata, varToEP_2[:,2], varToEI_2[:,2], maskdata,
          varToEA_2[:,1], varToEP_2[:,1], varToEI_2[:,1], maskdata, varToEA_2[:,3], maskdata, varToEP_2[:,4]]
data2 = data2[::-1]


labels = ['','','','','Indian','Pacific','Atlantic','','Indian','Pacific','','Indian','Pacific','Atlantic']
N = 15
ind = np.arange(1,N)
width = 0.25

fig, ax = plt.subplots(figsize=(7,12))

# ToE 2 1%CO2 vs. PiControl boxes
boxes1 = ax.boxplot(data2, vert=0, positions=ind-width, widths=width, whis=0)
for box in boxes1['boxes']:
    box.set(color='#004f82', linewidth=2) # #0072bb
for whisker in boxes1['whiskers']:
    whisker.set(color='#004f82', linestyle='-', linewidth=1)
for cap in boxes1['caps']:
    cap.set(color='#004f82', linewidth=1)
for flier in boxes1['fliers']:
    flier.set(color='#004f82')
for median in boxes1['medians']:
    median.set(color='#004f82', linewidth=2)


ax.set_xlim([0,141])
ax.set_xlabel('Years', fontweight='bold')
ax.yaxis.set_tick_params(left='off', right='off', labelright='on', labelleft='off', pad=7)
xmajorLocator = MultipleLocator(20)
xminorLocator = AutoMinorLocator(2)
ax.xaxis.set_major_locator(xmajorLocator)
ax.xaxis.set_minor_locator(xminorLocator)
ax.xaxis.set_tick_params(which='major',width=2)

ax2 = ax.twiny()
# ToE 1 1%CO2 vs. PiControl boxes
boxes2 = ax2.boxplot(data1, vert=0, positions=ind+width, widths=width, whis=0)
for box in boxes2['boxes']:
    box.set(color='#4c9ccf', linewidth=2)
for whisker in boxes2['whiskers']:
    whisker.set(color='#4c9ccf', linestyle='-', linewidth=1)
for cap in boxes2['caps']:
    cap.set(color='#4c9ccf', linewidth=1)
for flier in boxes2['fliers']:
    flier.set(color='#4c9ccf')
for median in boxes2['medians']:
    median.set(color='#4c9ccf', linewidth=2)


ax2.set_xlim([0,141])
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
ax2.text(-12,ind[0], 'North \n Pac', ha='center', va='center', fontweight='bold', fontsize=12)
ax2.text(-12,ind[2], 'North \n Atl', ha='center', va='center', fontweight='bold', fontsize=12)
ax2.text(-12,ind[5], 'Southern \n Ocean', ha='center', va='center', fontweight='bold', fontsize=12)
ax2.text(-12,ind[8]+width, 'Northern \n ST', ha='center', fontweight='bold', fontsize=12)
ax2.text(-12,ind[12], 'Southern \n ST', ha='center', va='center', fontweight='bold', fontsize=12)


plotTitle = 'ToE distribution for '+legVar+ ' in different regions \n 1%CO2 vs. PiControl ('+str(len(modelspiC))+' runs)'
ax.set_title(plotTitle, y=1.08, fontweight='bold', va='center')
ax2.text(0.35,1.045, 'ToE [>1std]', color='#4c9ccf',
         va='center', ha='center',transform=ax2.transAxes, fontweight='bold')
ax2.text(0.65,1.045, 'ToE [>2std]', color='#004f82',
         va='center', ha='center',transform=ax2.transAxes, fontweight='bold')

# Date
now = datetime.datetime.now()
date = now.strftime("%Y-%m-%d")

# Text at the bottom of the figure
plt.figtext(.8,.01,'Computed by : boxplot_ToE_CO2.py, '+date, fontsize=8, ha='center')
plt.figtext(.2,.01,'Noise: %s' %(method_noise_piC), fontsize=8, ha='center')

plotName = 'ToE_boxplot_1pctCO2_' + method_noise_piC


if outfmt == 'view':
    plt.show()
else:
    plt.savefig('/home/ysilvy/Density_bining/Yona_analysis/figures/models/ToE/boxplots/'+plotName+'.png')

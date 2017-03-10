#!/bin/env python
# -*- coding: utf-8 -*-

"""
Python matplotlib
Make whisker plots of ToE (method 2) hist vs. histNat and 1%CO2 vs. PiControl in all 5 domains
"""

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset as open_ncfile
from maps_matplot_lib import defVarmme
from modelsDef import defModels, defModelsCO2piC
from matplotlib.ticker import AutoMinorLocator


# ----- Work -----

# Directory
indir_hhn = '/home/ysilvy/Density_bining/Yona_analysis/data/toe_histNat_average_signal/'
indir_CO2piC = '/home/ysilvy/Density_bining/Yona_analysis/data/toe_1pctCO2vsPiC_average_signal/'


models = defModels()
modelspiC = defModelsCO2piC()

domains = ['Southern ST', 'SO', 'Northern ST', 'North Atlantic', 'North Pacific']

varname = defVarmme('salinity'); v = 'S'

method = 'average_signal'

# ----- Variables ------
var = varname['var_zonal_w/bowl']
legVar = varname['legVar']
unit = varname['unit']


# ----- Read ToE for each model ------

# == Historical vs. historicalNat ==

nruns = 0 # Initialize total number of runs
nrunmax = 100
nMembers = np.ma.empty(len(models)) # Initialize array for keeping number of members per model

# -- Initialize varToE containing ToE of all runs
varToEA = np.ma.masked_all((nrunmax, len(domains)))
varToEP = np.ma.masked_all((nrunmax, len(domains)))
varToEI = np.ma.masked_all((nrunmax, len(domains)))

# -- Loop over models
for i, model in enumerate(models):

    file_toe = 'cmip5.' + model['name'] + '.toe_histNat_method2.nc'
    ftoe = open_ncfile(indir_hhn + file_toe, 'r')

    # Read ToE (members, basin, domain)
    toe2read = ftoe.variables[var + 'ToE2'][:]
    nMembers[i] = toe2read.shape[0]
    print '- Reading ToE of',model['name'], 'with', nMembers[i], 'members'
    nruns1 = nruns + nMembers[i]

    # Save ToE
    varToEA[nruns:nruns1,:] = toe2read[:,1,:]
    varToEP[nruns:nruns1,:] = toe2read[:,2,:]
    varToEI[nruns:nruns1,:] = toe2read[:,3,:]

    nruns = nruns1


print 'Total number of runs:', nruns
varToEA = varToEA[0:nruns,:]
varToEP = varToEP[0:nruns,:]
varToEI = varToEI[0:nruns,:]

nruns = int(nruns)


# == 1%CO2 vs. Pi Control ==

# -- Initialize varToE containing ToE
varToEA_CO2 = np.ma.masked_all((len(modelspiC),len(domains)))
varToEP_CO2 = np.ma.masked_all((len(modelspiC),len(domains)))
varToEI_CO2 = np.ma.masked_all((len(modelspiC),len(domains)))

for i, model in enumerate(modelspiC):

    print '- Reading', model['name']

    # Read file
    file_CO2piC = 'cmip5.' + model['name'] + '.toe_1pctCO2vsPiControl_method2.nc'
    fpiC = open_ncfile(indir_CO2piC + file_CO2piC, 'r')

    # Read ToE (basin, domain)
    toe2read = fpiC.variables[var + 'ToE2'][:]

    # Save ToE
    varToEA_CO2[i,:] = toe2read[1,:]
    varToEP_CO2[i,:] = toe2read[2,:]
    varToEI_CO2[i,:] = toe2read[3,:]


# ----- Turn masked data into nans -----

varToEA[np.ma.getmask(varToEA)] = np.nan
varToEP[np.ma.getmask(varToEP)] = np.nan
varToEI[np.ma.getmask(varToEI)] = np.nan
varToEA_CO2[np.ma.getmask(varToEA_CO2)] = np.nan
varToEP_CO2[np.ma.getmask(varToEP_CO2)] = np.nan
varToEI_CO2[np.ma.getmask(varToEI_CO2)] = np.nan

# ----- Plot ------

maskdata  = np.nan

# ToE hist vs. histNat
data1 = [varToEA[:,0], varToEP[:,0], varToEI[:,0], maskdata, varToEP[:,2], varToEI[:,2], maskdata,
          varToEA[:,1], varToEP[:,1], varToEI[:,1], maskdata, varToEA[:,3], maskdata, varToEP[:,4]]
data1 = data1[::-1]
# ToE 1%CO2 vs. PiControl
data2 = [varToEA_CO2[:,0], varToEP_CO2[:,0], varToEI_CO2[:,0], maskdata, varToEP_CO2[:,2], varToEI_CO2[:,2], maskdata,
          varToEA_CO2[:,1], varToEP_CO2[:,1], varToEI_CO2[:,1], maskdata, varToEA_CO2[:,3], maskdata, varToEP_CO2[:,4]]
data2 = data2[::-1]

labels = ['','','','','Indian','Pacific','Atlantic','','Indian','Pacific','','Indian','Pacific','Atlantic']
N = 15
ind = np.arange(1,N)
width = 0.25

fig, ax = plt.subplots(figsize=(10,12))

# ToE Hist vs. HistNat boxes
boxes1 = ax.boxplot(data1, vert=0, positions=ind-width, widths=width)
for box in boxes1['boxes']:
    box.set(color='#c90016', linewidth=2)
for whisker in boxes1['whiskers']:
    whisker.set(color='#c90016', linestyle='-', linewidth=1)
for cap in boxes1['caps']:
    cap.set(color='#c90016', linewidth=1)
for flier in boxes1['fliers']:
    flier.set(marker=None)
for median in boxes1['medians']:
    median.set(color='#ff2052', linewidth=2)


ax.set_xlim([1860,2010])
ax.set_xlabel('Years', fontweight='bold')
plotTitle = 'ToE distribution for '+legVar+ ' in different regions'
ax.set_title(plotTitle, y=1.08, fontweight='bold', va='bottom')
ax.yaxis.set_tick_params(left='off', right='off', labelright='on', labelleft='off', pad=7)
xminorLocator = AutoMinorLocator(2)
ax.xaxis.set_minor_locator(xminorLocator)

ax2 = ax.twiny()
# ToE 1%CO2 vs. PiControl
boxes2 = ax2.boxplot(data2, vert=0, positions=ind+width, widths=width)
for box in boxes2['boxes']:
    box.set(color='#0072bb', linewidth=2)
for whisker in boxes2['whiskers']:
    whisker.set(color='#0072bb', linestyle='-', linewidth=1)
for cap in boxes2['caps']:
    cap.set(color='#0072bb', linewidth=1)
for flier in boxes2['fliers']:
    flier.set(marker=None)
for median in boxes2['medians']:
    median.set(color='#a1caf1', linewidth=2)


ax2.set_xlim([0,150])
ax2.set_yticks(ind)
ax2.set_yticklabels(labels, fontweight='bold')
ax2.yaxis.set_tick_params(left='off', right='off')
ax2.set_ylim([0,15])
xminorLocator2 = AutoMinorLocator(2)
ax2.xaxis.set_minor_locator(xminorLocator2)

plt.setp(ax.get_yticklabels(), visible=True)

ax2.axhline(y=ind[1], color='black', ls='--')
ax2.axhline(y=ind[3], color='black', ls='--')
ax2.axhline(y=ind[7], color='black', ls='--')
ax2.axhline(y=ind[10], color='black', ls='--')

# Domain labels
ax2.text(-12,ind[0], 'North \n Pac', ha='center', va='center', fontweight='bold', fontsize=13)
ax2.text(-12,ind[2], 'North \n Atl', ha='center', va='center', fontweight='bold', fontsize=13)
ax2.text(-12,ind[5], 'Southern \n Ocean', ha='center', va='center', fontweight='bold', fontsize=13)
ax2.text(-12,ind[8]+width, 'Northern \n ST', ha='center', fontweight='bold', fontsize=13)
ax2.text(-12,ind[12], 'Southern \n ST', ha='center', va='center', fontweight='bold', fontsize=13)

# Legend
legendlabel = 'Hist vs. HistNat ('+str(nruns)+' runs) \n 1%CO2 vs. PiControl ('+str(len(modelspiC))+' runs)'
ax2.text(0.5,1.045, 'Hist vs. HistNat ('+str(nruns)+' runs)', color='#c90016',
         va='center', ha='center',transform=ax2.transAxes, fontweight='bold')
ax2.text(0.5,1.065, '1%CO2 vs. PiControl ('+str(len(modelspiC))+' runs)', color='#0072bb',
         va='center', ha='center',transform=ax2.transAxes, fontweight='bold')


plt.figtext(.8,.01,'Computed by : boxplot_ToE.py', fontsize=8, ha='center')
plt.figtext(.2,.01,'Method : %s' %(method,), fontsize=8, ha='center')

plotName = 'ToE_boxplot'

plt.show()
#plt.savefig('/home/ysilvy/Density_bining/Yona_analysis/figures/models/ToE/'+plotName+'.png')

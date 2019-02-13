#!/bin/env python
# -*- coding: utf-8 -*-

"""
Python matplotlib
Read noise estimate (computed in compute_noise_estimate.py) of historical, historicalNat and PiControl from file,
and make bar plots showing each basin and domain
"""

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset as open_ncfile
from maps_matplot_lib import defVarmme
from modelsDef import defModels, defModelsCO2piC


# ----- Work -----

# Directory
indir_noise = '/home/ysilvy/Density_bining/Yona_analysis/data/noise_estimate/'

models = defModels()
modelspiC = defModelsCO2piC()

domains = ['Southern ST', 'SO', 'Northern ST', 'North Atlantic', 'North Pacific']
domains1 = ['Southern ST', 'Northern ST'] # First bar chart
domains2 = ['SO', 'North Atlantic', 'North Pacific'] # Second bar chart

varname = defVarmme('salinity'); v = 'S'

# ----- Variables ------
var = varname['var_zonal_w/bowl']
legVar = varname['legVar']
unit = varname['unit']


# ----- Read noise for each model ------

# == Historical and historicalNat ==

nruns_h = 0 # Initialize total number of hist runs
nruns_hn = 0 # Initialize total number of histNat runs
nrunmax = 100
nMembers_h = np.ma.empty(len(models)) # Initialize array for keeping number of members per model
nMembers_hn = np.ma.empty(len(models))

# -- Initialize varnoise hist and histNat containing std of all runs for each basin
varnoise_ha = np.ma.masked_all((nrunmax, len(domains)))
varnoise_hp = np.ma.masked_all((nrunmax, len(domains)))
varnoise_hi = np.ma.masked_all((nrunmax, len(domains)))
varnoise_hna = np.ma.masked_all((nrunmax, len(domains)))
varnoise_hnp = np.ma.masked_all((nrunmax, len(domains)))
varnoise_hni = np.ma.masked_all((nrunmax, len(domains)))

# -- Loop over models

for i, model in enumerate(models):

    # Read file
    file = 'cmip5.' + model['name'] + '.noise_domains_hist_histNat.nc'
    f = open_ncfile(indir_noise + file, 'r')

    # Read noise (members, basin, domain)
    varstdh = f.variables[var+'stdh'][:]
    varstdhn = f.variables[var+'stdhn'][:]
    nMembers_h[i] = varstdh.shape[0]
    nMembers_hn[i] = varstdhn.shape[0]
    nruns1_h = nruns_h + nMembers_h[i]
    nruns1_hn = nruns_hn + nMembers_hn[i]
    print '- Reading', model['name'], 'with', nMembers_h[i], 'hist members and', nMembers_hn[i], 'histNat members'
    print ''

    # Save noise
    varnoise_ha[nruns_h:nruns1_h,:] = varstdh[:,1,:]
    varnoise_hp[nruns_h:nruns1_h,:] = varstdh[:,2,:]
    varnoise_hi[nruns_h:nruns1_h,:] = varstdh[:,3,:]
    varnoise_hna[nruns_hn:nruns1_hn,:] = varstdhn[:,1,:]
    varnoise_hnp[nruns_hn:nruns1_hn,:] = varstdhn[:,2,:]
    varnoise_hni[nruns_hn:nruns1_hn,:] = varstdhn[:,3,:]

    nruns_h = nruns1_h
    nruns_hn = nruns1_hn

varnoise_ha = varnoise_ha[0:nruns_h,:]
varnoise_hp = varnoise_hp[0:nruns_h,:]
varnoise_hi = varnoise_hi[0:nruns_h,:]
varnoise_hna = varnoise_hna[0:nruns_hn,:]
varnoise_hnp = varnoise_hnp[0:nruns_hn,:]
varnoise_hni = varnoise_hni[0:nruns_hn,:]


# == Pi Control ==

# -- Initialize varnoise Pi Control containing std for each basin
varnoise_piCa = np.ma.masked_all((len(modelspiC),len(domains)))
varnoise_piCp = np.ma.masked_all((len(modelspiC),len(domains)))
varnoise_piCi = np.ma.masked_all((len(modelspiC),len(domains)))

# -- Loop over models

for i, model in enumerate(modelspiC):

    print '- Reading', model['name']

    # Read file
    file_piC = 'cmip5.' + model['name'] + '.noise_domains_PiControl.nc'
    fpiC = open_ncfile(indir_noise + file_piC, 'r')

    # Read noise (basin,domain)
    varstdpiC = fpiC.variables[var+'stdpiC'][:]

    varstdpiC[varstdpiC>10] = np.ma.masked # Mask potential errors

    # Save noise
    varnoise_piCa[i,:] = varstdpiC[1,:]
    varnoise_piCp[i,:] = varstdpiC[2,:]
    varnoise_piCi[i,:] = varstdpiC[3,:]

# -- Compute inter-model std of noise
stdnoise_ha = np.ma.std(varnoise_ha, axis=0)
stdnoise_hp = np.ma.std(varnoise_hp, axis=0)
stdnoise_hi = np.ma.std(varnoise_hi, axis=0)
stdnoise_hna = np.ma.std(varnoise_hna, axis=0)
stdnoise_hnp = np.ma.std(varnoise_hnp, axis=0)
stdnoise_hni = np.ma.std(varnoise_hni, axis=0)
stdnoise_piCa = np.ma.std(varnoise_piCa, axis=0)
stdnoise_piCp = np.ma.std(varnoise_piCp, axis=0)
stdnoise_piCi = np.ma.std(varnoise_piCi, axis=0)

# -- Compute average of all runs/models of noise
meannoise_ha = np.ma.average(varnoise_ha, axis=0)
meannoise_hp = np.ma.average(varnoise_hp, axis=0)
meannoise_hi = np.ma.average(varnoise_hi, axis=0)
meannoise_hna = np.ma.average(varnoise_hna, axis=0)
meannoise_hnp = np.ma.average(varnoise_hnp, axis=0)
meannoise_hni = np.ma.average(varnoise_hni, axis=0)
meannoise_piCa = np.ma.average(varnoise_piCa, axis=0)
meannoise_piCp = np.ma.average(varnoise_piCp, axis=0)
meannoise_piCi = np.ma.average(varnoise_piCi, axis=0)

# For displaying purposes, '.0' remove decimal
nruns_h = int(nruns_h)
nruns_hn = int(nruns_hn)

# ----- Plot ------

# Legend
labels = ['PiControl ('+str(len(modelspiC))+' runs)', 'HistNat ('+str(nruns_hn)+' runs)', 'Hist ('+str(nruns_h)+' runs)']

# == Graph 1 : Southern ST and Northern ST ==

maskdata  = np.ma.masked

piC = [meannoise_piCa[0], meannoise_piCp[0], meannoise_piCi[0], maskdata, meannoise_piCp[2], meannoise_piCi[2]]
piCstd = [stdnoise_piCa[0], stdnoise_piCp[0], stdnoise_piCi[0], maskdata, stdnoise_piCp[2], stdnoise_piCi[2]]
histNat = [meannoise_hna[0], meannoise_hnp[0], meannoise_hni[0], maskdata, meannoise_hnp[2], meannoise_hni[2]]
histNatstd = [stdnoise_hna[0], stdnoise_hnp[0], stdnoise_hni[0], maskdata, stdnoise_hnp[2], stdnoise_hni[2]]
hist = [meannoise_ha[0], meannoise_hp[0], meannoise_hi[0], maskdata, meannoise_hp[2], meannoise_hi[2]]
histstd = [stdnoise_ha[0], stdnoise_hp[0], stdnoise_hi[0], maskdata, stdnoise_hp[2], stdnoise_hi[2]]

N = 6
ind = np.arange(N) # x locations for groups
width = 0.25

fig, ax = plt.subplots()

# PiControl
piCbars = ax.bar(ind-width, piC, width, color='#83b2d0', yerr=piCstd, ecolor='black')
# HistNat
histNatbars = ax.bar(ind, histNat, width, color='#10c390', yerr=histNatstd, ecolor='black')
# Hist
histbars = ax.bar(ind+width, hist, width, color='#f04900', yerr=histstd, ecolor='black')


ax.set_ylabel('Variability (%s)' %(unit,), fontweight='bold')
plotTitle = 'Variability in the Subtropics for ' + legVar
ax.set_title(plotTitle, fontweight='bold', va='bottom')
ax.set_xticks(ind+width/2)
ax.set_xticklabels(['Atlantic','Pacific','Indian','','Pacific','Indian'])
ax.xaxis.set_tick_params(bottom='off', top='off')

ax.text(ind[1]+width/2,-0.01, 'Southern ST', ha='center', fontweight='bold', fontsize='13')
ax.text(ind[-1]-3*width/2,-0.01, 'Northern ST', ha='center', fontweight='bold', fontsize='13')

ax.legend([piCbars[0],histNatbars[0],histbars[0]],labels, fontsize=12, loc='upper left')

plt.figtext(.5,.01,'Computed by : noise_estimate_barplot.py',fontsize=8,ha='center')

plotName = 'variability_salinity_barchart1'

plt.savefig('/home/ysilvy/Density_bining/Yona_analysis/figures/models/noise_estimate/'+plotName+'.png', bbox_inches='tight')



# == Graph 2 : Southern Ocean, North Atlantic and North Pacific ==

maskdata  = np.ma.masked

piC = [meannoise_piCa[1], meannoise_piCp[1], meannoise_piCi[1], maskdata, meannoise_piCa[3], maskdata, meannoise_piCp[4]]
piCstd = [stdnoise_piCa[1], stdnoise_piCp[1], stdnoise_piCi[1], maskdata, stdnoise_piCa[3], maskdata, stdnoise_piCp[4]]
histNat = [meannoise_hna[1], meannoise_hnp[1], meannoise_hni[1], maskdata, meannoise_hna[3], maskdata, meannoise_hnp[4]]
histNatstd = [stdnoise_hna[1], stdnoise_hnp[1], stdnoise_hni[1], maskdata, stdnoise_hna[3], maskdata, stdnoise_hnp[4]]
hist = [meannoise_ha[1], meannoise_hp[1], meannoise_hi[1], maskdata, meannoise_ha[3], maskdata, meannoise_hp[4]]
histstd = [stdnoise_ha[1], stdnoise_hp[1], stdnoise_hi[1], maskdata, stdnoise_ha[3], maskdata, stdnoise_hp[4]]

N = 7
ind = np.arange(N)
width = 0.25

fig, ax = plt.subplots()

# PiControl
piCbars = ax.bar(ind-width, piC, width, color='#83b2d0', yerr=piCstd, ecolor='black')
# HistNat
histNatbars = ax.bar(ind, histNat, width, color='#10c390', yerr=histNatstd, ecolor='black')
# Hist
histbars = ax.bar(ind+width, hist, width, color='#f04900', yerr=histstd, ecolor='black')

ax.set_ylabel('Variability (%s)' %(unit,), fontweight='bold')
plotTitle = 'Variability in polar/subpolar regions for ' + legVar
ax.set_title(plotTitle, fontweight='bold', va='bottom')
ax.set_xticks(ind+width/2)
ax.set_xticklabels(['Atlantic','Pacific','Indian','','','',''])
ax.xaxis.set_tick_params(bottom='off', top='off')

ax.text(ind[1]+width/2,-0.015, 'Southern Ocean', ha='center', fontweight='bold', fontsize='13')
ax.text(ind[4]+width/2,-0.015, 'Nort Atl', ha='center', fontweight='bold', fontsize='13')
ax.text(ind[6]+width/2,-0.015, 'Nort Pac', ha='center', fontweight='bold', fontsize='13')

ax.legend([piCbars[0],histNatbars[0],histbars[0]],labels, fontsize=12, loc='upper left')

plt.figtext(.5,.01,'Computed by : noise_estimate_barplot.py',fontsize=8,ha='center')

plotName = 'variability_salinity_barchart2'
plt.savefig('/home/ysilvy/Density_bining/Yona_analysis/figures/models/noise_estimate/'+plotName+'.png', bbox_inches='tight')

#plt.show()
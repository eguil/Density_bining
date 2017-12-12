#!/bin/env python
# -*- coding: utf-8 -*-

"""
Plot Southern Subtropics time series hist+rcp8.5 vs. histNat
One subplot per model, one curve per run

"""

import numpy as np
from netCDF4 import Dataset as open_ncfile
import matplotlib.pyplot as plt
from maps_matplot_lib import defVarmme, averageDom
from modelsDef import defModels
from libToE import ToEdomainhistvshistNat
from matplotlib.ticker import AutoMinorLocator, MultipleLocator
import glob


# ----- Workspace ------

indir_histrcp85 = '/home/ysilvy/Density_bining/Yona_analysis/data/hist_rcp85/'
indir_histNat = '/data/ericglod/Density_binning/Prod_density_april15/mme_histNat/'
indir_piC = '/data/ericglod/Density_binning/Prod_density_april15/mme_piControl/'
indir_noise = '/home/ysilvy/Density_bining/Yona_analysis/data/noise_estimate/'

models = defModels()

# ----- Work ------

varname = defVarmme('salinity'); v = 'S'

method = 'average_signal' # Average signal and noise in the box, then compute ToE

# -- Choose which 'noise' to use
# method_noise = 'average_std' # Average the standard deviation of histNat or PiC in the specified domains
method_noise = 'average_histNat' # Average histNat or PiC in the specified domains then determine the std
# of this averaged value

domains = ['Southern ST', 'SO', 'Northern ST', 'North Atlantic', 'North Pacific']
domain_name = 'Southern ST'
idomain = 0

# ibasin = 1 ; basin_name = 'Atlantic' # Atlantic southern ST
ibasin = 2 ; basin_name = 'Pacific' # Pacific southern ST

use_piC = False # Over projection period, signal = RCP-average(histNat), noise = std(histNat)
# use_piC = True # Over projection period, signal = RCP-average(PiControl), noise = std(PiControl)

# Choose random file to read only the basic variables and properties common to all files
file = 'cmip5.' + models[1]['name'] + '.historicalNat.ensm.an.ocn.Omon.density.ver-' + \
       models[1]['file_end_histNat'] + '_zon2D.nc'
f = open_ncfile(indir_histNat + file,'r')

lat = f.variables['latitude'][:]; latN = lat.size
density = f.variables['lev'][:]; levN = density.size
timN = 240
time_label = np.arange(1860,2100)
var = varname['var_zonal_w/bowl']
basinN = 4

# Define variable properties
legVar = varname['legVar']


# ----- Initialize plot ------

if use_piC == True:
    nmodels = 8
else :
    nmodels = 11
if nmodels % 2 == 0:
    nrows = nmodels/2
else:
    nrows = nmodels/2 +1
fig, axes = plt.subplots(nrows = nrows, ncols=2, sharex = True, figsize=(12,15))

# ----- Average signal and noise for each run ------
j = 0 # Count models used

for i, model in enumerate(models):

    if (use_piC != True) or (use_piC == True and model['name'] != 'GISS-E2-R' and model['name'] != 'FGOALS-g2'
                             and model['name'] != 'MIROC-ESM'):

        # Define ax
        if j < nrows :
            ax = axes[j,0]
        else :
            ax = axes[j-nrows,1]

        # Read hist+rcp8.5 files
        listruns = glob.glob(indir_histrcp85 + 'cmip5.' + model['name'] + '.' + '*zon2D.nc')
        nruns = len(listruns)
        if nruns != 0:
            print('Working on %s'%(model['name'],))

            # Index of common time interval
            tstart = model['props'][2]
            tend = model['props'][3]

            # Read histNat ensemble mean
            filehn = 'cmip5.' + model['name'] + '.historicalNat.ensm.an.ocn.Omon.density.ver-' \
                     + model['file_end_histNat'] + '_zon2D.nc'
            fhn = open_ncfile(indir_histNat + filehn,'r')
            # Read var histNat
            varhn = fhn.variables[var][:,ibasin,:,:].squeeze()
            # Read std of histNat (max std of all runs for each model)
            if method_noise == 'average_std':
                varstd = fhn.variables[var+'Std'][ibasin,:,:].squeeze()
            # Compute time average of the whole histNat series (signal over projection = RCP - mean(histNat))
            meanvarhn = np.ma.mean(varhn, axis=0)

            if use_piC == True:
                # Read and Compute time average of PiControl over last 95 years + std of PiControl
                filepiC = glob.glob(indir_piC + 'cmip5.' + model['name'] + '.' + '*zon2D.nc')[0]
                fpiC = open_ncfile(filepiC,'r')
                varpiC = fpiC.variables[var][-95:,ibasin,:,:].squeeze()
                meanvarpiC = np.ma.mean(varpiC, axis=0)
                stdvarpiC = np.ma.std(varpiC, axis=0)

            # Initialize varsignal containing averaged signal for each run
            varsignal = np.ma.masked_all((timN,nruns))

            if method_noise == 'average_histNat':
                # Here we read the std of the averaged histNat for all runs, then take the max as our noise
                filenoise = 'cmip5.' + model['name'] + '.noise_domains_hist_histNat.std_of_average.nc'
                fnoise = open_ncfile(indir_noise + filenoise,'r')
                varstd = fnoise.variables[var+'stdhn'][:,ibasin,idomain].squeeze()
                varnoise = np.ma.max(varstd)

            # Select domain to average
            domain = ToEdomainhistvshistNat(model['name'], domain_name)[0]

            if method_noise == 'average_std':
                # Average noise
                if domain[basin_name] != None:
                    varnoise = averageDom(varstd, 2, domain[basin_name], lat, density)

            # Average noise over projection period if we want to use PiControl
            if use_piC == True:
                if method_noise == 'average_std':
                    varnoise2 = averageDom(stdvarpiC, 2, domain[basin_name], lat, density)
                else:
                    varnoise2 = np.ma.std(averageDom(varpiC, 3, domain[basin_name], lat, density), axis=0)

            # Loop over number of runs
            for k in range(nruns):
                # Read file
                fhrcp = open_ncfile(listruns[k],'r')
                # Read var hist + rcp85
                varh = fhrcp.variables[var][tstart:tend,ibasin,:,:].squeeze()
                varrcp = fhrcp.variables[var][-95:,ibasin,:,:].squeeze()

                # Average signal hist - histNat over historical period,
                # rcp85 - mean(histNat) or rcp85 - mean(PiC) over projection period
                varsignal[0:145,k] = averageDom(varh-varhn, 3, domain[basin_name], lat, density)
                if use_piC != True:
                    varsignal[145:,k] = averageDom(varrcp-meanvarhn, 3, domain[basin_name], lat, density)
                else:
                    varsignal[145:,k] = averageDom(varrcp-meanvarpiC, 3, domain[basin_name], lat, density)

                # Plot run k
                ax.plot(time_label,varsignal[:,k])

            # Plot noise
            if use_piC == False:
                ax.axhline(y=2*varnoise,color='black',ls='--')
                ax.axhline(y=-2*varnoise,color='black',ls='--')
            else :
                ax.hlines(y=2*varnoise,xmin=1860,xmax=2005,colors='black',linestyles='--')
                ax.hlines(y=-2*varnoise,xmin=1860,xmax=2005,colors='black',linestyles='--')
                ax.hlines(y=2*varnoise2,xmin=2005,xmax=2100,colors='black',linestyles='--')
                ax.hlines(y=-2*varnoise2,xmin=2005,xmax=2100,colors='black',linestyles='--')
            ax.axhline(y=0,color='black',ls=':')
            ax.axvline(x=2005,color='black',ls=':')

            ax.set_xlim([1860,2100])
            ax.tick_params(axis='x',which='both',top='off')
            ax.set_xticks(np.arange(1860,2101,40))
            # xmajorLocator = MultipleLocator(40)
            xminorLocator = AutoMinorLocator(2)
            # ax.xaxis.set_major_locator(xmajorLocator)
            ax.xaxis.set_minor_locator(xminorLocator)

            if nruns>1:
                subplot_title = '%s (%d members)'%(model['name'],nruns)
            else :
                subplot_title = '%s (%d member)'%(model['name'],nruns)
            ax.set_title(subplot_title,va='center',fontsize=12,fontweight='bold')

            j = j+1

plt.subplots_adjust(left=0.08, right=0.97) #hspace=.0001, wspace=0.05,

plt.setp(ax.get_xticklabels(), visible=True)

# Add a big axes, hide frame, for common labels
ax1=fig.add_subplot(111, frameon=False)
# Hide tick and tick label of the big axes
plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
plt.ylabel('Salinity change (PSU)', fontweight='bold')
ax1.set_yticks([])
ax1.yaxis.labelpad = 40

# Delete last subplot in case of odd number
if nmodels % 2 != 0:
    fig.delaxes(axes[-1,-1])

if method_noise == 'average_histNat':
    method_noise = 'std_of_average'
if use_piC == True:
    title = 'Evolution of salinity change in the %s Southern Subtropics \n ' \
            'Hist vs. HistNat + RCP8.5 vs. PiControl'%(basin_name,)
    end_name = 'use_piC'
    end_noise = 'RCP8.5 vs. PiControl'
    end_title = 'RCP85vsPiC'
else :
    title = 'Evolution of salinity change in the %s Southern Subtropics \n' \
            'Hist + RCP8.5 vs. HistNat'%(basin_name,)
    end_name = 'use_histNat'
    end_noise = 'RCP8.5 vs. HistNat'
    end_title = 'RCP85vsHistNat'

plt.suptitle(title, fontweight='bold', fontsize=14, verticalalignment='top')

# Foot notes
plt.figtext(.8,.01,'Computed by : plot_southern_ST.py', fontsize=8, ha='center')
plt.figtext(.2,.01,'Method: %s  Noise: %s %s' %(method, method_noise, end_noise), fontsize=8, ha='center')

plotName = 'Salinitychange_SouthernST_'+basin_name+'_'+end_title

plt.show()
# plt.savefig('/home/ysilvy/Density_bining/Yona_analysis/figures/models/time_series/'+plotName+'.png')

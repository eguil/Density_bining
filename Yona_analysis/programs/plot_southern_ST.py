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
from libToE import ToEdomainrcp85vshistNat
from matplotlib.ticker import AutoMinorLocator, MultipleLocator
import glob
import datetime


# ----- Workspace ------

indir_histrcp85 = '/home/ysilvy/Density_bining/Yona_analysis/data/hist_rcp85/'
indir_histNat = '/data/ericglod/Density_binning/Prod_density_april15/mme_histNat/'
indir_piC = '/data/ericglod/Density_binning/Prod_density_april15/mme_piControl/'
indir_noise = '/home/ysilvy/Density_bining/Yona_analysis/data/noise_estimate/RCP85vshistNat_domains/'

models = defModels()

# ----- Work ------

#varname = defVarmme('salinity'); v = 'S'
varname = defVarmme('depth'); v = 'Z'

method = 'average_signal' # Average signal and noise in the box, then compute ToE

# -- Choose which 'noise' to use
# method_noise = 'average_std' # Average the standard deviation of histNat or PiC in the specified domains
method_noise = 'average_histNat' # Average histNat or PiC in the specified domains then determine the std
# of this averaged value

domains = ['Southern ST', 'SO', 'Northern ST', 'North Atlantic', 'North Pacific']
idomain = 0
domain_name = domains[idomain]
#signal_domain = 'fresher'

#ibasin = 1 ; basin_name = 'Atlantic' # Atlantic 
ibasin = 2 ; basin_name = 'Pacific' # Pacific 
#ibasin = 3 ; basin_name = 'Indian' # Indian

use_piC = False # Over projection period, signal = RCP-average(histNat), noise = std(histNat)
#use_piC = True # Over projection period, signal = RCP-average(PiControl), noise = std(PiControl)

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
unit = varname['unit']

# ----- Initialize plot ------

if use_piC == True:
    nmodels = 8
else :
    nmodels = 11
if nmodels % 2 == 0:
    nrows = nmodels/2
else:
    nrows = nmodels/2 +1
fig, axes = plt.subplots(nrows = nrows, ncols=2, sharex = True, sharey=True, figsize=(12,15))

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

            if use_piC == False:
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

            else:
                # Read and Compute time average of PiControl over last 240 years + std of PiControl
                filepiC = glob.glob(indir_piC + 'cmip5.' + model['name'] + '.' + '*zon2D.nc')[0]
                fpiC = open_ncfile(filepiC,'r')
                varpiC = fpiC.variables[var][-240:,ibasin,:,:].squeeze()
                meanvarpiC = np.ma.mean(varpiC, axis=0)
                stdvarpiC = np.ma.std(varpiC, axis=0)

            # Initialize varsignal containing averaged signal for each run
            varsignal = np.ma.masked_all((timN,nruns))
            
            # Select domain to average
            domain = ToEdomainrcp85vshistNat(model['name'], domain_name)[0]

            if use_piC == False:
                if method_noise == 'average_histNat':
                    # Here we read the std of the averaged histNat for all runs, then take the max as our noise
                    if v == 'S': # Salinity
                        filenoise = 'cmip5.' + model['name'] + '.noise_domains_hist_histNat.std_of_average.nc'
                    else: # Depth
                        filenoise = 'cmip5.' + model['name'] + '.Depth_noise_domains_hist_histNat.std_of_average.nc'
                    fnoise = open_ncfile(indir_noise + filenoise,'r')
                    varstd = fnoise.variables[var+'stdhn'][:,ibasin,idomain].squeeze()
                    varnoise = np.ma.max(varstd)

                if method_noise == 'average_std':
                    # Average noise
                    if domain[basin_name] != None:
                        varnoise = averageDom(varstd, 2, domain[basin_name], lat, density)

            # Average noise if we want to use PiControl
            else:
                if method_noise == 'average_std':
                    varnoise = averageDom(stdvarpiC, 2, domain[basin_name], lat, density)
                else:
                    varnoise = np.ma.std(averageDom(varpiC, 3, domain[basin_name], lat, density), axis=0)

            # Loop over number of runs
            for k in range(nruns):
                # Read file
                fhrcp = open_ncfile(listruns[k],'r')
                # Read var hist + rcp85
                varh = fhrcp.variables[var][tstart:tend,ibasin,:,:].squeeze()
                varrcp = fhrcp.variables[var][tend:tend+95,ibasin,:,:].squeeze()

                # Average signal hist - histNat or hist - mean(PiC) over historical period
                # rcp85 - mean(histNat) or rcp85 - mean(PiC) over projection period
                if use_piC != True:
                    varsignal[0:145,k] = averageDom(varh-varhn, 3, domain[basin_name], lat, density)
                    varsignal[145:,k] = averageDom(varrcp-meanvarhn, 3, domain[basin_name], lat, density)
                else:
                    varsignal[0:145,k] = averageDom(varh-meanvarpiC, 3, domain[basin_name], lat, density)
                    varsignal[145:,k] = averageDom(varrcp-meanvarpiC, 3, domain[basin_name], lat, density)

                # # Don't plot runs where the signal is of opposite sign than expected
                # if np.ma.mean(varsignal[-5:,k],axis=0) <= 2*varnoise:

                # Plot signal run k
                ax.plot(time_label,varsignal[:,k],zorder=3)

            # Plot noise
            #ax.axhline(y=2*varnoise,color='black',ls='--')
            #ax.axhline(y=-2*varnoise,color='black',ls='--')
            ax.fill_between(time_label,-2*varnoise,-varnoise,facecolor='0.3',edgecolor=None,alpha=0.3,zorder=1)
            ax.fill_between(time_label,2*varnoise,varnoise,facecolor='0.3',edgecolor=None,alpha=0.3,zorder=1)
            #ax.axhline(y=varnoise,color='grey',ls='--')
            #ax.axhline(y=-varnoise,color='grey',ls='--')
            ax.fill_between(time_label,-varnoise,varnoise,facecolor='0.7',edgecolor=None,alpha=0.3,zorder=2)
           
            ax.axhline(y=0,color='black',ls=':')
            ax.axvline(x=2005,color='black',ls=':')

            ax.set_xlim([1860,2100])
            ax.tick_params(axis='x',which='both',top='off')
            ax.set_xticks(np.arange(1860,2101,40))
            # xmajorLocator = MultipleLocator(40)
            xminorLocator = AutoMinorLocator(2)
            # ax.xaxis.set_major_locator(xmajorLocator)
            ax.xaxis.set_minor_locator(xminorLocator)
            
            for tk in ax.get_yticklabels():
                tk.set_visible(True)
            for tk in ax.get_xticklabels():
                tk.set_visible(True)

            if nruns>1:
                subplot_title = '%s (%d members)'%(model['name'],nruns)
            else :
                subplot_title = '%s (%d member)'%(model['name'],nruns)
            ax.set_title(subplot_title,va='center',fontsize=12,fontweight='bold')

            j = j+1

# if domain_name == 'Southern ST' or domain_name == 'Northern ST':
#     ax.set_ylim([-0.5,0.2])

plt.subplots_adjust(left=0.08, right=0.97, bottom=0.05, top=0.93) #hspace=.0001, wspace=0.05,

plt.setp(ax.get_xticklabels(), visible=True)

# Add a big axes, hide frame, for common labels
ax1=fig.add_subplot(111, frameon=False)
# Hide tick and tick label of the big axes
plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
plt.ylabel(legVar+' change ('+unit+')', fontweight='bold')
ax1.set_yticks([])
ax1.yaxis.labelpad = 40

# Delete last subplot in case of odd number
if nmodels % 2 != 0:
    fig.delaxes(axes[-1,-1])

if method_noise == 'average_histNat':
    method_noise = 'std_of_average'
if use_piC == True:
    title = 'Evolution of '+legVar+' change signal in the %s Southern Subtropics \n ' \
            'Hist + RCP8.5 vs. PiControl'%(basin_name,)
    end_name = 'use_piC'
    end_noise = 'RCP8.5 vs. PiControl'
    end_title = 'RCP85vsPiC'
else :
    title = 'Evolution of '+legVar+' change signal in the %s Southern Subtropics \n' \
            'Hist + RCP8.5 vs. HistNat'%(basin_name,)
    end_name = 'use_histNat'
    end_noise = 'RCP8.5 vs. HistNat'
    end_title = 'RCP85vsHistNat'

plt.suptitle(title, fontweight='bold', fontsize=14, verticalalignment='top')

# Date
now = datetime.datetime.now()
date = now.strftime("%Y-%m-%d")

# Foot notes
plt.figtext(.8,.01,'Computed by : plot_southern_ST.py  '+date, fontsize=8, ha='center')
plt.figtext(.2,.01,'Method: %s  Noise: %s %s' %(method, method_noise, end_noise), fontsize=8, ha='center')

plotName = legVar+'change_SouthernST_'+basin_name+'_'+end_title

#plt.show()
plt.savefig('/home/ysilvy/figures/models/time_series/'+plotName+'.pdf')

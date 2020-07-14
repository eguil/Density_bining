#!/bin/env python
# -*- coding: utf-8 -*-

"""
Python matplotlib
Make whisker plots of ToE hist+RCP85 vs. histNat (or Picontrol) in all 5 domains, with x axis being the GSAT anomaly
>1std and >2std
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
from functions_ToE import read_toe_rcp85, read_gsat_rcp85, read_toe_1pctCO2, read_gsat_1pctCO2


# ----- Work -----

# Directory
indir_rcphn = '/home/ysilvy/Density_bining/Yona_analysis/data/toe_rcp85_histNat_average_signal/'
indir_rcppiC = '/home/ysilvy/Density_bining/Yona_analysis/data/toe_rcp85_PiControl_average_signal/'
indir_CO2piC = '/home/ysilvy/Density_bining/Yona_analysis/data/toe_1pctCO2vsPiC_average_signal/'
indir_gsat = '/home/ysilvy/Density_bining/Yona_analysis/data/gsat/'

models = defModels()
modelsCO2 = defModelsCO2piC()

domains = ['Southern ST', 'SO', 'Northern ST', 'North Atlantic', 'North Pacific']

varname = defVarmme('salinity'); v = 'S'

method = 'average_signal' # Average signal and noise in the box, then compute ToE
method_noise_histNat = 'average_histNat' # Average histNat (or PiControl) in the specified domains then determine the std of this averaged value
method_noise_piC = 'average_piC'

# === INPUTS ===

# -- Choose 2 datasets to plot/compare on the figure
# work = 'rcp85_histNat_1_2std'
work = 'rcp85_histNat_PiControl_2std'
# work = 'rcp85_histNat_1pctCO2_2std'

# output format
#outfmt = 'view'
outfmt = 'save'

# ===========

# ----- Variables ------
var = varname['var_zonal_w/bowl']
legVar = varname['legVar']
unit = varname['unit']

timN = 240

degree_sign= u'\N{DEGREE SIGN}'

# ------ Define directories, plot names, etc.. according to work -----

if work == 'rcp85_histNat_1_2std':
    indir1 = indir_rcphn
    indir2 = indir_rcphn
    varread1 = var+'ToE2'
    varread2 = var+'ToE1'
    ignore1 = [] # Models to ignore for some reason
    ignore2 = []
    color = '#e47f8a' # Color boxes for 2nd set of data
    title1 = 'Hist + RCP8.5 vs. histNat [>2std]'
    title2 = 'Hist + RCP8.5 vs. histNat [>1std]'
    plotName = 'GSATatE_boxplot_' + work +'_'+method_noise_histNat
elif work == 'rcp85_histNat_PiControl_2std':
    indir1 = indir_rcphn
    indir2 = indir_rcppiC
    varread1 = var+'ToE2'
    varread2 = var+'ToE2'
    ignore1 = ['GISS-E2-R','FGOALS-g2','MIROC-ESM']
    ignore2 = ['GISS-E2-R','FGOALS-g2','MIROC-ESM']
    color = '#3D8383'
    title1 = 'Hist + RCP8.5 vs. histNat [>2std]'
    title2 = 'Hist + RCP8.5 vs. PiControl [>2std]'
    plotName = 'GSATatE_boxplot_' + work +'_'+method_noise_histNat+'_'+method_noise_piC
else:
    indir1 = indir_rcphn
    indir2 = indir_CO2piC
    varread1 = var+'ToE2'
    varread2 = var+'ToE2'
    ignore1 = []
    ignore2 = []
    color = '#004f82'
    title1 = 'Hist + RCP8.5 vs. histNat [>2std]'
    title2 = '1pctCO2 vs. PiControl [>2std]'
    plotName = 'GSATatE_boxplot_' + work +'_'+method_noise_histNat+'_'+method_noise_piC

# ----- Read ToE and tas for each model ------

listfiles1 = sorted(glob.glob(indir1 + method_noise_histNat + '/*'+legVar+'_toe_rcp_histNat*.nc'))
nmodels1 = len(listfiles1)-len(ignore1)

if work != 'rcp85_histNat_1_2std':
    listfiles2 = sorted(glob.glob(indir2 + method_noise_piC + '/*.nc'))
    nmodels2=len(listfiles2)
else:
    listfiles2 = listfiles1
    nmodels2 = len(listfiles2)-len(ignore2)


# Read ToE
varToEA_1, varToEP_1, varToEI_1, nMembers1 = read_toe_rcp85(varread1,listfiles1,ignore1,len(domains))
if work != 'rcp85_histNat_1pctCO2_2std':
    varToEA_2, varToEP_2, varToEI_2, nMembers2 = read_toe_rcp85(varread2,listfiles2,ignore2,len(domains))
else:
    varToEA_2, varToEP_2, varToEI_2, nMembers2 = read_toe_1pctCO2(varread2, indir2+'average_piC',modelsCO2,ignore2,len(domains))
    nmodels2 = nmodels2-len(ignore2)

# Read GSAT
gsat_anom1 = read_gsat_rcp85(indir_gsat+'hist-rcp85/',listfiles1,ignore1)
if work != 'rcp85_histNat_1pctCO2_2std':
    gsat_anom2 = gsat_anom1
else:
    gsat_anom2 = read_gsat_1pctCO2(indir_gsat,modelsCO2,ignore2) # Anomaly relative to the last 100 years of piControl

nruns1 = np.sum(nMembers1)
nruns1 = int(nruns1)
nruns2 = int(np.sum(nMembers2))


# ---- Associate ToE to GSAT anomaly ----

maskdata  = np.nan
time1 = np.arange(1861,2101)
if work == 'rcp85_histNat_1pctCO2_2std':
    time2 = np.arange(1,141)
else:
    time2 = time1

# Make new data associating each ToE to its corresponding temperature anomaly
def associate_gsat_ToE(varToE,time,gsat_anom):
    """ Make new data array associating each ToE to its corresponding temperature anomaly """
    ndomains = varToE.shape[1]
    nruns = varToE.shape[0]
    varToE_gsat = np.ma.empty((nruns,ndomains))

    for idomain in range(ndomains): # Loop on regions
        for irun in range(nruns): # Loop on all realizations
            toe = varToE[irun,idomain] # Read ToE
            if not np.ma.is_masked(toe):
                iyear = np.argwhere(time==toe)[0][0] # Read index of said ToE in time vector
                varToE_gsat[irun,idomain] = gsat_anom[iyear,irun] # Fill new data array with gsat anomaly
            else:
                varToE_gsat[irun,idomain] = np.ma.masked

    return varToE_gsat

varToEA_1_gsat = associate_gsat_ToE(varToEA_1,time1,gsat_anom1)
varToEP_1_gsat = associate_gsat_ToE(varToEP_1,time1,gsat_anom1)
varToEI_1_gsat = associate_gsat_ToE(varToEI_1,time1,gsat_anom1)
varToEA_2_gsat = associate_gsat_ToE(varToEA_2,time2,gsat_anom2)
varToEP_2_gsat = associate_gsat_ToE(varToEP_2,time2,gsat_anom2)
varToEI_2_gsat = associate_gsat_ToE(varToEI_2,time2,gsat_anom2)

# -- Turn masked data into nans

#varToEA_1_gsat[np.ma.getmask(varToEA_1_gsat)] = np.nan
#varToEP_1_gsat[np.ma.getmask(varToEP_1_gsat)] = np.nan
#varToEI_1_gsat[np.ma.getmask(varToEI_1_gsat)] = np.nan
#varToEA_2_gsat[np.ma.getmask(varToEA_2_gsat)] = np.nan
#varToEP_2_gsat[np.ma.getmask(varToEP_2_gsat)] = np.nan
#varToEI_2_gsat[np.ma.getmask(varToEI_2_gsat)] = np.nan

# -- Organize data

# New domain labels
new_domains = ['SO subpolar', 'SO subtropics', 'NH subtropics', 'subpolar North Pacific']
# regroup previous "North Atlantic" with NH subtropics

# tas ToE reference hist+rcp8.5 vs. histNat [>2std]
data1 = [varToEA_1_gsat[:,1], varToEP_1_gsat[:,1], varToEI_1_gsat[:,1], maskdata, varToEA_1_gsat[:,0], varToEP_1_gsat[:,0], varToEI_1_gsat[:,0], maskdata, varToEA_1_gsat[:,3], varToEP_1_gsat[:,2], maskdata, varToEP_1_gsat[:,4]]

# tas ToE other case
data2 = [varToEA_2_gsat[:,1], varToEP_2_gsat[:,1], varToEI_2_gsat[:,1], maskdata, varToEA_2_gsat[:,0], varToEP_2_gsat[:,0], varToEI_2_gsat[:,0], maskdata, varToEA_2_gsat[:,3], varToEP_2_gsat[:,2], maskdata, varToEP_2_gsat[:,4]]


# ----- Make pseudo-time vector for gsat_anom1 multi-model mean -----

def movingaverage(interval, window_size):
    window = np.ones(int(window_size))/float(window_size)
    return np.convolve(interval, window, 'same')

# Take gsat_anom multi-model mean
mean_gsat_anom1 = np.ma.average(gsat_anom1,axis=1)

# Smooth to make it a bijection
mean_gsat_anom = mean_gsat_anom1
mean_gsat_anom[0:50] = np.ma.average(mean_gsat_anom[0:50])
gsat_anom_smooth = [0]*240
gsat_anom_smooth[50:] = movingaverage(mean_gsat_anom[50:],60)
gsat_anom_smooth[0:50] = mean_gsat_anom[0:50]
gsat_anom_smooth[-50:] = mean_gsat_anom[-50:]
#plt.plot(mean_gsat_anom1,'k-',lw=1.5,label='gsat anom mme')
#plt.plot(mean_gsat_anom,'b-',lw=1.5)
#plt.plot(gsat_anom_smooth,'r--',lw=1.5,label='with smoothing')
#plt.legend()

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

# Option 1 : align time labels with each half degree, from 0.5 degree to 5.5 degree
newticklocations = np.arange(0.,5.1,0.5) #np.arange(0.0,15.0)
newticknames = ['']*len(newticklocations)
for itick,tickloc in enumerate(newticklocations):
    idx = find_nearest(gsat_anom_smooth,tickloc)
    newticknames[itick] = '%d'%(time1[idx],)

# Option 2 : choose which years to show and find corresponding gsat anomaly
xgsat = np.arange(-1,6.01,0.02)
tickvalues = [1900,1980,2000,2020,2040,2060,2080,2100]
newticklocations2 = [0]*len(tickvalues)
for itick, timeval in enumerate(tickvalues):
    idx = np.argwhere(time1==timeval)[0][0]
    gsatval = gsat_anom_smooth[idx]
    igsattick = find_nearest(xgsat,gsatval)
    newticklocations2[itick] = xgsat[igsattick]
newticknames2 = ['%d' % t for t in tickvalues]


# ----- Plot ------

y1 = 1850
y2 = 1900

labels = ['Atlantic','Pacific','Indian','','Atlantic','Pacific','Indian','','Atlantic','Pacific','','']
N = 13
ind = np.arange(1,N)
width = 0.25

fig, ax = plt.subplots(figsize=(11,13))

ax.axvline(x=1.5, color='black', ls=':')
ax.axvline(x=2, color='black', ls=':')
ax.axvline(x=0, color='black', ls=':')

red_crosses = dict(markeredgecolor='#c90016', marker='+',linewidth=0.5)
# ToE reference boxes
boxes1 = ax.boxplot(data1, vert=0, positions=ind-width, widths=width, whis=0,flierprops=red_crosses)
for box in boxes1['boxes']:
    box.set(color='#c90016', linewidth=2) #c90016 #ad3c48
for whisker in boxes1['whiskers']:
    whisker.set(color='#c90016', linestyle='-', linewidth=1)
for cap in boxes1['caps']:
    cap.set(color='#c90016', linewidth=1)
#for flier in boxes1['fliers']:
#    flier.set(color='#c90016')
for median in boxes1['medians']:
    median.set(color='#c90016', linewidth=2)


ax.set_xlim([-1,6.01])
ax.set_xlabel('GSAT anomaly ('+degree_sign+'C) relative to '+str(y1)+'-'+str(y2), fontweight='bold')
ax.tick_params(axis='y',left=False, right=False, labelright=True, labelleft=False, pad=7)
xmajorLocator = MultipleLocator(0.5)
xminorLocator = AutoMinorLocator(2)
ax.xaxis.set_major_locator(xmajorLocator)
ax.xaxis.set_minor_locator(xminorLocator)
ax.xaxis.set_tick_params(which='major',width=2)

ax2 = ax.twiny()
color_crosses = dict(markeredgecolor=color, marker='+',linewidth=0.5)
# ToE other case boxes
boxes2 = ax2.boxplot(data2, vert=0, positions=ind+width, widths=width, whis=0,flierprops=color_crosses)
for box in boxes2['boxes']:
    box.set(color=color, linewidth=2)
for whisker in boxes2['whiskers']:
    whisker.set(color=color, linestyle='-', linewidth=1)
for cap in boxes2['caps']:
    cap.set(color=color, linewidth=1)
for flier in boxes2['fliers']:
    flier.set(color=color)
for median in boxes2['medians']:
    median.set(color=color, linewidth=2)


ax2.set_xlim([-1,6.01])
ax2.set_yticks(ind)
ax2.set_yticklabels(labels) #, fontweight='bold')
ax2.tick_params(axis='y',labelleft=False,left=False,right=False,labelright=True)
ax2.set_ylim([0,N])
xmajorLocator2 = MultipleLocator(0.5)
xminorLocator2 = AutoMinorLocator(2)
ax2.xaxis.set_major_locator(xmajorLocator2)
ax2.xaxis.set_minor_locator(xminorLocator2)

plt.setp(ax.get_yticklabels(), visible=True)
plt.setp(ax.get_xticklabels(), fontweight='bold')

ax2.axhline(y=ind[3], color='black', ls='--')
ax2.axhline(y=ind[7], color='black', ls='--')
ax2.axhline(y=ind[10], color='black', ls='--')

# Domain labels
ax2.text(-1-0.6,ind[1], 'SO \n subpolar', ha='center', va='center', fontweight='bold', fontsize=13)
ax2.text(-1-0.6,ind[5], 'SO \n subtropics', ha='center', va='center', fontweight='bold', fontsize=13)
ax2.text(-1-0.6,ind[8]+0.5, 'NH \n subtropics', ha='center', va='center', fontweight='bold', fontsize=13)
ax2.text(-1-0.6,ind[11], 'Subpolar \n North Pacific', ha='center', va='center',fontweight='bold', fontsize=13)

#plotTitle = 'Distribution of GSAT at emergence for '+legVar+ ' in different regions'
#ax.set_title(plotTitle, y=1.08, fontweight='bold', va='center')
ax2.text(0.5,1.04, title1 + ' ('+str(nmodels1)+' models, '+str(nruns1)+' runs)', color='#c90016',
         va='center', ha='center',transform=ax2.transAxes, fontweight='bold')
ax2.text(0.5,1.062, title2 + ' ('+str(nmodels2)+' models, '+str(nruns2)+' runs)', color=color,
         va='center', ha='center',transform=ax2.transAxes, fontweight='bold')

# Date
now = datetime.datetime.now()
date = now.strftime("%Y-%m-%d")

# Text at the bottom of the figure
#plt.figtext(.8,.005,'Computed by : boxplot_tas_ToE_rcp85_renamelabels.py, '+date, fontsize=8, ha='center')
#if work == 'rcp85_histNat_1_2std':
#    plt.figtext(.2,.01,'Noise: %s' %(method_noise_histNat), fontsize=8, ha='center')
#else:
#    plt.figtext(.2,.01,'Noise: %s %s' %(method_noise_histNat, method_noise_piC), fontsize=8, ha='center')
#if work == 'rcp85_histNat_PiControl_2std':
#    plt.figtext(.5,.01,'PiControl : mean(last_240_years)',fontsize=8,ha='center')

# -- Make pseudo-time x axis below the original one
ax3 = ax.twiny()
ax3.set_xlim([-1,6.01])
# Add some extra space for the second axis at the bottom
fig.subplots_adjust(bottom=0.1)
# Move twinned axis ticks and label from top to bottom
ax3.xaxis.set_ticks_position("bottom")
ax3.xaxis.set_label_position("bottom")
# Offset the twin axis below the host
ax3.spines["bottom"].set_position(("axes", -0.06))
# Turn on the frame for the twin axis, but then hide all but the bottom spine
ax3.set_frame_on(True)
ax3.patch.set_visible(False)
#for sp in ax3.spines.itervalues():
#    sp.set_visible(False)
ax3.spines["bottom"].set_visible(True)
ax3.set_xticks(newticklocations2)
ax3.set_xticklabels(newticknames2,fontweight='bold',fontsize=10)
ax3.set_xlabel('Pseudo-Years (historical+RCP8.5 mean)',fontweight='bold',fontsize=10)
ax3.xaxis.set_tick_params(which='major',width=2)

if outfmt == 'view':
    plt.show()
else:
    plt.savefig('/home/ysilvy/figures/models/ToE/boxplots/'
                +plotName+'_'+str(y1)+'_'+str(y2)+'_newlabels_paper.png')

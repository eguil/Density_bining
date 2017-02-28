#!/bin/env python
# -*- coding: utf-8 -*-

"""
Python matplotlib
Make time of emergence PDF bar plots from 1%CO2 vs. piControl
Choose which variable and domain to work on

"""

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset as open_ncfile
from maps_matplot_lib import defVarmme, averageDom
from modelsDef import defModelsCO2piC
from libToE import findToE, ToEdomain1pctCO2vsPiC

# ----- Workspace ------

indir_1pctCO2 = '/data/ericglod/Density_binning/Prod_density_april15/mme_1pctCO2/'
indir_piC = '/data/ericglod/Density_binning/Prod_density_april15/mme_piControl/'

models = defModelsCO2piC()

# ----- Work ------

varname = defVarmme('salinity'); v = 'S'
#varname = defVarmme('temp'); v = 'T'
#varname = defVarmme('depth'); v = 'Z'

# -- Choose method for computing ToE
method = 'average_ToE' # Determine 2D lat/rho ToE then average in the box
#method = 'average_signal' # Average signal and noise in the box, then compute ToE (much faster)


domain_name = 'Northern ST'
# 'Southern ST', 'SO', 'Northern ST', 'North Atlantic, 'North Pacific'
print domain_name

multStd = 2. # detect ToE at multStd std dev of piControl

labBowl = ['piControl', '1pctCO2']

valmask = 1.e20

iniyear = 0
finalyear = 140
deltay = 10.

# density domain
rhomin = 21
rhomid = 26
rhomax = 28
domrho = [rhomin, rhomid, rhomax]

# ----- Variables ------

# Choose random file to read only the basic variables and properties common to all files
file = 'cmip5.' + models[0]['name'] + '.1pctCO2.ensm.an.ocn.Omon.density.ver-' + models[0]['file_end_CO2'] + '_zon2D.nc'
f = open_ncfile(indir_1pctCO2 + file,'r')

lat = f.variables['latitude'][:]; latN = lat.size
density = f.variables['lev'][:]; levN = density.size
time = f.variables['time'][:]; timN = time.size
var = varname['var_zonal']

# Define variable properties
minmax = varname['minmax_zonal']
legVar = varname['legVar']
unit = varname['unit']


# ----- Compute ToE for each model ------

# -- Initialize toe
toe_a = np.ma.ones((len(models), levN, latN))*1.
toe_p = np.ma.ones((len(models), levN, latN))*1.
toe_i = np.ma.ones((len(models), levN, latN))*1.
# -- Initialize averaged toe
varToEA = np.ma.masked_all(len(models))
varToEP = np.ma.masked_all(len(models))
varToEI = np.ma.masked_all(len(models))

if method == 'average_ToE' :

    for i, model in enumerate(models):
        print '- Computing ToE of',model['name']

        # -- Read 1pctCO2 file and piControl file
        file_1pctCO2 = 'cmip5.' + model['name'] + '.1pctCO2.ensm.an.ocn.Omon.density.ver-' + model['file_end_CO2'] + '_zon2D.nc'
        file_piC = 'cmip5.' + model['name'] + '.piControl.ensm.an.ocn.Omon.density.ver-' + model['file_end_piC'] + '_zon2D.nc'

        fCO2 = open_ncfile(indir_1pctCO2 + file_1pctCO2,'r')
        fpiC = open_ncfile(indir_piC + file_piC,'r')

        # -- Read var 1pctCO2
        varCO2_a = fCO2.variables[var][:,1,:,:].squeeze()
        varCO2_p = fCO2.variables[var][:,2,:,:].squeeze()
        varCO2_i = fCO2.variables[var][:,3,:,:].squeeze()
        # -- Read var piControl
        varpiC_a = fpiC.variables[var][-140:,1,:,:].squeeze()
        varpiC_p = fpiC.variables[var][-140:,2,:,:].squeeze()
        varpiC_i = fpiC.variables[var][-140:,3,:,:].squeeze()

        # -- Compute significance of difference when diff within 1 stddev of piControl variability (in the MME sense)
        varams = np.ma.std(varpiC_a, axis=0)
        varpms = np.ma.std(varpiC_p, axis=0)
        varims = np.ma.std(varpiC_i, axis=0)

        # -- reorganise i,j dims in single dimension data (speeds up loops)
        varCO2_a  = np.reshape(varCO2_a, (timN,levN*latN))
        varpiC_a = np.reshape(varpiC_a,(timN,levN*latN))
        varams  = np.reshape(varams, (levN*latN))
        varCO2_p  = np.reshape(varCO2_p, (timN,levN*latN))
        varpiC_p = np.reshape(varpiC_p,(timN,levN*latN))
        varpms  = np.reshape(varpms, (levN*latN))
        varCO2_i  = np.reshape(varCO2_i, (timN,levN*latN))
        varpiC_i = np.reshape(varpiC_i,(timN,levN*latN))
        varims  = np.reshape(varims, (levN*latN))

        # -- Compute ToE as last date when diff 1pctCO2 - piControl is larger than mult * stddev
        toei_a = np.reshape(findToE(varCO2_a-varpiC_a, varams, multStd),(levN,latN))
        toei_p = np.reshape(findToE(varCO2_p-varpiC_p, varpms, multStd),(levN,latN))
        toei_i = np.reshape(findToE(varCO2_i-varpiC_i, varims, multStd),(levN,latN))

        # -- Save
        toe_a[i,:,:] = toei_a
        toe_p[i,:,:] = toei_p
        toe_i[i,:,:] = toei_i

        # -- Select domain to average for each model
        domain = ToEdomain1pctCO2vsPiC(model['name'], domain_name)[0]
        domain_char = ToEdomain1pctCO2vsPiC(model['name'], domain_name)[1]

        # -- Average toe
        if domain['Atlantic'] != None:
            varToEA[i] = np.ma.around(averageDom(toe_a[i,:,:], 2, domain['Atlantic'], lat, density))
            print 'ToE Atlantic:',varToEA[i]
        if domain['Pacific'] != None:
            varToEP[i] = np.ma.around(averageDom(toe_p[i,:,:], 2, domain['Pacific'], lat, density))
            print 'ToE Pacific:',varToEP[i]
        if domain['Indian'] != None:
            varToEI[i] = np.ma.around(averageDom(toe_i[i,:,:], 2, domain['Indian'], lat, density))
            print 'ToE Indian:',varToEI[i]


if method == 'average_signal':

    for i, model in enumerate(models):
        print '- Computing ToE of',model['name']

        # -- Read 1pctCO2 file and piControl file
        file_1pctCO2 = 'cmip5.' + model['name'] + '.1pctCO2.ensm.an.ocn.Omon.density.ver-' + model['file_end_CO2'] + '_zon2D.nc'
        file_piC = 'cmip5.' + model['name'] + '.piControl.ensm.an.ocn.Omon.density.ver-' + model['file_end_piC'] + '_zon2D.nc'

        fCO2 = open_ncfile(indir_1pctCO2 + file_1pctCO2,'r')
        fpiC = open_ncfile(indir_piC + file_piC,'r')

        # -- Select domain to average for each model
        domain = ToEdomain1pctCO2vsPiC(model['name'], domain_name)[0]
        domain_char = ToEdomain1pctCO2vsPiC(model['name'], domain_name)[1]

        # -- Read var 1pctCO2
        varCO2_a = fCO2.variables[var][:,1,:,:].squeeze()
        varCO2_p = fCO2.variables[var][:,2,:,:].squeeze()
        varCO2_i = fCO2.variables[var][:,3,:,:].squeeze()
        # -- Read var piControl
        varpiC_a = fpiC.variables[var][-140:,1,:,:].squeeze()
        varpiC_p = fpiC.variables[var][-140:,2,:,:].squeeze()
        varpiC_i = fpiC.variables[var][-140:,3,:,:].squeeze()

        # -- Compute significance of difference when diff within 1 stddev of piControl variability (in the MME sense)
        varstda = np.ma.std(varpiC_a, axis=0)
        varstdp = np.ma.std(varpiC_p, axis=0)
        varstdi = np.ma.std(varpiC_i, axis=0)

        # -- Take difference 1pctCO2 - PiControl
        varsignal_a = varCO2_a - varpiC_a
        varsignal_p = varCO2_p - varpiC_p
        varsignal_i = varCO2_i - varpiC_i

        # -- Average signal and noise
        if domain['Atlantic'] != None:
            varsignal_a = averageDom(varsignal_a, 3, domain['Atlantic'], lat, density)
            varnoise_a = averageDom(varstda, 2, domain['Atlantic'], lat, density)
        if domain['Pacific'] != None:
            varsignal_p = averageDom(varsignal_p, 3, domain['Pacific'], lat, density)
            varnoise_p = averageDom(varstdp, 2, domain['Pacific'], lat, density)
        if domain['Indian'] != None:
            varsignal_i = averageDom(varsignal_i, 3, domain['Indian'], lat, density)
            varnoise_i = averageDom(varstdi, 2, domain['Indian'], lat, density)

        # -- Compute ToE of averaged domain
        if domain['Atlantic'] != None and np.ma.is_masked(varnoise_a) == False:
            toei_a = findToE(varsignal_a, varnoise_a, multStd)
            varToEA[i] = toei_a
        if domain['Pacific'] != None and np.ma.is_masked(varnoise_p) == False:
            toei_p = findToE(varsignal_p, varnoise_p, multStd)
            varToEP[i] = toei_p
        if domain['Indian'] != None and np.ma.is_masked(varnoise_i) == False:
            toei_i = findToE(varsignal_i, varnoise_i, multStd)
            varToEI[i] = toei_i



# -- Compute median ToE
medToEA = np.ma.around(np.ma.median(varToEA))
print(medToEA)
medToEP = np.ma.around(np.ma.median(varToEP))
print(medToEP)
medToEI = np.ma.around(np.ma.median(varToEI))
print(medToEI)


# # ----- Compute ToE for MME ------
#
# if work_toe == 'mme' :
#     print '- Computing MME ToE'
#
#     file_1pctCO2 = 'cmip5.multimodel_piCtl.1pctCO2.ensm.an.ocn.Omon.density_zon2D.nc'
#     file_piC = 'cmip5.multimodel_1pct.piControl.ensm.an.ocn.Omon.density_zon2D.nc'
#
#     fCO2 = open_ncfile(indir_1pctCO2 + file_1pctCO2,'r')
#     fpiC = open_ncfile(indir_piC + file_piC,'r')
#
#     # -- Read var 1pctCO2
#     varCO2_a = fCO2.variables[var][:,1,:,:].squeeze()
#     varCO2_p = fCO2.variables[var][:,2,:,:].squeeze()
#     varCO2_i = fCO2.variables[var][:,3,:,:].squeeze()
#     # -- Read var piControl
#     varpiC_a = fpiC.variables[var][-140:,1,:,:].squeeze()
#     varpiC_p = fpiC.variables[var][-140:,2,:,:].squeeze()
#     varpiC_i = fpiC.variables[var][-140:,3,:,:].squeeze()
#
#     # -- Compute significance of difference when diff within 1 stddev of piControl variability (in the MME sense)
#     varams = np.ma.std(varpiC_a, axis=0)
#     varpms = np.ma.std(varpiC_p, axis=0)
#     varims = np.ma.std(varpiC_i, axis=0)
#
#     # -- reorganise i,j dims in single dimension data (speeds up loops)
#     varCO2_a  = np.reshape(varCO2_a, (timN,levN*latN))
#     varpiC_a = np.reshape(varpiC_a,(timN,levN*latN))
#     varams  = np.reshape(varams, (levN*latN))
#     varCO2_p  = np.reshape(varCO2_p, (timN,levN*latN))
#     varpiC_p = np.reshape(varpiC_p,(timN,levN*latN))
#     varpms  = np.reshape(varpms, (levN*latN))
#     varCO2_i  = np.reshape(varCO2_i, (timN,levN*latN))
#     varpiC_i = np.reshape(varpiC_i,(timN,levN*latN))
#     varims  = np.reshape(varims, (levN*latN))
#
#     # -- Compute ToE as last date when diff 1pctCO2 - piControl is larger than mult * stddev
#     toemme_a = np.reshape(findToE(varCO2_a-varpiC_a, varams, multStd),(levN,latN))
#     toemme_p = np.reshape(findToE(varCO2_p-varpiC_p, varpms, multStd),(levN,latN))
#     toemme_i = np.reshape(findToE(varCO2_i-varpiC_i, varims, multStd),(levN,latN))
#
#
# varmmeToEA = np.around(averageDom(toemme_a, 2, DomToEA['domain'], lat, density))
# varmmeToEP = np.around(averageDom(toemme_p, 2, DomToEA['domain'], lat, density))
# varmmeToEI = np.around(averageDom(toemme_i, 2, DomToEA['domain'], lat, density))


# ----- Plot ------

ndecades = int((finalyear - iniyear)/deltay)
yearbins = np.arange(ndecades+1)*10+iniyear
width = 4

ToEA_bars, bin_edges = np.histogram(varToEA, yearbins)
ToEP_bars, bin_edges = np.histogram(varToEP, yearbins)
ToEI_bars, bin_edges = np.histogram(varToEI, yearbins)
center = (bin_edges[:-1] + bin_edges[1:]) / 2

medToEA_bars, bin_edges = np.histogram(medToEA, yearbins)
medToEP_bars, bin_edges = np.histogram(medToEP, yearbins)
medToEI_bars, bin_edges = np.histogram(medToEI, yearbins)


# -- Function for attaching a label above each bar indicating the number of models
def autolabel(rects, axis):
    """
    Attach a text label above each bar displaying its height
    """
    for rect in rects:
        height = rect.get_height()
        if height > 0:
            axis.text(rect.get_x() + rect.get_width()/2., 1.05*height,
                    '%d' % int(height),
                    ha='center', va='bottom')


nb_basins = domain_char['nb_basins']

fig, ax = plt.subplots(nrows=nb_basins, ncols=1, sharex=True, sharey=True, frameon=False)

# -- Plot according to number of basins for the chosen domain
if nb_basins == 3:

    rects1 = ax[0].bar(center-width, ToEA_bars, width, color ='#87cefa')
    rect1med = ax[0].bar(center, medToEA_bars, width, color ='#f08080')
    ax[0].set_title('Atlantic', fontweight='bold', fontsize=13)
    autolabel(rects1, ax[0])

    rects2 = ax[1].bar(center-width, ToEP_bars, width, color='#87cefa')
    rect2med = ax[1].bar(center, medToEP_bars, width, color ='#f08080')
    ax[1].set_title('Pacific', fontweight='bold', fontsize=13)
    plt.setp(ax[1].get_xticklabels(), visible=True)
    ax[1].xaxis.set_tick_params(width=2)
    autolabel(rects2, ax[1])

    rects3 = ax[2].bar(center-width, ToEI_bars, width, color='#87cefa')
    rect3med = ax[2].bar(center, medToEI_bars, width, color ='#f08080')
    ax[2].set_title('Indian', fontweight='bold', fontsize=13)
    ax[2].xaxis.set_tick_params(width=2)
    autolabel(rects3, ax[2])

    ax = ax[0]

elif nb_basins == 2:
    rects1 = ax[0].bar(center-width, ToEP_bars, width, color ='#87cefa')
    rect1med = ax[0].bar(center, medToEP_bars, width, color ='#f08080')
    ax[0].set_title('Pacific', fontweight='bold', fontsize=13)
    autolabel(rects1, ax[0])

    rects2 = ax[1].bar(center-width, ToEI_bars, width, color='#87cefa')
    rect2med = ax[1].bar(center, medToEI_bars, width, color ='#f08080')
    ax[1].set_title('Indian', fontweight='bold', fontsize=13)
    plt.setp(ax[1].get_xticklabels(), visible=True)
    ax[1].xaxis.set_tick_params(width=2)
    autolabel(rects2, ax[1])

    ax = ax[0]

else :
    if domain_char['Atlantic'] == True:
        rects1 = ax.bar(center-width, ToEA_bars, width, color ='#87cefa')
        rect1med = ax.bar(center, medToEA_bars, width, color ='#f08080')
    else:
        rects1 = ax.bar(center-width, ToEP_bars, width, color ='#87cefa')
        rect1med = ax.bar(center, medToEP_bars, width, color ='#f08080')
    autolabel(rects1, ax)


ax.set_xlim([0,140])
ax.set_ylim([0,7])
ax.set_xticks([0,10,20,30,40,50,60,70,80,90,100,110,120,130,140])
plt.setp(ax.get_xticklabels(), visible=True)
ax.xaxis.set_tick_params(width=2)
ax.legend((rects1[0], rect1med[0]), ('Models', 'Median'), loc='upper left', fontsize=12)

# Keep bottom axis only
for i, axis in enumerate(fig.axes):
    axis.spines['top'].set_visible(False)
    axis.xaxis.set_ticks_position('bottom')
    axis.axes.get_yaxis().set_visible(False)
    axis.spines['left'].set_visible(False); axis.spines['right'].set_visible(False)

# Add a big axes, hide frame, for common labels
fig.add_subplot(111, frameon=False)
# hide tick and tick label of the big axes
plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
plt.xlabel('Years', fontweight='bold')
plt.ylabel('Nb of models per decade', fontweight='bold')


plotTitle = 'ToE (' + legVar + ') in '+domain_name
if nb_basins >1 :
    plt.suptitle(plotTitle, fontweight='bold', fontsize=14, verticalalignment='top')
    plt.subplots_adjust(hspace=.5, wspace=5)
else:
    plt.suptitle(plotTitle, fontweight='bold', fontsize=14)

plt.figtext(.8,.02,'Computed by : bar_toe_1pctCO2vsPiC.py',fontsize=9,ha='center')

#plt.show()

plotName = 'ToE_pdf_' + domain_name + '_' + legVar
#plt.savefig('/home/ysilvy/Density_bining/Yona_analysis/figures/models/ToE/1pctCO2vsPiC/Average_ToE/'+plotName+'.png', bbox_inches='tight')

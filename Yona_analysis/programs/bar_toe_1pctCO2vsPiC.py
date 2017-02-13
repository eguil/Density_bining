#!/bin/env python
# -*- coding: utf-8 -*-

"""
Python matplotlib
Make time of emergence PDF bar plots from 1%CO2 vs. piControl

"""

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset as open_ncfile
from maps_matplot_lib import defVarmme, averageDom
from modelsDef import defModelsCO2piC
from libToE import findToE, ToEdomain

# ----- Workspace ------

indir_1pctCO2 = '/data/ericglod/Density_binning/Prod_density_april15/mme_1pctCO2/'
indir_piC = '/data/ericglod/Density_binning/Prod_density_april15/mme_piControl/'

models = defModelsCO2piC()

# ----- Work ------

varname = defVarmme('salinity'); v = 'S'
#varname = defVarmme('temp'); v = 'T'
#varname = defVarmme('depth'); v = 'Z'

#average = 'same_all'
average = 'each_model'

domain_name = 'Southern ST' #'SO'

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
#clevsm = varname['clevsm_zonal']
#clevsm_bold = varname['clevsm_bold']
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

    if average == 'each_model' :

        # -- Select domain to average for each model
        domain = ToEdomain(model['name'], domain_name)

        # -- Average toe
        if domain['Atlantic'] != None:
            varToEA[i] = np.around(averageDom(toe_a[i,:,:], 2, domain['Atlantic'], lat, density))
            print 'ToE Atlantic:',varToEA[i]
        if domain['Pacific'] != None:
            varToEP[i] = np.around(averageDom(toe_p[i,:,:], 2, domain['Pacific'], lat, density))
            print 'ToE Pacific:',varToEP[i]
        if domain['Indian'] != None:
            varToEI[i] = np.around(averageDom(toe_i[i,:,:], 2, domain['Indian'], lat, density))
            print 'ToE Indian:',varToEI[i]


if 'average' == 'same_all' :

    DomToEA = {'domain': [-40., -20, 25.75, 26.6], 'name': 'Southern ST'}
    DomToEP = {'domain': [-15, -10, 26, 26.3]   , 'name': 'Southern ST'}
    DomToEI = {'domain': [-40, -15, 25.6, 26.8] , 'name': 'Southern ST'}

    # -- Average toe
    varToEA = np.around(averageDom(toe_a, 3, DomToEA['domain'], lat, density))
    print min(varToEA),max(varToEA), np.around(np.average(varToEA)), np.median(varToEA),np.std(varToEA)
    varToEP = np.around(averageDom(toe_p, 3, DomToEP['domain'], lat, density))
    print min(varToEP),max(varToEP), np.around(np.average(varToEP)), np.median(varToEP),np.std(varToEP)
    varToEI = np.around(averageDom(toe_i, 3, DomToEI['domain'], lat, density))
    print min(varToEI),max(varToEI), np.around(np.average(varToEI)), np.median(varToEI),np.std(varToEI)


# -- Compute median ToE
medToEA = np.around(np.ma.median(varToEA))
print(medToEA)
medToEP = np.around(np.ma.median(varToEP))
print(medToEP)
medToEI = np.around(np.ma.median(varToEI))
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


fig, ax = plt.subplots(nrows=3, ncols=1, sharex=True, sharey=True)

rects1 = ax[0].bar(center-width, ToEA_bars, width, color ='#87cefa')
rect1med = ax[0].bar(center, medToEA_bars, width, color ='#f08080')
ax[0].set_xlim([30,140])
ax[0].set_ylim([0,7])
ax[0].set_xticks([30,40,50,60,70,80,90,100,110,120,130,140])
ax[0].set_title('Atlantic', fontweight='bold', fontsize=13)
plt.setp(ax[0].get_xticklabels(), visible=True)
ax[0].xaxis.set_tick_params(width=2)
ax[0].legend((rects1[0], rect1med[0]), ('Models', 'Median'), loc='upper left', fontsize=12)
autolabel(rects1, ax[0])

rects2 = ax[1].bar(center-width, ToEP_bars, width, color='#87cefa')
rect2med = ax[1].bar(center, medToEP_bars, width, color ='#f08080')
ax[1].set_ylabel('Nb of models per decade', fontweight='bold')
ax[1].set_title('Pacific', fontweight='bold', fontsize=13)
plt.setp(ax[1].get_xticklabels(), visible=True)
ax[1].xaxis.set_tick_params(width=2)
autolabel(rects2, ax[1])

rects3 = ax[2].bar(center-width, ToEI_bars, width, color='#87cefa')
rect3med = ax[2].bar(center, medToEI_bars, width, color ='#f08080')
ax[2].set_xlabel('Years', fontweight='bold')
ax[2].set_title('Indian', fontweight='bold', fontsize=13)
ax[2].xaxis.set_tick_params(width=2)
autolabel(rects3, ax[2])

if 'average' == 'same_all':
    method = 'same domain for all models'
else :
    method = 'different domain for each model'
plotTitle = 'ToE (' + legVar + ') in '+domain_name
plt.suptitle(plotTitle, fontweight='bold', fontsize=14, verticalalignment='top')
plt.subplots_adjust(hspace=.5, wspace=5)

plt.show()

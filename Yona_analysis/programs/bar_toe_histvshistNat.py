#!/bin/env python
# -*- coding: utf-8 -*-

"""
Python matplotlib
Make time of emergence PDF bar plots from historical vs. historicalNat
Choose which variable and domain to work on

"""

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset as open_ncfile
from maps_matplot_lib import defVarmme, averageDom
from modelsDef import defModels
from libToE import findToE, ToEdomainhistvshistNat

# ----- Workspace ------

indir_toe = '/data/ericglod/Density_binning/Prod_density_april15/toe_histNat/'

models = defModels()

# ----- Work ------

varname = defVarmme('salinity'); v = 'S'
#varname = defVarmme('temp'); v = 'T'
#varname = defVarmme('depth'); v = 'Z'

# -- Choose method for computing ToE
method = 'average_ToE' # Determine 2D lat/rho ToE then average in the box
#method = 'average_signal' # Average signal and noise in the box, then compute ToE (much faster)


domain_name = 'SO'
# 'Southern ST', 'SO', 'Northern ST', 'North Atlantic, 'North Pacific'
print domain_name

multStd = 2. # detect ToE at multStd std dev of histNat

labBowl = ['histNat', 'hist']

valmask = 1.e20

iniyear = 1860
finalyear = 2005
deltay = 10.

# density domain
rhomin = 21
rhomid = 26
rhomax = 28
domrho = [rhomin, rhomid, rhomax]

# ----- Variables ------

# Choose random file to read only the basic variables and properties common to all files
file = 'cmip5.' + models[0]['name'] + '.historical.ensm.an.ocn.Omon.density.ver-' + models[0]['file_end_hist'] + '_zon2D.nc'
f = open_ncfile(indir_toe + file,'r')

lat = f.variables['latitude'][:]; latN = lat.size
density = f.variables['lev'][:]; levN = density.size
time = f.variables['time'][:]; timN = time.size
var = varname['var_zonal_w/bowl']

# Define variable properties
legVar = varname['legVar']


# ----- Read ToE for each model ------
nruns = 0 # Initialize total number of runs
nrunmax = 100
nMembers = np.ma.empty(len(models)) # Initialize array for keeping nb of members per model
# -- Initialize ToE arrays containing ToE of all runs
toe_a = np.ma.ones((nrunmax, levN, latN))*1.
toe_p = np.ma.ones((nrunmax, levN, latN))*1.
toe_i = np.ma.ones((nrunmax, levN, latN))*1.
# -- Initialize averaged ToE containing ToE of all runs
varToEA = np.ma.masked_all(nrunmax)
varToEP = np.ma.masked_all(nrunmax)
varToEI = np.ma.masked_all(nrunmax)
# -- Initialize median ToE containing medians of each model
medmodelsToEA = np.ma.masked_all(len(models))
medmodelsToEP = np.ma.masked_all(len(models))
medmodelsToEI = np.ma.masked_all(len(models))

# -- Loop over models
for i, model in enumerate(models):

    file_toe = 'cmip5.' + models[i]['name'] + '.historical.ensm.an.ocn.Omon.density.ver-' + models[i]['file_end_hist'] + '_zon2D.nc'

    ftoe = open_ncfile(indir_toe + file_toe, 'r')

    # -- Read ToE (members, basin, density, latitude)
    toe2read = ftoe.variables[var + 'ToE2'][:]
    nMembers[i] = toe2read.shape[0]
    print '- Computing ToE of',model['name'], 'with', nMembers[i], 'members'

    nruns1 = nruns + nMembers[i]
    toe_a[nruns:nruns1,:,:] = toe2read[:,1,:,:]
    toe_p[nruns:nruns1,:,:] = toe2read[:,2,:,:]
    toe_i[nruns:nruns1,:,:] = toe2read[:,3,:,:]

    # -- Select domain to average for each model
    domain = ToEdomainhistvshistNat(model['name'], domain_name)[0]
    domain_char = ToEdomainhistvshistNat(model['name'], domain_name)[1]

    # -- Average toe and determine median of the model
    if domain['Atlantic'] != None:
        varToEA[nruns:nruns1] = np.ma.around(averageDom(toe_a[nruns:nruns1,:,:], 3, domain['Atlantic'], lat, density)) + iniyear
        medmodelsToEA[i] = np.ma.around(np.ma.median(varToEA[nruns:nruns1]))
        print 'Median ToE Atlantic:', medmodelsToEA[i]
    if domain['Pacific'] != None:
        varToEP[nruns:nruns1] = np.ma.around(averageDom(toe_p[nruns:nruns1,:,:], 3, domain['Pacific'], lat, density)) + iniyear
        medmodelsToEP[i] = np.ma.around(np.ma.median(varToEP[nruns:nruns1]))
        print 'Median ToE Pacific:', medmodelsToEP[i]
    if domain['Indian'] != None:
        varToEI[nruns:nruns1] = np.ma.around(averageDom(toe_i[nruns:nruns1,:,:], 3, domain['Indian'], lat, density)) + iniyear
        medmodelsToEI[i] = np.ma.around(np.ma.median(varToEI[nruns:nruns1]))
        print 'Median ToE Indian:', medmodelsToEI[i]


    nruns = nruns1

print 'Total number of runs :', nruns
varToEA = varToEA[0:nruns]
varToEP = varToEP[0:nruns]
varToEI = varToEI[0:nruns]

# -- Determine the median ToE of median ToEs of models (from their different members)
medmedmodelsToEA = np.ma.around(np.ma.median(medmodelsToEA))
medmedmodelsToEP = np.ma.around(np.ma.median(medmodelsToEP))
medmedmodelsToEI = np.ma.around(np.ma.median(medmodelsToEI))

# -- Determine the median of all runs
globalmedToEA = np.ma.around(np.ma.median(varToEA))
globalmedToEP = np.ma.around(np.ma.median(varToEP))
globalmedToEI = np.ma.around(np.ma.median(varToEI))

# ----- Plot ------

# -- Build plot variables (histograms)

ndecades = int((finalyear - iniyear)/deltay)
yearbins = np.arange(ndecades+1)*10+iniyear
yearbins = np.append(yearbins, [finalyear, 2010])
width = 2

ToEA_bars, bin_edges = np.histogram(varToEA, yearbins)
ToEP_bars, bin_edges = np.histogram(varToEP, yearbins)
ToEI_bars, bin_edges = np.histogram(varToEI, yearbins)
center = (bin_edges[:-1] + bin_edges[1:]) / 2

medmodelsToEA_bars, bin_edges = np.histogram(medmodelsToEA, yearbins)
medmodelsToEP_bars, bin_edges = np.histogram(medmodelsToEP, yearbins)
medmodelsToEI_bars, bin_edges = np.histogram(medmodelsToEI, yearbins)

medmedmodelsToEA_bar, bin_edges = np.histogram(medmedmodelsToEA, yearbins)
medmedmodelsToEP_bar, bin_edges = np.histogram(medmedmodelsToEP, yearbins)
medmedmodelsToEI_bar, bin_edges = np.histogram(medmedmodelsToEI, yearbins)

globalmedToEA_bar, bin_edges = np.histogram(globalmedToEA, yearbins)
globalmedToEP_bar, bin_edges = np.histogram(globalmedToEP, yearbins)
globalmedToEI_bar, bin_edges = np.histogram(globalmedToEI, yearbins)


# -- Create variable bundles
varAtl = {'basin': 'Atlantic', 'ToE_bars': ToEA_bars, 'medmodelsToE_bars': medmodelsToEA_bars,
          'medmedmodelsToE_bar': medmedmodelsToEA_bar, 'globalmedToE_bar': globalmedToEA_bar}
varPac = {'basin': 'Pacifc', 'ToE_bars': ToEP_bars, 'medmodelsToE_bars': medmodelsToEP_bars,
          'medmedmodelsToE_bar': medmedmodelsToEP_bar, 'globalmedToE_bar': globalmedToEP_bar}
varInd = {'basin': 'Indian', 'ToE_bars': ToEI_bars, 'medmodelsToE_bars': medmodelsToEI_bars,
          'medmedmodelsToE_bar': medmedmodelsToEI_bar, 'globalmedToE_bar': globalmedToEI_bar}

bundles = [varAtl, varPac, varInd]

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

for i, axis in enumerate(fig.axes):
    print i
    if nb_basins == 3:
        varBasin = bundles[i]

    elif nb_basins == 2:
        varBasin = bundles[i+1]
    elif nb_basins == 1 and domain_name == 'North Pacific':
        varBasin = bundles[1]
    else :
        varBasin = bundles[0]

    # Set median bar to a bigger value for displaying purposes (1 is too small compared to the other bars' height)
    varBasin['medmedmodelsToE_bar'][np.nonzero(varBasin['medmedmodelsToE_bar'])] = 5
    varBasin['globalmedToE_bar'][np.nonzero(varBasin['globalmedToE_bar'])] = 5


    rects = axis.bar(center[:-2]-width*2, varBasin['ToE_bars'][:-2], width, color ='#87cefa')
    rectsbis = axis.bar(center[-2:]-width, varBasin['ToE_bars'][-2:], width, color ='#87cefa')
    globalmed = axis.bar(center[:-2]+width*2, varBasin['globalmedToE_bar'][:-2], width, color ='#1fffaf')
    globalmedbis = axis.bar(center[-2:]-width/2, varBasin['globalmedToE_bar'][-2:], width, color ='#1fffaf')
    medmed = axis.bar(center[:-2]+width*2, varBasin['medmedmodelsToE_bar'][:-2], width, color ='#fafa8c')
    medmedbis = axis.bar(center[-2:], varBasin['medmedmodelsToE_bar'][-2:], width, color ='#fafa8c')
    med = axis.bar(center[:-2]-width, varBasin['medmodelsToE_bars'][:-2], width, color ='#f08080')
    medbis = axis.bar(center[-2:]+width/2, varBasin['medmodelsToE_bars'][-2:], width, color ='#f08080')
    axis.set_title(varBasin['basin'], fontweight='bold', fontsize=13)
    autolabel(rects, axis); autolabel(rectsbis, axis)
    autolabel(med, axis); autolabel(medbis, axis)

    # Set x axis limits, ticks, ticklabels
    axis.set_xlim([1930,2010])
    axis.set_xticks([1930,1940,1950,1960,1970,1980,1990,2000,2005])
    axis.xaxis.set_tick_params(width=2)
    plt.setp(axis.get_xticklabels(), visible=True, rotation=20)

    # Show bottom axis only
    axis.spines['top'].set_visible(False)
    axis.xaxis.set_ticks_position('bottom')
    axis.axes.get_yaxis().set_visible(False)
    axis.spines['left'].set_visible(False); axis.spines['right'].set_visible(False)

    # Set legend
    if i == 0:
        axis.legend((rects[0], med[0], medmed, globalmed), ('Nb of members','Nb of model medians','Median of medians','Median of all members')
                    , loc='upper left', fontsize=11)



# Add a big axes, hide frame, for common labels
fig.add_subplot(111, frameon=False)
# hide tick and tick label of the big axes
plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
plt.xlabel('Years', fontweight='bold')
plt.ylabel('Nb per decade', fontweight='bold')


plotTitle = 'ToE (' + legVar + ') in '+domain_name+ ' (hist vs. histNat)'
if nb_basins >1 :
    plt.suptitle(plotTitle, fontweight='bold', fontsize=14, verticalalignment='top')
    plt.subplots_adjust(hspace=.5, wspace=5)
else:
    plt.suptitle(plotTitle, fontweight='bold', fontsize=14)

plt.figtext(.8,.02,'Computed by : bar_toe_histvshistNat.py',fontsize=9,ha='center')

#plt.show()

plotName = 'ToE_pdf_' + domain_name + '_' + legVar
#plt.savefig('/home/ysilvy/Density_bining/Yona_analysis/figures/models/ToE/histvshistNat/Average_ToE/'+plotName+'.png', bbox_inches='tight')





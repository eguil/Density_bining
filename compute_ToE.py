#!/bin/env python
# -*- coding: utf-8 -*-
import os,glob
from libDensity import defModels,mmeAveMsk2D,mmeAveMsk1D
from string import replace
import warnings

warnings.filterwarnings("ignore")
'''
# ----------------------------------------------------------------------------
#
# Compute ToE 'pythoncd -W ignore compute_ToE.py' (cdms python on mac)
#
# ----------------------------------------------------------------------------

 Compute ToE as last date for which signal = diff hist-histNat > noise = mult * stddev_histNat
    - per model
    - compute the hist-histNat diff for each hist member [rho,lat,time, Nh]
    - compute std dev of histNat: stddev_histNat [rho,lat]
        - average of all Nhn std dev ? TRY THIS FIRST: max of Nhn Stddev ?
    - S/N is (diff hist-histNat) / (mult * stddev_histNat) : ToE[rho,lat,Nh]
    - ToE may be latest date for which S/N > 1 (we then get a PDF and a mean)
    - ToE may be where model agreement on S/N is over 66%
    - How do these two computations compare on the mean ToE ?

Author: Eric Guilyardi - April 2016

'''

# define all models
models = defModels()
nmodels = len(models)

# perform a selection of a few models (for testing or updating)?
modelSel = range(nmodels)
# modelSel = [3,10,18,19,25,27,28]
# modelSel = [22,23,24]

experh = 'historical'
indirh  = '/Users/ericg/Projets/Density_bining/Prod_density_april15/historical'
experhn = 'historicalNat'
indirhn  = '/Users/ericg/Projets/Density_bining/Prod_density_april15/historicalNat'

outdir = '/Users/ericg/Projets/Density_bining/Prod_density_april15/ToE'

print
print '-------------------------------------------------------------------------'
print ' Enter compute_ToE.py to compute ToE'

for i in modelSel:
    # model properties and number of ensembles
    mod = models[i]['name']
    nensh = models[i]['props'][0]
    nenshn = models[i]['props'][1]
    years = [models[i]['props'][2],models[i]['props'][3]]
    if years[1] <> 0: # do not ignore model
        if nensh > 0: # only if 1 member or more
            listf  = glob.glob('cmip5.'+mod+'.*zon2D*')
            listf1 = glob.glob('cmip5.'+mod+'.*zon1D*')
            start = listf[0].find(exper)+len(exper)
            end = listf[0].find('.an.')
            rip = listf[0][start:end]
            outFile = replace(listf[0],rip,'.ensmToE')
            outFile1 = replace(outFile,'2D','1D')
            print ' -> working on: ', i,mod, 'slice', years, nensh, 'members'
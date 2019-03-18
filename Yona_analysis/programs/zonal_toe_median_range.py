#!/bin/env python
# -*- coding: utf-8 -*-

"""
Python matplotlib
Plot density/latitude Time of Emergence (median and distribution range)

"""

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset as open_ncfile
from maps_matplot_lib import defVarmme, zonal_2D
from modelsDef import defModels, defModelsCO2piC
import glob
import os
import colormaps as cmaps
import cmocean
import palettable.cartocolors.sequential
import datetime


# ----- Workspace ------

indir_toe_rcphn = '/home/ysilvy/Density_bining/Yona_analysis/data/toe_zonal/toe_rcp85_histNat/'
indir_toe_rcppiC = '/home/ysilvy/Density_bining/Yona_analysis/data/toe_zonal/toe_rcp85_PiControl/'
indir_mme_rcp85 = '/data/ericglod/Density_binning/Prod_density_april15/mme_rcp85/'
indir_mme_hn = '/data/ericglod/Density_binning/Prod_density_april15/mme_histNat/'
indir_toe_CO2piC = '/home/ysilvy/Density_bining/Yona_analysis/data/toe_zonal/toe_1pctCO2_piC/'
indir_mme_CO2 = '/data/ericglod/Density_binning/Prod_density_april15/mme_1pctCO2/'
indir_mme_piC = '/data/ericglod/Density_binning/Prod_density_april15/mme_piControl/'

# ----- Work ------

# === INPUTS ===
work = 'RCP85'
# work = 'CO2'

figure = 'median'
# figure = 'range'

# output format
outfmt = 'view'
# outfmt = 'save'

use_piC = False # Signal = (hist-histNat) + RCP8.5-average(histNat), noise = std(histNat)
# use_piC = True # Signal = hist + RCP8.5 - average(PiControl), noise = std(PiControl)

# runs_rcp = 'same' # Same runs (30 runs) for hist+RCP8.5 vs. histNat as for hist+RCP8.5 vs. PiControl
runs_rcp = 'all' # All runs (35)

varname = defVarmme('salinity'); v = 'S'

multstd = 2 # detect ToE at multstd std dev of histnat/PiControl

nb_outliers = 5 # Nb of outlier runs allowed to compute the ToE distribution for hist + RCP8.5
nb_outliers_CO2 = 1 # Nb of outlier runs allowed to compute the ToE distribution for 1%CO2 vs. PiControl

# =====

if work == 'RCP85':
    iniyear = 1860
    finalyear = 2100
    models = defModels()
    min = 1960 # Start year for colorbar
    crit_min_runs = 5 # min number of runs with unmasked values = nruns-crit_min_runs
else:
    iniyear = 0
    finalyear = 140
    models = defModelsCO2piC()
    min = 20 # Start year for colorbar
    crit_min_runs = 1 # min number of runs with unmasked values = nruns-crit_min_runs
deltat = 10.

# density domain
rhomin = 21
rhomid = 26
rhomax = 28
domrho = [rhomin, rhomid, rhomax]

# ----- Variables ------
var = varname['var_zonal_w/bowl']
legVar = varname['legVar']
unit = varname['unit']

# Read latitude and density from original file
fileh_2d = '/data/ericglod/Density_binning/Prod_density_april15/mme_hist/' \
       'cmip5.multimodel_Nat.historical.ensm.an.ocn.Omon.density_zon2D.nc'
fh2d = open_ncfile(fileh_2d, 'r')
lat = fh2d.variables['latitude'][:]; latN = len(lat)
density = fh2d.variables['lev'][:]; levN = len(density)


# ----- Read ToE for each model ------

if work == 'RCP85':
    # == Historical + RCP8.5 vs. historicalNat or vs. PiControl ==

    nruns = 0 # Initialize total number of runs
    nrunmax = 100
    nMembers = np.ma.empty(len(models)) # Initialize array for keeping number of members per model

    # -- Initialize varToE containing ToE of all runs
    varToEA = np.ma.masked_all((nrunmax, levN, latN))
    varToEP = np.ma.masked_all((nrunmax, levN, latN))
    varToEI = np.ma.masked_all((nrunmax, levN, latN))

    # -- Initialize varsignal (essential for knowing the sign of the emergence)
    varsignal_a = np.ma.masked_all((nrunmax, levN, latN))
    varsignal_p = np.ma.masked_all((nrunmax, levN, latN))
    varsignal_i = np.ma.masked_all((nrunmax, levN, latN))

    # Loop over models
    if use_piC == False:
        listfiles = glob.glob(indir_toe_rcphn + '/*.nc')
    else:
        listfiles = glob.glob(indir_toe_rcppiC + '/*.nc')
    nmodels = len(listfiles)

    for i in range(nmodels):

        file_toe = listfiles[i]
        ftoe = open_ncfile(file_toe, 'r')
        name = os.path.basename(file_toe).split('.')[1]

        # If use same runs in vs. histNat as in vs. PiControl, take out deficient models
        if (runs_rcp == 'all') or (runs_rcp =='same' and name != 'GISS-E2-R' and name != 'FGOALS-g2' and name != 'MIROC-ESM'):

            # Read ToE (members, basin, density, latitude)
            if multstd == 1:
                toeread = ftoe.variables[var + 'ToE1'][:]
            else:
                toeread = ftoe.variables[var + 'ToE2'][:]
            nMembers[i] = toeread.shape[0]
            print('- Reading ToE of %s with %d members'%(name,nMembers[i]))
            nruns1 = int(nruns + nMembers[i])

            # Save ToE
            varToEA[nruns:nruns1,:,:] = toeread[:,1,:,:]
            varToEP[nruns:nruns1,:,:] = toeread[:,2,:,:]
            varToEI[nruns:nruns1,:,:] = toeread[:,3,:,:]

            # Read signal
            signalread = ftoe.variables[var + '_change'][:]

            # Save signal
            varsignal_a[nruns:nruns1,:,:] = signalread[:,1,:,:]
            varsignal_p[nruns:nruns1,:,:] = signalread[:,2,:,:]
            varsignal_i[nruns:nruns1,:,:] = signalread[:,3,:,:]

            nruns = nruns1

    print('Total number of runs:', nruns)
    varToEA = varToEA[0:nruns,:,:]
    varToEP = varToEP[0:nruns,:,:]
    varToEI = varToEI[0:nruns,:,:]
    varsignal_a = varsignal_a[0:nruns,:,:]
    varsignal_p = varsignal_p[0:nruns,:,:]
    varsignal_i = varsignal_i[0:nruns,:,:]

    nruns = int(nruns)

    if runs_rcp == 'same':
        nmodels=nmodels-3

else:
    # == 1%CO2 vs. PiControl ==

    listfiles = glob.glob(indir_toe_CO2piC + '/*.nc')
    nmodels = len(listfiles)

    # -- Initialize varToE containing ToE of all runs
    varToEA = np.ma.masked_all((nmodels, levN, latN))
    varToEP = np.ma.masked_all((nmodels, levN, latN))
    varToEI = np.ma.masked_all((nmodels, levN, latN))

    # -- Initialize varsignal (essential for knowing the sign of the emergence)
    varsignal_a = np.ma.masked_all((nmodels, levN, latN))
    varsignal_p = np.ma.masked_all((nmodels, levN, latN))
    varsignal_i = np.ma.masked_all((nmodels, levN, latN))

    # Loop over models
    for i in range(nmodels):

        file_toe = listfiles[i]
        ftoe = open_ncfile(file_toe, 'r')
        name = os.path.basename(file_toe).split('.')[1]
        # Read ToE (basin, density, latitude)
        if multstd == 1:
            toeread = ftoe.variables[var + 'ToE1'][:]
        else:
            toeread = ftoe.variables[var + 'ToE2'][:]
        print('- Reading ToE of %s'%(name,))

        # Save ToE
        varToEA[i,:,:] = toeread[1,:,:]
        varToEP[i,:,:] = toeread[2,:,:]
        varToEI[i,:,:] = toeread[3,:,:]

        # Read signal
        signalread = ftoe.variables[var + '_change'][:]

        # Save signal
        varsignal_a[i,:,:] = signalread[1,:,:]
        varsignal_p[i,:,:] = signalread[2,:,:]
        varsignal_i[i,:,:] = signalread[3,:,:]

    print('Total number of models:', nmodels)

    nruns = nmodels


# # -- Compute median and range -- old method
# # Median
# medianToEA = np.ma.around(np.ma.median(varToEA, axis=0)) + iniyear
# medianToEP = np.ma.around(np.ma.median(varToEP, axis=0)) + iniyear
# medianToEI = np.ma.around(np.ma.median(varToEI, axis=0)) + iniyear
# # 16th percentile
# percentile16ToEA = np.percentile(varToEA, 16, axis=0)
# percentile16ToEP = np.percentile(varToEP, 16, axis=0)
# percentile16ToEI = np.percentile(varToEI, 16, axis=0)
# # 84th percentile
# percentile84ToEA = np.percentile(varToEA, 84, axis=0)
# percentile84ToEP = np.percentile(varToEP, 84, axis=0)
# percentile84ToEI = np.percentile(varToEI, 84, axis=0)
# # 16-84% range
# rangeToEA = np.ma.around(percentile84ToEA - percentile16ToEA)
# rangeToEP = np.ma.around(percentile84ToEP - percentile16ToEP)
# rangeToEI = np.ma.around(percentile84ToEI - percentile16ToEI)
#
#
# # -- Mask points
# rangeToEA[rangeToEA == 0] = np.ma.masked
# rangeToEP[rangeToEP == 0] = np.ma.masked
# rangeToEI[rangeToEI == 0] = np.ma.masked
#
# rangeToEA[medianToEA > finalyear-20] = np.ma.masked # Mask points where median hasn't emerged

# ----- Compute distribution according to Lyu et. al method ------

groups_a = np.ma.masked_all((nruns, levN, latN))
groups_p = np.ma.masked_all((nruns, levN, latN))
groups_i = np.ma.masked_all((nruns, levN, latN))
# Groups values : 0 for no emergence, 1 emergence with positive change, -1 emergence with negative change

# -- Check if emergence

groups_a = np.ma.where(varsignal_a>0,1,-1) # 1 if positive change, -1 if negative change
groups_a[np.ma.where(varToEA > finalyear-20)] = 0 # Overwrite with 0 if no emergence
groups_p = np.ma.where(varsignal_p>0,1,-1)
groups_p[np.ma.where(varToEP > finalyear-20)] = 0
groups_i = np.ma.where(varsignal_i>0,1,-1)
groups_i[np.ma.where(varToEI > finalyear-20)] = 0

nruns_a = np.ma.count(varsignal_a,axis=0) # Nb of runs with unmasked values
nruns_p = np.ma.count(varsignal_p,axis=0)
nruns_i = np.ma.count(varsignal_i,axis=0)

# -- Initialize median and percentiles
medianToEA = np.ma.masked_all((levN, latN))
medianToEP = np.ma.masked_all((levN, latN))
medianToEI = np.ma.masked_all((levN, latN))
percentile25ToEA = np.ma.masked_all((levN, latN))
percentile25ToEP = np.ma.masked_all((levN, latN))
percentile25ToEI = np.ma.masked_all((levN, latN))
percentile75ToEA = np.ma.masked_all((levN, latN))
percentile75ToEP = np.ma.masked_all((levN, latN))
percentile75ToEI = np.ma.masked_all((levN, latN))

# -- Initialize where there is no agreement between models on the sign of the signal
noagree_a = np.ma.masked_all((levN, latN))
noagree_p = np.ma.masked_all((levN, latN))
noagree_i = np.ma.masked_all((levN, latN))


# -- Compute distribution if conditions are met 

varToEA_clean = np.ma.masked_all((nruns,levN,latN))
varToEP_clean = np.ma.masked_all((nruns,levN,latN))
varToEI_clean = np.ma.masked_all((nruns,levN,latN))

# TODO maybe use reshape instead of double loop ? Check feasability --> idem TOE computation
for ilat in range(len(lat)):
    for isig in range(len(density)):

        # -- Atlantic
        if nruns_a[isig,ilat] >= nruns - crit_min_runs :
            # Positive emergence
            if ((groups_a[:,isig,ilat]==1).sum() > nruns_a[isig,ilat]/2) and ((groups_a[:,isig,ilat]==-1).sum() < nb_outliers):
                varToEA_clean[np.where(groups_a[:,isig,ilat]!=-1),isig,ilat] = varToEA[np.where(groups_a[:,isig,ilat]!=-1),isig,ilat]
                # print(varToEA_clean[:,isig,ilat])
            # Negative emergence
            elif ((groups_a[:,isig,ilat]==-1).sum() > nruns_a[isig,ilat]/2) and ((groups_a[:,isig,ilat]==1).sum() < nb_outliers):
                 varToEA_clean[np.where(groups_a[:,isig,ilat]!=1),isig,ilat] = varToEA[np.where(groups_a[:,isig,ilat]!=1),isig,ilat]
                 # print(varToEA_clean[:,isig,ilat])
            # No emergence
            elif ((groups_a[:,isig,ilat]==0).sum() > nruns_a[isig,ilat]/2) and ((groups_a[:,isig,ilat]==1).sum()<=nb_outliers or
                    (groups_a[:,isig,ilat]==-1).sum()<=nb_outliers):
                if (groups_a[:,isig,ilat]==1).sum() > (groups_a[:,isig,ilat]==-1).sum() :
                    varToEA_clean[np.where(groups_a[:,isig,ilat]!=-1),isig,ilat] = varToEA[np.where(groups_a[:,isig,ilat]!=-1),isig,ilat]
                else:
                    varToEA_clean[np.where(groups_a[:,isig,ilat]!=1),isig,ilat] = varToEA[np.where(groups_a[:,isig,ilat]!=1),isig,ilat]
            # No agreement between models
            else:
                noagree_a[isig,ilat] = 1
                varToEA_clean[:,isig,ilat] = varToEA[:,isig,ilat]

        # -- Pacific
        if nruns_p[isig,ilat] >= nruns - crit_min_runs :
            # Positive emergence
            if ((groups_p[:,isig,ilat]==1).sum() > nruns_p[isig,ilat]/2) and ((groups_p[:,isig,ilat]==-1).sum() < nb_outliers):
                varToEP_clean[np.where(groups_p[:,isig,ilat]!=-1),isig,ilat] = varToEP[np.where(groups_p[:,isig,ilat]!=-1),isig,ilat]
            # Negative emergence
            elif ((groups_p[:,isig,ilat]==-1).sum() > nruns_p[isig,ilat]/2) and ((groups_p[:,isig,ilat]==1).sum() < nb_outliers):
                 varToEP_clean[np.where(groups_p[:,isig,ilat]!=1),isig,ilat] = varToEP[np.where(groups_p[:,isig,ilat]!=1),isig,ilat]
            # No emergence
            elif ((groups_p[:,isig,ilat]==0).sum() > nruns_p[isig,ilat]/2) and ((groups_p[:,isig,ilat]==1).sum()<=nb_outliers or
                    (groups_p[:,isig,ilat]==-1).sum()<=nb_outliers):
                if (groups_p[:,isig,ilat]==1).sum() > (groups_p[:,isig,ilat]==-1).sum() :
                    varToEP_clean[np.where(groups_p[:,isig,ilat]!=-1),isig,ilat] = varToEP[np.where(groups_p[:,isig,ilat]!=-1),isig,ilat]
                else:
                    varToEP_clean[np.where(groups_p[:,isig,ilat]!=1),isig,ilat] = varToEP[np.where(groups_p[:,isig,ilat]!=1),isig,ilat]
            # No agreement between models
            else:
                noagree_p[isig,ilat] = 1
                varToEP_clean[:,isig,ilat] = varToEP[:,isig,ilat]

        # -- Indian
        if nruns_i[isig,ilat] >= nruns - crit_min_runs :
            # Positive emergence
            if ((groups_i[:,isig,ilat]==1).sum() > nruns_i[isig,ilat]/2) and ((groups_i[:,isig,ilat]==-1).sum() < nb_outliers):
                varToEI_clean[np.where(groups_i[:,isig,ilat]!=-1),isig,ilat] = varToEI[np.where(groups_i[:,isig,ilat]!=-1),isig,ilat]
            # Negative emergence
            elif ((groups_i[:,isig,ilat]==-1).sum() > nruns_i[isig,ilat]/2) and ((groups_i[:,isig,ilat]==1).sum() < nb_outliers):
                 varToEI_clean[np.where(groups_i[:,isig,ilat]!=1),isig,ilat] = varToEI[np.where(groups_i[:,isig,ilat]!=1),isig,ilat]
            # No emergence
            elif ((groups_i[:,isig,ilat]==0).sum() > nruns_i[isig,ilat]/2) and ((groups_i[:,isig,ilat]==1).sum()<=nb_outliers or
                    (groups_i[:,isig,ilat]==-1).sum()<=nb_outliers):
                if (groups_i[:,isig,ilat]==1).sum() > (groups_i[:,isig,ilat]==-1).sum() :
                    varToEI_clean[np.where(groups_i[:,isig,ilat]!=-1),isig,ilat] = varToEI[np.where(groups_i[:,isig,ilat]!=-1),isig,ilat]
                else:
                    varToEI_clean[np.where(groups_i[:,isig,ilat]!=1),isig,ilat] = varToEI[np.where(groups_i[:,isig,ilat]!=1),isig,ilat]
            # No agreement between models
            else:
                noagree_i[isig,ilat] = 1
                varToEI_clean[:,isig,ilat] = varToEI[:,isig,ilat]

# THOUGHT : deep Atlantic blank due to not enough unmasked runs ?

print('Loops done')
medianToEA = np.ma.around(np.ma.median(varToEA_clean, axis=0)) + iniyear
percentile25ToEA = np.ma.around(np.percentile(varToEA_clean, 25, axis=0)) + iniyear
percentile75ToEA = np.ma.around(np.percentile(varToEA_clean, 75, axis=0)) + iniyear
medianToEP = np.ma.around(np.ma.median(varToEP_clean, axis=0)) + iniyear
percentile25ToEP = np.ma.around(np.percentile(varToEP_clean, 25, axis=0)) + iniyear
percentile75ToEP = np.ma.around(np.percentile(varToEP_clean, 75, axis=0)) + iniyear
medianToEI = np.ma.around(np.ma.median(varToEI_clean, axis=0)) + iniyear
percentile25ToEI = np.ma.around(np.percentile(varToEI_clean, 25, axis=0)) + iniyear
percentile75ToEI = np.ma.around(np.percentile(varToEI_clean, 75, axis=0)) + iniyear


# for ilat in range(len(lat)):
#     for isig in range(len(density)):
#
#         # Atlantic
#         if nruns_a[isig,ilat] >= nruns - 5 : # Criteria for RCP8.5 (lower for 1%CO2)
#             if ((groups_a[:,isig,ilat]==1).sum() > nruns_a[isig,ilat]/2) and ((groups_a[:,isig,ilat]==-1).sum() < nb_outliers):
#                 varToEA_clean = varToEA[np.where(groups_a[:,isig,ilat]!=-1),isig,ilat]
#                 medianToEA[isig,ilat] = np.ma.around(np.ma.median(varToEA_clean)) + iniyear
#                 percentile16ToEA[isig,ilat] = np.ma.around(np.percentile(varToEA_clean, 16)) + iniyear
#                 percentile84ToEA[isig,ilat] = np.ma.around(np.percentile(varToEA_clean, 84)) + iniyear
#             elif ((groups_a[:,isig,ilat]==-1).sum() > nruns_a[isig,ilat]/2) and ((groups_a[:,isig,ilat]==1).sum() < nb_outliers):
#                 varToEA_clean = varToEA[np.where(groups_a[:,isig,ilat]!=1),isig,ilat]
#                 medianToEA[isig,ilat] = np.ma.around(np.ma.median(varToEA_clean)) + iniyear
#                 percentile16ToEA[isig,ilat] = np.ma.around(np.percentile(varToEA_clean, 16)) + iniyear
#                 percentile84ToEA[isig,ilat] = np.ma.around(np.percentile(varToEA_clean, 84)) + iniyear
#             elif ((groups_a[:,isig,ilat]==0).sum() > nruns_a[isig,ilat]/2) and ((groups_a[:,isig,ilat]==1).sum()<=nb_outliers or
#                     (groups_a[:,isig,ilat]==-1).sum()<=-1):
#                 if (groups_a[:,isig,ilat]==1).sum() > (groups_a[:,isig,ilat]==-1).sum() :
#                     varToEA_clean = varToEA[np.where(groups_a[:,isig,ilat]!=-1),isig,ilat]
#                 else:
#                     varToEA_clean = varToEA[np.where(groups_a[:,isig,ilat]!=1),isig,ilat]
#                 medianToEA[isig,ilat] = np.ma.around(np.ma.median(varToEA_clean)) + iniyear
#                 percentile16ToEA[isig,ilat] = np.ma.around(np.percentile(varToEA_clean, 16)) + iniyear
#                 percentile84ToEA[isig,ilat] = np.ma.around(np.percentile(varToEA_clean, 84)) + iniyear
#             else:
#                 noagree_a[isig,ilat] = 1
#                 medianToEA[isig,ilat] = np.ma.around(np.ma.median(varToEA)) + iniyear
#                 percentile16ToEA[isig,ilat] = np.ma.around(np.percentile(varToEA, 16)) + iniyear
#                 percentile84ToEA[isig,ilat] = np.ma.around(np.percentile(varToEA, 84)) + iniyear
#
#         # Pacific
#         if nruns_p[isig,ilat] >= nruns - 5:
#             if ((groups_p[:,isig,ilat]==1).sum() > nruns_p[isig,ilat]/2) and ((groups_p[:,isig,ilat]==-1).sum() < nb_outliers):
#                 varToEP_clean = varToEP[np.where(groups_p[:,isig,ilat]!=-1),isig,ilat]
#                 medianToEP[isig,ilat] = np.ma.around(np.ma.median(varToEP_clean)) + iniyear
#                 percentile16ToEP[isig,ilat] = np.ma.around(np.percentile(varToEP_clean, 16)) + iniyear
#                 percentile84ToEP[isig,ilat] = np.ma.around(np.percentile(varToEP_clean, 84)) + iniyear
#             elif ((groups_p[:,isig,ilat]==-1).sum() > nruns_p[isig,ilat]/2) and ((groups_p[:,isig,ilat]==1).sum() < nb_outliers):
#                 varToEP_clean = varToEP[np.where(groups_p[:,isig,ilat]!=1),isig,ilat]
#                 medianToEP[isig,ilat] = np.ma.around(np.ma.median(varToEP_clean)) + iniyear
#                 percentile16ToEP[isig,ilat] = np.ma.around(np.percentile(varToEP_clean, 16)) + iniyear
#                 percentile84ToEP[isig,ilat] = np.ma.around(np.percentile(varToEP_clean, 84)) + iniyear
#             elif ((groups_p[:,isig,ilat]==0).sum() > nruns_p[isig,ilat]/2) and ((groups_p[:,isig,ilat]==1).sum()<=nb_outliers or
#                     (groups_p[:,isig,ilat]==-1).sum()<=-1):
#                 if (groups_p[:,isig,ilat]==1).sum() > (groups_p[:,isig,ilat]==-1).sum() :
#                     varToEP_clean = varToEP[np.where(groups_p[:,isig,ilat]!=-1),isig,ilat]
#                 else:
#                     varToEP_clean = varToEP[np.where(groups_p[:,isig,ilat]!=1),isig,ilat]
#                 medianToEP[isig,ilat] = np.ma.around(np.ma.median(varToEP_clean)) + iniyear
#                 percentile16ToEP[isig,ilat] = np.ma.around(np.percentile(varToEP_clean, 16)) + iniyear
#                 percentile84ToEP[isig,ilat] = np.ma.around(np.percentile(varToEP_clean, 84)) + iniyear
#             else:
#                 noagree_p[isig,ilat] = 1
#                 medianToEP[isig,ilat] = np.ma.around(np.ma.median(varToEP)) + iniyear
#                 percentile16ToEP[isig,ilat] = np.ma.around(np.percentile(varToEP, 16)) + iniyear
#                 percentile84ToEP[isig,ilat] = np.ma.around(np.percentile(varToEP, 84)) + iniyear
#
#         # Indian
#         if nruns_i[isig,ilat] >= nruns - 5:
#             if ((groups_i[:,isig,ilat]==1).sum() > nruns_i[isig,ilat]/2) and ((groups_i[:,isig,ilat]==-1).sum() < nb_outliers):
#                 varToEI_clean = varToEI[np.where(groups_i[:,isig,ilat]!=-1),isig,ilat]
#                 medianToEI[isig,ilat] = np.ma.around(np.ma.median(varToEI_clean)) + iniyear
#                 percentile16ToEI[isig,ilat] = np.ma.around(np.percentile(varToEI_clean, 16)) + iniyear
#                 percentile84ToEI[isig,ilat] = np.ma.around(np.percentile(varToEI_clean, 84)) + iniyear
#             elif ((groups_i[:,isig,ilat]==-1).sum() > nruns_i[isig,ilat]/2) and ((groups_i[:,isig,ilat]==1).sum() < nb_outliers):
#                 varToEI_clean = varToEI[np.where(groups_i[:,isig,ilat]!=1),isig,ilat]
#                 medianToEI[isig,ilat] = np.ma.around(np.ma.median(varToEI_clean)) + iniyear
#                 percentile16ToEI[isig,ilat] = np.ma.around(np.percentile(varToEI_clean, 16)) + iniyear
#                 percentile84ToEI[isig,ilat] = np.ma.around(np.percentile(varToEI_clean, 84)) + iniyear
#             elif ((groups_i[:,isig,ilat]==0).sum() > nruns_i[isig,ilat]/2) and ((groups_i[:,isig,ilat]==1).sum()<=nb_outliers or
#                     (groups_i[:,isig,ilat]==-1).sum()<=-1):
#                 if (groups_i[:,isig,ilat]==1).sum() > (groups_i[:,isig,ilat]==-1).sum() :
#                     varToEI_clean = varToEI[np.where(groups_i[:,isig,ilat]!=-1),isig,ilat]
#                 else:
#                     varToEI_clean = varToEI[np.where(groups_i[:,isig,ilat]!=1),isig,ilat]
#                 medianToEI[isig,ilat] = np.ma.around(np.ma.median(varToEI_clean)) + iniyear
#                 percentile16ToEI[isig,ilat] = np.ma.around(np.percentile(varToEI_clean, 16)) + iniyear
#                 percentile84ToEI[isig,ilat] = np.ma.around(np.percentile(varToEI_clean, 84)) + iniyear
#             else:
#                 noagree_i[isig,ilat] = 1
#                 medianToEI[isig,ilat] = np.ma.around(np.ma.median(varToEI)) + iniyear
#                 percentile16ToEI[isig,ilat] = np.ma.around(np.percentile(varToEI, 16)) + iniyear
#                 percentile84ToEI[isig,ilat] = np.ma.around(np.percentile(varToEI, 84)) + iniyear
#
# print('Loops done')

# 16-84% range
rangeToEA = percentile75ToEA - percentile25ToEA
rangeToEP = percentile75ToEP - percentile25ToEP
rangeToEI = percentile75ToEI - percentile25ToEI

rangeToEA[percentile75ToEA>finalyear-20] = np.ma.masked
norangeA = np.where(percentile75ToEA>finalyear-20,1,0)
rangeToEP[percentile75ToEP>finalyear-20] = np.ma.masked
norangeP = np.where(percentile75ToEP>finalyear-20,1,0)
rangeToEI[percentile75ToEI>finalyear-20] = np.ma.masked
norangeI = np.where(percentile75ToEI>finalyear-20,1,0)

# Mask
rangeToEA[rangeToEA==0] = np.ma.masked
rangeToEP[rangeToEP==0] = np.ma.masked
rangeToEI[rangeToEI==0] = np.ma.masked

# ----- Save to file ------
# TODO ?

# ----- Read bowl position and mask points above ------

if work == 'RCP85':
    # Read files
    file_rcp85 = 'cmip5.multimodel_Nat.rcp85.ensm.an.ocn.Omon.density_zon1D.nc'
    file_hn = 'cmip5.multimodel_Nat.historicalNat.ensm.an.ocn.Omon.density_zon1D.nc'
    file_piC = 'cmip5.multimodel_1pct.piControl.ensm.an.ocn.Omon.density_zon1D.nc'
    f2 = open_ncfile(indir_mme_rcp85+file_rcp85,'r')
    if use_piC:
        f1 = open_ncfile(indir_mme_piC+file_piC,'r')
        labBowl = ['PiControl','RCP8.5']
    else:
        f1 = open_ncfile(indir_mme_hn+file_hn,'r')
        labBowl = ['histNat','RCP8.5']

else:
    file_CO2 = 'cmip5.multimodel_piCtl.1pctCO2.ensm.an.ocn.Omon.density_zon1D.nc'
    file_piC = 'cmip5.multimodel_1pct.piControl.ensm.an.ocn.Omon.density_zon1D.nc'
    f2 = open_ncfile(indir_mme_CO2+file_CO2,'r')
    f1 = open_ncfile(indir_mme_piC+file_piC,'r')
    labBowl = ['PiControl','4*CO2']

# Read bowl position
bowl2 = f2.variables['ptopsigma'][-5:,:,:]
bowl1 = f1.variables['ptopsigma'][-5:,:,:]
bowl2 = np.ma.average(bowl2, axis=0)
bowl1 = np.ma.average(bowl1, axis=0)

# Mask points above RCP8.5/1%CO2 bowl
for ilat in range(len(lat)):
    if np.ma.is_masked(bowl2[1,ilat]) == False :
        inda = np.ma.nonzero(bowl2[1,ilat]>=density)
        medianToEA[inda,ilat] = np.ma.masked
        rangeToEA[inda,ilat] = np.ma.masked
        noagree_a[inda,ilat] = np.ma.masked
        norangeA[inda,ilat] = np.ma.masked
    if np.ma.is_masked(bowl2[2,ilat]) == False :
        indp = np.ma.nonzero(bowl2[2,ilat]>=density)
        medianToEP[indp,ilat] = np.ma.masked
        rangeToEP[indp,ilat] = np.ma.masked
        noagree_p[indp,ilat] = np.ma.masked
        norangeP[indp,ilat] = np.ma.masked
    if np.ma.is_masked(bowl2[3,ilat]) == False :
        indi = np.ma.nonzero(bowl2[3,ilat]>=density)
        medianToEI[indi,ilat] = np.ma.masked
        rangeToEI[indi,ilat] = np.ma.masked
        noagree_i[indi,ilat] = np.ma.masked
        norangeI[indi,ilat] = np.ma.masked


# -- Create variable bundles
varAtlmedian = {'name': 'Atlantic', 'ToE': medianToEA, 'bowl2': bowl2[1,:], 'bowl1': bowl1[1,:], 'labBowl': labBowl}
varPacmedian = {'name': 'Pacific', 'ToE': medianToEP, 'bowl2': bowl2[2,:], 'bowl1': bowl1[2,:], 'labBowl': labBowl}
varIndmedian = {'name': 'Indian', 'ToE': medianToEI, 'bowl2': bowl2[3,:], 'bowl1': bowl1[3,:], 'labBowl': labBowl}

varAtlrange = {'name': 'Atlantic', 'ToE': rangeToEA, 'bowl2': bowl2[1,:], 'bowl1': bowl1[1,:], 'labBowl': labBowl}
varPacrange = {'name': 'Pacific', 'ToE': rangeToEP, 'bowl2': bowl2[2,:], 'bowl1': bowl1[2,:], 'labBowl': labBowl}
varIndrange = {'name': 'Indian', 'ToE': rangeToEI, 'bowl2': bowl2[3,:], 'bowl1': bowl1[3,:], 'labBowl': labBowl}


# ----- Plot ToE ------
lat2d, density2d = np.meshgrid(lat,density)

# Date
now = datetime.datetime.now()
date = now.strftime("%Y-%m-%d")

if figure == 'median':
    # -- Median

    fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(17,5))

    minmax = [min, finalyear-20 +0.1, deltat]
    unit = 'ToE'
    cmap = cmaps.viridis
    levels = np.arange(minmax[0], minmax[1], minmax[2])

    cnplot = zonal_2D(plt, 'ToE', axes[0, 0], axes[1, 0], 'left', lat, density, varAtlmedian, domrho, cmap, levels)
    cnplot[0].cmap.set_over('0.8')
    cnplot[1].cmap.set_over('0.8')
    axes[0,0].contourf(lat2d, density2d, noagree_a, levels=[0.25,0.5,1.5], colors='None',
                            hatches=['','\\\\'], edgecolor='0.3', linewidth=0.0) # No agreement
    axes[1,0].contourf(lat2d, density2d, noagree_a, levels=[0.25,0.5,1.5], colors='None',
                            hatches=['','\\\\'], edgecolor='0.3', linewidth=0.0)

    # cb = fig.colorbar(cnplot[0],ax=axes.ravel().tolist(), ticks = levels, fraction=0.015, shrink=2.0, pad=0.05)
    cb = fig.colorbar(cnplot[1],ax=axes.ravel().tolist(), ticks = levels, fraction=0.015, shrink=2.0, pad=0.05)
    cb.set_label('%s' % (unit,), fontweight='bold')

    cnplot = zonal_2D(plt, 'ToE', axes[0, 1], axes[1, 1], 'mid', lat, density, varPacmedian, domrho, cmap, levels)
    cnplot[0].cmap.set_over('0.8')
    cnplot[1].cmap.set_over('0.8')
    axes[0,1].contourf(lat2d, density2d, noagree_p, levels=[0.25,0.5,1.5], colors='None',
                            hatches=['','\\\\'], edgecolor='0.3', linewidth=0.0) # No agreement
    axes[1,1].contourf(lat2d, density2d, noagree_p, levels=[0.25,0.5,1.5], colors='None',
                            hatches=['','\\\\'], edgecolor='0.3', linewidth=0.0)
    #axes[0,1].contour(lat2d, density2d, medianToEA, levels=[34.5,35.5], linewidth=2, colors='white')
    #axes[1,1].contour(lat2d, density2d, medianToEA, levels=[34.5,35.5], linewidth=2, colors='white')

    cnplot = zonal_2D(plt, 'ToE', axes[0, 2], axes[1, 2], 'right', lat, density, varIndmedian, domrho, cmap, levels)
    cnplot[0].cmap.set_over('0.8')
    cnplot[1].cmap.set_over('0.8')
    axes[0,2].contourf(lat2d, density2d, noagree_i, levels=[0.25,0.5,1.5], colors='None',
                            hatches=['','\\\\'], edgecolor='0.3', linewidth=0.0) # No agreement
    axes[1,2].contourf(lat2d, density2d, noagree_i, levels=[0.25,0.5,1.5], colors='None',
                            hatches=['','\\\\'], edgecolor='0.3', linewidth=0.0)
    #axes[0,2].contour(lat2d, density2d, medianToEA, levels=[34.5,35.5], linewidth=2, colors='white')
    #axes[1,2].contour(lat2d, density2d, medianToEA, levels=[34.5,35.5], linewidth=2, colors='white')


    plt.subplots_adjust(hspace=.0001, wspace=0.05, left=0.04, right=0.86)

    # Bug with Atlantic colorbar - top panel, so repeat
    cnplot = zonal_2D(plt, 'ToE', axes[0, 0], axes[1, 0], 'left', lat, density, varAtlmedian, domrho, cmap, levels)
    cnplot[0].cmap.set_over('0.8')
    cnplot[1].cmap.set_over('0.8')
    axes[0,0].contourf(lat2d, density2d, noagree_a, levels=[0.25,0.5,1.5], colors='None',
                            hatches=['','\\\\'], edgecolor='0.3', linewidth=0.0) # No agreement
    axes[1,0].contourf(lat2d, density2d, noagree_a, levels=[0.25,0.5,1.5], colors='None',
                            hatches=['','\\\\'], edgecolor='0.3', linewidth=0.0)
    #axes[0,0].contour(lat2d, density2d, medianToEA, levels=[34.5,35.5], linewidth=2, colors='white')
    #axes[1,0].contour(lat2d, density2d, medianToEA, levels=[34.5,35.5], linewidth=2, colors='white')

    if work == 'RCP85':
        if use_piC == False:
            name = 'hist+RCP8.5 vs. histNat'
            if runs_rcp == 'all':
                plotName = 'median_ToE_rcp85vshistNat_'+ str(nb_outliers)+'_outliers_'+str(multstd)+'std'
            else:
                plotName = 'median_ToE_rcp85vshistNat_'+ str(nb_outliers)+'_outliers_'+str(multstd)+'std_samerunsvsPiC'
        else:
            name = 'hist+RCP8.5 vs. PiControl'
            plotName = 'median_ToE_rcp85vspiC_'+ str(nb_outliers)+'_outliers_'+str(multstd)+'std'
    else:
        name = '1pctCO2 vs. PiControl'
        plotName = 'median_ToE_1pctCO2vsPiC_'+ str(nb_outliers_CO2)+'_outliers_'+str(multstd)+'std'
        nruns = nmodels
    plotTitle = 'Multimodel ensemble median ToE for ' + legVar + ', ' + name + ' [> ' + str(multstd) + ' std]' \
        '\n %d models, %d runs '%(nmodels,nruns)


    plt.suptitle(plotTitle, fontweight='bold', fontsize=14, verticalalignment='top')
    plt.figtext(.006,.5,'Density',rotation='vertical',horizontalalignment='left',fontsize=12,fontweight='bold')

    plt.figtext(.5,.01,'Computed by : zonal_toe_median_range.py '+date,fontsize=9,ha='center')
    if use_piC and work =='RCP85':
        plt.figtext(.2,.01,'PiControl : mean(last_240_years)',fontsize=9,ha='center')
    if use_piC == False and work =='RCP85':
        plt.figtext(.2,.01,'Runs : '+runs_rcp,fontsize=9,ha='center')


else:
    # -- 25-75% inter-model range

    fig2, axes = plt.subplots(nrows=2, ncols=3, figsize=(17,5))

    minmax = [0, 91, deltat]
    unit = 'Years'
    cmap = 'YlOrRd' #palettable.cartocolors.sequential.OrYel.mpl_colormap #cmocean.cm.tempo
    levels = np.arange(minmax[0], minmax[1], minmax[2])

    cnplot = zonal_2D(plt, 'ToE', axes[0, 0], axes[1, 0], 'left', lat, density, varAtlrange, domrho, cmap, levels)
    for i in range(2):
        axes[i,0].contourf(lat2d, density2d, norangeA, levels=[0.25,0.5,1.5], colors='0.8') # No emergence
        axes[i,0].contourf(lat2d, density2d, noagree_a, levels=[0.25,0.5,1.5], colors='None',
                                hatches=['','\\\\'], edgecolor='0.3', linewidth=0.0) # No agreement

    cnplot = zonal_2D(plt, 'ToE', axes[0, 1], axes[1, 1], 'mid', lat, density, varPacrange, domrho, cmap, levels)
    for i in range(2):
        axes[i,1].contourf(lat2d, density2d, norangeP, levels=[0.25,0.5,1.5], colors='0.8') # No emergence
        axes[i,1].contourf(lat2d, density2d, noagree_p, levels=[0.25,0.5,1.5], colors='None',
                                hatches=['','\\\\'], edgecolor='0.3', linewidth=0.0) # No agreement

    cnplot = zonal_2D(plt, 'ToE', axes[0, 2], axes[1, 2], 'right', lat, density, varIndrange, domrho, cmap, levels)
    for i in range(2):
        axes[i,2].contourf(lat2d, density2d, norangeI, levels=[0.25,0.5,1.5], colors='0.8') # No emergence
        axes[i,2].contourf(lat2d, density2d, noagree_i, levels=[0.25,0.5,1.5], colors='None',
                                hatches=['','\\\\'], edgecolor='0.3', linewidth=0.0) # No agreement

    plt.subplots_adjust(hspace=.0001, wspace=0.05, left=0.04, right=0.86)

    cb = fig2.colorbar(cnplot[1], ax=axes.ravel().tolist(), ticks = levels, fraction=0.015, shrink=2.0, pad=0.05)
    cb.set_label('%s' % (unit,), fontweight='bold')

    if work == 'RCP85':
        if use_piC == False :
            name = 'hist+RCP8.5 vs. histNat'
            if runs_rcp == 'all':
                plotName = 'range_ToE_rcp85vshistNat_'+ str(nb_outliers)+'_outliers_'+str(multstd)+'std'
            else:
                plotName = 'range_ToE_rcp85vshistNat_'+ str(nb_outliers)+'_outliers_'+str(multstd)+'std_samerunsvsPiC'
        else:
            name = 'hist+RCP8.5 vs. PiControl'
            plotName = 'range_ToE_rcp85vspiC_'+ str(nb_outliers)+'_outliers_'+str(multstd)+'std'

    else:
        name = '1pctCO2 vs. PiControl'
        plotName = 'range_ToE_1pctCO2vsPiC_'+ str(nb_outliers_CO2)+'_outliers_'+str(multstd)+'std'
        nruns = nmodels

    plotTitle = '25-75% multimodel ensemble range of the ToE for ' + legVar + ', ' + name + ' [> ' + str(multstd) + ' std]' \
        '\n %d models, %d runs '%(nmodels,nruns)

    plt.suptitle(plotTitle, fontweight='bold', fontsize=14, verticalalignment='top')
    plt.figtext(.006,.5,'Density',rotation='vertical',horizontalalignment='left',fontsize=12,fontweight='bold')

    plt.figtext(.5,.01,'Computed by : zonal_toe_median_range.py '+date,fontsize=9,ha='center')
    if use_piC and work =='RCP85':
        plt.figtext(.2,.01,'PiControl : mean(last_240_years)',fontsize=9,ha='center')
    if use_piC == False and work == 'RCP85':
        plt.figtext(.2,.01,'Runs : '+runs_rcp,fontsize=9,ha='center')


# -- Show or save figure ---

if work == 'RCP85':
    if use_piC:
        path_end = 'rcp85_PiControl/'
    else:
        path_end = 'rcp85_histNat/'
else:
    path_end = '1pctCO2_piC/'

if outfmt == 'view':
    plt.show()
else:
    plt.savefig('/home/ysilvy/Density_bining/Yona_analysis/figures/models/zonal_ys/ToE/'+path_end+plotName+'.png', bbox_inches='tight')
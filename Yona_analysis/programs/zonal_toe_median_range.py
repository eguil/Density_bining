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

# ----- Workspace ------

indir_toe_rcphn = '/home/ysilvy/Density_bining/Yona_analysis/data/toe_zonal/toe_rcp85_histNat/'
indir_mme_rcp85 = '/data/ericglod/Density_binning/Prod_density_april15/mme_rcp85/'
indir_mme_hn = '/data/ericglod/Density_binning/Prod_density_april15/mme_histNat/'
indir_toe_CO2piC = '/home/ysilvy/Density_bining/Yona_analysis/data/toe_zonal/toe_1pctCO2_piC/'
indir_mme_CO2 = '/data/ericglod/Density_binning/Prod_density_april15/mme_1pctCO2/'
indir_mme_piC = '/data/ericglod/Density_binning/Prod_density_april15/mme_piControl/'

# ----- Work ------

work = 'RCP85'
# work = 'CO2'

# figure = 'median'
figure = 'range'

use_piC = False # Over projection period, signal = RCP-average(histNat), noise = std(histNat)
# use_piC = True # Over projection period, signal = RCP-average(PiControl), noise = std(PiControl)

varname = defVarmme('salinity'); v = 'S'

multStd = 2. # detect ToE at multStd std dev of histnat/PiControl

nb_outliers = 5 # Nb of outlier runs allowed to compute the ToE distribution for hist + RCP8.5


if work == 'RCP85':
    iniyear = 1860
    finalyear = 2100
    models = defModels()
    min = 1960
else:
    iniyear = 0
    finalyear = 140
    models = defModelsCO2piC()
    min = 0
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
    # == Historical vs historicalNat + RCP8.5 vs. historicalNat or RCP8.5 vs. PiControl ==

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
    listfiles = glob.glob(indir_toe_rcphn + '/*.nc')
    nmodels = len(listfiles)

    for i in range(nmodels):

        file_toe = listfiles[i]
        ftoe = open_ncfile(file_toe, 'r')
        name = os.path.basename(file_toe).split('.')[1]
        # Read ToE (members, basin, density, latitude)
        toe2read = ftoe.variables[var + 'ToE2'][:]
        nMembers[i] = toe2read.shape[0]
        print('- Reading ToE of %s with %d members'%(name,nMembers[i]))
        nruns1 = int(nruns + nMembers[i])

        # Save ToE
        varToEA[nruns:nruns1,:,:] = toe2read[:,1,:,:]
        varToEP[nruns:nruns1,:,:] = toe2read[:,2,:,:]
        varToEI[nruns:nruns1,:,:] = toe2read[:,3,:,:]

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

else:
    # == 1%CO2 vs. PiControl ==

    listfiles = glob.glob(indir_toe_CO2piC + '/*.nc')
    nmodels = len(listfiles)

    # -- Initialize varToE containing ToE of all runs
    varToEA = np.ma.masked_all((nmodels, levN, latN))
    varToEP = np.ma.masked_all((nmodels, levN, latN))
    varToEI = np.ma.masked_all((nmodels, levN, latN))

    # Loop over models
    for i in range(nmodels):

        file_toe = listfiles[i]
        ftoe = open_ncfile(file_toe, 'r')
        name = os.path.basename(file_toe).split('.')[1]
        # Read ToE (basin, density, latitude)
        toe2read = ftoe.variables[var + 'ToE2'][:]
        print('- Reading ToE of %s'%(name,))

        # Save ToE
        varToEA[i,:,:] = toe2read[1,:,:]
        varToEP[i,:,:] = toe2read[2,:,:]
        varToEI[i,:,:] = toe2read[3,:,:]

    print('Total number of models:', nmodels)


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
# Done for hist + RCP8.5 only

groups_a = np.ma.masked_all((nruns, levN, latN))
groups_p = np.ma.masked_all((nruns, levN, latN))
groups_i = np.ma.masked_all((nruns, levN, latN))
# Groups values : 0 for no emergence, 1 emergence with positive change, -1 emergence with negative change

# -- Check if emergence

groups_a = np.ma.where(varsignal_a>0,1,-1) # 1 if positive change, -1 if negative change
groups_a[np.ma.where(varToEA > 220)] = 0 # Overwrite with 0 if no emergence
groups_p = np.ma.where(varsignal_p>0,1,-1)
groups_p[np.ma.where(varToEP > 220)] = 0
groups_i = np.ma.where(varsignal_i>0,1,-1)
groups_i[np.ma.where(varToEI > 220)] = 0

nruns_a = np.ma.count(varsignal_a,axis=0) # Nb of runs with unmasked values
nruns_p = np.ma.count(varsignal_p,axis=0)
nruns_i = np.ma.count(varsignal_i,axis=0)

# -- Initialize median and percentiles
medianToEA = np.ma.masked_all((levN, latN))
medianToEP = np.ma.masked_all((levN, latN))
medianToEI = np.ma.masked_all((levN, latN))
percentile16ToEA = np.ma.masked_all((levN, latN))
percentile16ToEP = np.ma.masked_all((levN, latN))
percentile16ToEI = np.ma.masked_all((levN, latN))
percentile84ToEA = np.ma.masked_all((levN, latN))
percentile84ToEP = np.ma.masked_all((levN, latN))
percentile84ToEI = np.ma.masked_all((levN, latN))

# -- Initialize where there is no agreement
noagree_a = np.ma.masked_all((levN, latN))
noagree_p = np.ma.masked_all((levN, latN))
noagree_i = np.ma.masked_all((levN, latN))


# -- Compute distribution if conditions are met

varToEA_clean = np.ma.masked_all((nruns,levN,latN))
varToEP_clean = np.ma.masked_all((nruns,levN,latN))
varToEI_clean = np.ma.masked_all((nruns,levN,latN))

for ilat in range(len(lat)):
    for isig in range(len(density)):

        # Atlantic
        if nruns_a[isig,ilat] >= nruns - 5 : # Criteria for RCP8.5 (lower for 1%CO2)
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
                    (groups_a[:,isig,ilat]==-1).sum()<=-1):
                if (groups_a[:,isig,ilat]==1).sum() > (groups_a[:,isig,ilat]==-1).sum() :
                    varToEA_clean[np.where(groups_a[:,isig,ilat]!=-1),isig,ilat] = varToEA[np.where(groups_a[:,isig,ilat]!=-1),isig,ilat]
                else:
                    varToEA_clean[np.where(groups_a[:,isig,ilat]!=1),isig,ilat] = varToEA[np.where(groups_a[:,isig,ilat]!=1),isig,ilat]
            # No agreement between models
            else:
                noagree_a[isig,ilat] = 1
                varToEA_clean[:,isig,ilat] = varToEA[:,isig,ilat]

        # Pacific
        if nruns_p[isig,ilat] >= nruns - 5 : # Criteria for RCP8.5 (lower for 1%CO2)
            # Positive emergence
            if ((groups_p[:,isig,ilat]==1).sum() > nruns_p[isig,ilat]/2) and ((groups_p[:,isig,ilat]==-1).sum() < nb_outliers):
                varToEP_clean[np.where(groups_p[:,isig,ilat]!=-1),isig,ilat] = varToEP[np.where(groups_p[:,isig,ilat]!=-1),isig,ilat]
            # Negative emergence
            elif ((groups_p[:,isig,ilat]==-1).sum() > nruns_p[isig,ilat]/2) and ((groups_p[:,isig,ilat]==1).sum() < nb_outliers):
                 varToEP_clean[np.where(groups_p[:,isig,ilat]!=1),isig,ilat] = varToEP[np.where(groups_p[:,isig,ilat]!=1),isig,ilat]
            # No emergence
            elif ((groups_p[:,isig,ilat]==0).sum() > nruns_p[isig,ilat]/2) and ((groups_p[:,isig,ilat]==1).sum()<=nb_outliers or
                    (groups_p[:,isig,ilat]==-1).sum()<=-1):
                if (groups_p[:,isig,ilat]==1).sum() > (groups_p[:,isig,ilat]==-1).sum() :
                    varToEP_clean[np.where(groups_p[:,isig,ilat]!=-1),isig,ilat] = varToEP[np.where(groups_p[:,isig,ilat]!=-1),isig,ilat]
                else:
                    varToEP_clean[np.where(groups_p[:,isig,ilat]!=1),isig,ilat] = varToEP[np.where(groups_p[:,isig,ilat]!=1),isig,ilat]
            # No agreement between models
            else:
                noagree_p[isig,ilat] = 1
                varToEP_clean[:,isig,ilat] = varToEP[:,isig,ilat]

        # Indian
        if nruns_i[isig,ilat] >= nruns - 5 : # Criteria for RCP8.5 (lower for 1%CO2)
            # Positive emergence
            if ((groups_i[:,isig,ilat]==1).sum() > nruns_i[isig,ilat]/2) and ((groups_i[:,isig,ilat]==-1).sum() < nb_outliers):
                varToEI_clean[np.where(groups_i[:,isig,ilat]!=-1),isig,ilat] = varToEI[np.where(groups_i[:,isig,ilat]!=-1),isig,ilat]
            # Negative emergence
            elif ((groups_i[:,isig,ilat]==-1).sum() > nruns_i[isig,ilat]/2) and ((groups_i[:,isig,ilat]==1).sum() < nb_outliers):
                 varToEI_clean[np.where(groups_i[:,isig,ilat]!=1),isig,ilat] = varToEI[np.where(groups_i[:,isig,ilat]!=1),isig,ilat]
            # No emergence
            elif ((groups_i[:,isig,ilat]==0).sum() > nruns_i[isig,ilat]/2) and ((groups_i[:,isig,ilat]==1).sum()<=nb_outliers or
                    (groups_i[:,isig,ilat]==-1).sum()<=-1):
                if (groups_i[:,isig,ilat]==1).sum() > (groups_i[:,isig,ilat]==-1).sum() :
                    varToEI_clean[np.where(groups_i[:,isig,ilat]!=-1),isig,ilat] = varToEI[np.where(groups_i[:,isig,ilat]!=-1),isig,ilat]
                else:
                    varToEI_clean[np.where(groups_i[:,isig,ilat]!=1),isig,ilat] = varToEI[np.where(groups_i[:,isig,ilat]!=1),isig,ilat]
            # No agreement between models
            else:
                noagree_i[isig,ilat] = 1
                varToEI_clean[:,isig,ilat] = varToEI[:,isig,ilat]

print('Loops done')
medianToEA = np.ma.around(np.ma.median(varToEA_clean, axis=0)) + iniyear
percentile16ToEA = np.ma.around(np.percentile(varToEA_clean, 16, axis=0)) + iniyear
percentile84ToEA = np.ma.around(np.percentile(varToEA_clean, 84, axis=0)) + iniyear
medianToEP = np.ma.around(np.ma.median(varToEP_clean, axis=0)) + iniyear
percentile16ToEP = np.ma.around(np.percentile(varToEP_clean, 16, axis=0)) + iniyear
percentile84ToEP = np.ma.around(np.percentile(varToEP_clean, 84, axis=0)) + iniyear
medianToEI = np.ma.around(np.ma.median(varToEI_clean, axis=0)) + iniyear
percentile16ToEI = np.ma.around(np.percentile(varToEI_clean, 16, axis=0)) + iniyear
percentile84ToEI = np.ma.around(np.percentile(varToEI_clean, 84, axis=0)) + iniyear


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
rangeToEA = percentile84ToEA - percentile16ToEA
rangeToEP = percentile84ToEP - percentile16ToEP
rangeToEI = percentile84ToEI - percentile16ToEI

rangeToEA[percentile84ToEA>finalyear-20] = np.ma.masked
norangeA = np.where(percentile84ToEA>finalyear-20,1,0)
rangeToEP[percentile84ToEP>finalyear-20] = np.ma.masked
norangeP = np.where(percentile84ToEP>finalyear-20,1,0)
rangeToEI[percentile84ToEI>finalyear-20] = np.ma.masked
norangeI = np.where(percentile84ToEI>finalyear-20,1,0)

# Mask
rangeToEA[rangeToEA==0] = np.ma.masked
rangeToEP[rangeToEP==0] = np.ma.masked
rangeToEI[rangeToEI==0] = np.ma.masked

# ----- Read bowl position and mask points above ------

if work == 'RCP85':
    # Read files
    file_rcp85 = 'cmip5.multimodel_Nat.rcp85.ensm.an.ocn.Omon.density_zon1D.nc'
    file_hn = 'cmip5.multimodel_Nat.historicalNat.ensm.an.ocn.Omon.density_zon1D.nc'
    f2 = open_ncfile(indir_mme_rcp85+file_rcp85,'r')
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

if figure == 'median':
    # -- Median

    fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(17,5))

    minmax = [min, finalyear-20 +0.1, deltat]
    unit = 'ToE'
    cmap = 'jet_r'
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

    cnplot = zonal_2D(plt, 'ToE', axes[0, 2], axes[1, 2], 'right', lat, density, varIndmedian, domrho, cmap, levels)
    cnplot[0].cmap.set_over('0.8')
    cnplot[1].cmap.set_over('0.8')
    axes[0,2].contourf(lat2d, density2d, noagree_i, levels=[0.25,0.5,1.5], colors='None',
                            hatches=['','\\\\'], edgecolor='0.3', linewidth=0.0) # No agreement
    axes[1,2].contourf(lat2d, density2d, noagree_i, levels=[0.25,0.5,1.5], colors='None',
                            hatches=['','\\\\'], edgecolor='0.3', linewidth=0.0)

    plt.subplots_adjust(hspace=.0001, wspace=0.05, left=0.04, right=0.86)


    if work == 'RCP85':
        name = 'hist+RCP8.5 vs. histNat'
    else:
        name = '1pctCO2 vs. PiControl'
        nruns = nmodels
    plotTitle = 'Multimodel ensemble median ToE for ' + legVar + ', ' + name + ' [> ' + str(multStd) + ' std]' \
        '\n %d models, %d runs '%(nmodels,nruns)


    plt.suptitle(plotTitle, fontweight='bold', fontsize=14, verticalalignment='top')

    plt.figtext(.5,.02,'Computed by : zonal_toe_median_range.py',fontsize=9,ha='center')


else:
    # -- 16-84% inter-model range

    fig2, axes = plt.subplots(nrows=2, ncols=3, figsize=(17,5))

    minmax = [0, 101, deltat]
    unit = 'Years'
    cmap = 'jet_r'
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
        name = 'hist+RCP8.5 vs. histNat'
    else:
        name = '1pctCO2 vs. PiControl'
        nruns = nmodels

    plotTitle = '16-84% multimodel ensemble range of the ToE for ' + legVar + ', ' + name + ' [> ' + str(multStd) + ' std]' \
        '\n %d models, %d runs '%(nmodels,nruns)

    plt.suptitle(plotTitle, fontweight='bold', fontsize=14, verticalalignment='top')

    plt.figtext(.5,.02,'Computed by : zonal_toe_median_range.py',fontsize=9,ha='center')


plt.show()
#!/bin/env python
# -*- coding: utf-8 -*-

"""
Python matplotlib
Plot pseudo-z/latitude Time of Emergence (median and distribution range)

"""

import os
import glob
from netCDF4 import Dataset as open_ncfile
import matplotlib.pyplot as plt
from maps_matplot_lib import defVarmme, zon_2Dz
from modelsDef import defModels
from lib_remapping import remaptoz
import numpy as np
import colormaps as cmaps
import cmocean
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

# ------------------------------------
# ----- Read ToE for each model ------
# ------------------------------------

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

# ----------------------------------------------------------------
# ----- Compute distribution according to Lyu et. al method ------
# ----------------------------------------------------------------

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

# ------------------------------        
# -------- Remapping -----------
# ------------------------------

# ---- Regroup basins before remapping ----
medianToEr = np.ma.masked_all((basinN,levN,latN))
medianToEr[1,:,:] = medianToEA
medianToEr[2,:,:] = medianToEP
medianToEr[3,:,:] = medianToEI

rangeToEr = np.ma.masked_all((basinN,levN,latN))
rangeToEr[1,:,:] = rangeToEA
rangeToEr[2,:,:] = rangeToEP
rangeToEr[3,:,:] = rangeToEI

noagreer = np.ma.masked_all((basinN,levN,latN))
noagreer[1,:,:] = noagree_a
noagreer[2,:,:] = noagree_p
noagreer[3,:,:] = noagree_i

noranger = np.ma.masked_all((basinN,levN,latN))
noranger[1,:,:] = norangeA
noranger[2,:,:] = norangeP
noranger[3,:,:] = norangeI

# --- Read reference pseudo-depth used for remapping ---
indir_z = '/home/ysilvy/Density_bining/Yona_analysis/data/remaptoz/'
file_z = ''
fz = open_ncfile(indir_z+file_z,'r')
pseudo_depth = fz.variables['pseudo_depth'][:]

# -- Target grid for remapping - WOA13 grid --
targetz = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80,
   85, 90, 95, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 375,
   400, 425, 450, 475, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950,
   1000, 1050, 1100, 1150, 1200, 1250, 1300, 1350, 1400, 1450, 1500, 1550,
   1600, 1650, 1700, 1750, 1800, 1850, 1900, 1950, 2000, 2100, 2200, 2300,
   2400, 2500, 2600, 2700, 2800, 2900, 3000, 3100, 3200, 3300, 3400, 3500,
   3600, 3700, 3800, 3900, 4000, 4100, 4200, 4300, 4400, 4500, 4600, 4700,
   4800, 4900, 5000, 5100, 5200, 5300, 5400, 5500]

# -- Remap --
medianToEz = np.ma.around(remaptoz(medianToEr,pseudo_depth,targetz))
print('Median remapping done')
rangeToEz = np.ma.around(remaptoz(rangeToEr,pseudo_depth,targetz))
print('Range remapping done')
noagreez = remaptoz(noagreer,pseudo_depth,targetz)
print('noagree remapping done')
norangez = remaptoz(noranger,pseudo_depth,targetz)
print('noagree remapping done')

# # Mask
# #rangeToEz[medianToEz > finalyear-20] = np.ma.masked # Mask points where median hasn't emerged
# #medianToEz[noagreez==1] = np.ma.masked
# medianToEz[medianToEz<iniyear] = np.ma.masked
# rangeToEz[rangeToEz==0] = np.ma.masked

# # ----- Read depth of bowl ------
# # Read files
# file1d_rcp85 = 'cmip5.multimodel_Nat.rcp85.ensm.an.ocn.Omon.density_zon1D.nc'
# f1drcp85 = open_ncfile(indir_mme_rcp85+file1d_rcp85,'r')
# file1d_hn = 'cmip5.multimodel_Nat.historicalNat.ensm.an.ocn.Omon.density_zon1D.nc'
# fhn = open_ncfile(indir_mme_hn+file1d_hn,'r')

# # Read variables
# bowl2z = f1drcp85.variables['ptopdepth'][-5:,:,:]
# bowl1z = fhn.variables['ptopdepth'][-5:,:,:]

# # Average
# bowl2z = np.ma.average(bowl2z, axis=0)
# bowl1z = np.ma.average(bowl1z, axis=0)

# labBowl = ['histNat', 'RCP8.5']
labBowl=[None,None]

# -- Make variable bundles for each basin
varAtlmedian = {'name': 'Atlantic', 'var_change': medianToEz[1,:,:], 'bowl1': bowl1z[1,:], 'bowl2': bowl2z[1,:], 'labBowl': labBowl}
varPacmedian = {'name': 'Pacific', 'var_change': medianToEz[2,:,:], 'bowl1': bowl1z[2,:], 'bowl2': bowl2z[2,:], 'labBowl': labBowl}
varIndmedian = {'name': 'Indian', 'var_change': medianToEz[3,:,:], 'bowl1': bowl1z[3,:], 'bowl2': bowl2z[3,:], 'labBowl': labBowl}

varAtlrange = {'name': 'Atlantic', 'var_change': rangeToEz[1,:,:], 'bowl1': bowl1z[1,:], 'bowl2': bowl2z[1,:], 'labBowl': labBowl}
varPacrange = {'name': 'Pacific', 'var_change': rangeToEz[2,:,:], 'bowl1': bowl1z[2,:], 'bowl2': bowl2z[2,:], 'labBowl': labBowl}
varIndrange = {'name': 'Indian', 'var_change': rangeToEz[3,:,:], 'bowl1': bowl1z[3,:], 'bowl2': bowl2z[3,:], 'labBowl': labBowl}


# ----- Plot -----

domzed = [0,500,5000]
# Create meshgrid
lat2d, lev2d = np.meshgrid(lat, targetz)

# Date
now = datetime.datetime.now()
date = now.strftime("%Y-%m-%d")

if figure == 'median':
    # -- Median

    # -- Create figure and axes instances
    fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(17, 5))

    # -- Color map
    cmap = cmaps.plasma
    # -- Unit
    unit = 'ToE'
    # -- Levels
    minmax = [min, finalyear-20 +0.1, deltat]
    levels = np.arange(minmax[0], minmax[1], minmax[2]) 
    ext_cmap = 'both'
    # -- Put everything into a dictionary
    contourDict = {'cmap':cmap, 'levels':levels, 'ext_cmap':ext_cmap}

    # -- Contourf
    # Atlantic
    cnplot = zon_2Dz(plt, axes[0,0], axes[1,0], 'left', lat, targetz, varAtlmedian,
                     contourDict, domzed)
    cnplot[0].cmap.set_over('0.8')
    cnplot[1].cmap.set_over('0.8')
    axes[0,0].contourf(lat2d, lev2d, noagreez[1,:,:], levels=[0.25,0.5,1.5], colors='None', hatches=['','\\\\'], edgecolor='0.3', linewidth=0.0) # No agreement
    axes[1,0].contourf(lat2d, lev2d, noagreez[1,:,:], levels=[0.25,0.5,1.5], colors='None',hatches=['','\\\\'], edgecolor='0.3', linewidth=0.0)

    # -- Add colorbar
    cb = fig.colorbar(cnplot[1], ax=axes.ravel().tolist(), ticks=levels, fraction=0.015, shrink=2.0, pad=0.05)
    cb.set_label('%s' % (unit,), fontweight='bold')

    # Pacific
    cnplot = zon_2Dz(plt, axes[0,1], axes[1,1], 'mid', lat, targetz, varPacmedian,
                     contourDict, domzed)
    cnplot[0].cmap.set_over('0.8')
    cnplot[1].cmap.set_over('0.8')
    axes[0,1].contourf(lat2d, lev2d, noagreez[2,:,:], levels=[0.25,0.5,1.5], colors='None',
                                hatches=['','\\\\'], edgecolor='0.3', linewidth=0.0) # No agreement
    axes[1,1].contourf(lat2d, lev2d, noagreez[2,:,:], levels=[0.25,0.5,1.5], colors='None',
                                hatches=['','\\\\'], edgecolor='0.3', linewidth=0.0)

    # Indian
    cnplot = zon_2Dz(plt, axes[0,2], axes[1,2], 'right', lat, targetz, varIndmedian,
                     contourDict, domzed)
    cnplot[0].cmap.set_over('0.8')
    cnplot[1].cmap.set_over('0.8')
    axes[0,2].contourf(lat2d, lev2d, noagreez[3,:,:], levels=[0.25,0.5,1.5], colors='None',
                                hatches=['','\\\\'], edgecolor='0.3', linewidth=0.0) # No agreement
    axes[1,2].contourf(lat2d, lev2d, noagreez[3,:,:], levels=[0.25,0.5,1.5], colors='None',
                                hatches=['','\\\\'], edgecolor='0.3', linewidth=0.0)

    plt.subplots_adjust(hspace=.0001, wspace=0.05, left=0.04, right=0.86)


    if work == 'RCP85':
        if use_piC == False:
            name = 'hist+RCP8.5 vs. histNat'
            if runs_rcp == 'all':
                plotName = 'remapping_median_ToE_rcp85vshistNat_'+ str(nb_outliers)+'_outliers_'+str(multstd)+'std'
            else:
                plotName = 'remapping_median_ToE_rcp85vshistNat_'+ str(nb_outliers)+'_outliers_'+str(multstd)+'std_samerunsvsPiC'
        else:
            name = 'hist+RCP8.5 vs. PiControl'
            plotName = 'remapping_median_ToE_rcp85vspiC_'+ str(nb_outliers)+'_outliers_'+str(multstd)+'std'
    else:
        name = '1pctCO2 vs. PiControl'
        plotName = 'remapping_median_ToE_1pctCO2vsPiC_'+ str(nb_outliers_CO2)+'_outliers_'+str(multstd)+'std'
        nruns = nmodels
    
    # -- Add title    
    plotTitle = 'Multimodel ensemble median ToE for ' + legVar + ', ' + name + ' [> ' + str(multstd) + ' std]' \
        '\n %d models, %d runs '%(nmodels,nruns)
    
    plt.suptitle(plotTitle, fontweight='bold', fontsize=14, verticalalignment='top')

    plt.figtext(.5,.02,'Computed by : remaptoz2_toe_median_range.py '+date,fontsize=9,ha='center')
    plt.figtext(.004,.65,'Pseudo-depth (m)',rotation='vertical',fontweight='bold')

    figureDir = 'models/zonal_remaptoz/'
    if fig == 'save':
        plt.savefig('/home/ysilvy/Density_bining/Yona_analysis/figures/'+figureDir+plotName+'.png', bbox_inches='tight')

else:
    # -- 25-75% range

    # -- Create figure and axes instances
    fig2, axes = plt.subplots(nrows=2, ncols=3, figsize=(17, 5))

    # -- Color map
    cmap = cmocean.cm.tempo
    # -- Unit
    unit = 'Years'
    # -- Levels
    minmax = [0, 101, deltat]
    levels = np.arange(minmax[0], minmax[1], minmax[2])

    ext_cmap = 'max'
    # -- Put everything into a dictionary
    contourDict = {'cmap':cmap, 'levels':levels, 'ext_cmap':ext_cmap}

    # -- Contourf
    cnplot2 = zon_2Dz(fig2, axes[0,0], axes[1,0], 'left', lat, targetz, varAtlrange,
                     contourDict, domzed)
    cnplot2 = zon_2Dz(fig2, axes[0,1], axes[1,1], 'mid', lat, targetz, varPacrange,
                     contourDict, domzed)
    cnplot2 = zon_2Dz(fig2, axes[0,2], axes[1,2], 'right', lat, targetz, varIndrange,
                     contourDict, domzed)

    for i in range(2):
        for ibasin in range(3):
            axes[i,ibasin].contourf(lat2d, lev2d, norangez[ibasin+1,:,:], levels=[0.25,0.5,1.5], colors='0.8') # No emergence
            axes[i,ibasin].contourf(lat2d, lev2d, noagreez[ibasin+1,:,:], levels=[0.25,0.5,1.5], colors='None', hatches=['','\\\\'], edgecolor='0.3', linewidth=0.0) # No agreement

    fig2.subplots_adjust(hspace=.0001, wspace=0.05, left=0.04, right=0.86)

    # -- Add colorbar
    cb = fig2.colorbar(cnplot2[1], ax=axes.ravel().tolist(), ticks=levels, fraction=0.015, shrink=2.0, pad=0.05)
    cb.set_label('%s' % (unit,), fontweight='bold')

    if work == 'RCP85':
        if use_piC == False :
            name = 'hist+RCP8.5 vs. histNat'
            if runs_rcp == 'all':
                plotName = 'remapping_range_ToE_rcp85vshistNat_'+ str(nb_outliers)+'_outliers_'+str(multstd)+'std'
            else:
                plotName = 'remapping_range_ToE_rcp85vshistNat_'+ str(nb_outliers)+'_outliers_'+str(multstd)+'std_samerunsvsPiC'
        else:
            name = 'hist+RCP8.5 vs. PiControl'
            plotName = 'remapping_range_ToE_rcp85vspiC_'+ str(nb_outliers)+'_outliers_'+str(multstd)+'std'

    else:
        name = '1pctCO2 vs. PiControl'
        plotName = 'remapping_range_ToE_1pctCO2vsPiC_'+ str(nb_outliers_CO2)+'_outliers_'+str(multstd)+'std'
        nruns = nmodels
    
    # -- Add title
    plotTitle = '25-75% multimodel ensemble range of the ToE for ' + legVar + ', ' + name + ' [> ' + str(multstd) + ' std]' \
        '\n %d models, %d runs '%(nmodels,nruns)


    plt.suptitle(plotTitle, fontweight='bold', fontsize=14, verticalalignment='top')

    plt.figtext(.5,.01,'Computed by : remaptoz2_toe_median_range.py '+date,fontsize=9,ha='center')
    plt.figtext(.004,.65,'Pseudo-depth (m)',rotation='vertical',fontweight='bold')

    if use_piC and work =='RCP85':
        plt.figtext(.2,.01,'PiControl : mean(last_240_years)',fontsize=9,ha='center')
    if use_piC == False and work == 'RCP85':
        plt.figtext(.2,.01,'Runs : '+runs_rcp,fontsize=9,ha='center')

figureDir = 'models/zonal_remaptoz/'

if fig == 'save':
    plt.savefig('/home/ysilvy/Density_bining/Yona_analysis/figures/'+figureDir+plotName+'.png', bbox_inches='tight')
else:
    plt.show()

#!/bin/env python
# -*- coding: utf-8 -*-

"""
Python matplotlib
Make scatter plots of ToE/noise hist+RCP85 vs. histNat or vs. PiControl in all 5 domains
"""

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset as open_ncfile
from maps_matplot_lib import defVarmme
from modelsDef import defModels, defModelsCO2piC
from matplotlib.ticker import AutoMinorLocator, MultipleLocator
import glob
import os

# ----- Workspace -----

# Directory
indir_rcphn = '/home/ysilvy/Density_bining/Yona_analysis/data/toe_rcp85_histNat_average_signal/'
indir_rcppiC = '/home/ysilvy/Density_bining/Yona_analysis/data/toe_rcp85_PiControl_average_signal/'
indir_noise = '/home/ysilvy/Density_bining/Yona_analysis/data/noise_estimate/'
indir_histNat = '/data/ericglod/Density_binning/Prod_density_april15/mme_histNat/'
indir_piC = '/data/ericglod/Density_binning/Prod_density_april15/mme_piControl/'


# ----- Work ------

varname = defVarmme('salinity'); v = 'S'

method = 'average_signal' # Average signal and noise in the box, then compute ToE

# -- Choose which 'noise' to use
# method_noise = 'average_std' # Average the standard deviation of histNat or PiC in the specified domains
method_noise = 'average_histNat' # Average histNat or PiC in the specified domains then determine the std
# of this averaged value

models = defModels()

domains = ['Southern ST', 'SO', 'Northern ST', 'North Atlantic', 'North Pacific']

markers = ['o','v','^','8','s','p','h','D','d','*','>']

var = varname['var_zonal_w/bowl']

# ----- Read ToE and noise for each run ------

# ToEhn, ToEpiC, noisehn, noisepiC
listfiles_rcphn = glob.glob(indir_rcphn + method_noise + '/*.nc')
nmodels = len(listfiles_rcphn)
model_names = np.empty(nmodels).astype('S20')

for i in range(nmodels):

    # Read ToE RCP8.5 vs. histNat (members, basin, domain)
    filetoehn = listfiles_rcphn[i]
    ftoehn = open_ncfile(filetoehn)
    name = os.path.basename(filetoehn).split('.')[1]
    model_names[i] = name
    print(name)
    toehnread = ftoehn.variables[var + 'ToE2'][:]

    # Read noise histNat
    if method_noise == 'average_histNat':
        # Here we read the std of the averaged histNat for all runs, then take the max as our noise
        filenoise = 'cmip5.' + name + '.noise_domains_hist_histNat.std_of_average.nc'
        fnoise = open_ncfile(indir_noise + filenoise,'r')
        varstdhn = fnoise.variables[var+'stdhn'][:] # Already averaged in the domains (members,basin,domain)
        varnoisehn = np.ma.max(varstdhn,axis=0)
    else:
        # Read histNat ensemble mean
        filehn = glob.glob(indir_histNat + 'cmip5.' + name + '.'+'*.zon2D.nc')[0]
        fhn = open_ncfile(filehn,'r')
        # Read std of histNat (max std of all runs for each model)
        varstdhn = fhn.variables[var+'Std'][:] # Not averaged yet

    # Read corresponding model for RCP8.5 vs. piC
    filetoepiC = glob.glob(indir_rcppiC + method_noise + '/cmip5.'+name+'.'+'*.nc')
    if len(filetoepiC) != 0 :
        ftoepiC = open_ncfile(filetoepiC[0])
        toepiCread = ftoepiC.variables[var + 'ToE2'][:]

        # Read noise PiControl
        if method_noise == 'average_histNat':
            varnoisepiC = fnoise.variables[var+'stdpiC'][:] # Already averaged in the domains (basin,domain)
        else:
            filepiC = glob.glob(indir_piC + 'cmip5.' + name + '.'+'*.zon2D.nc')[0]
            fpiC = open_ncfile(filepiC,'r')
            # Read std of PiControl
            varstdpiC = fpiC.variables[var+'Std'][:] # Not averaged yet









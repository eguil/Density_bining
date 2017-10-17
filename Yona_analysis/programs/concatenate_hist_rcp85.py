#!/bin/env python
# -*- coding: utf-8 -*-

""" Concatenate historical runs with RCP85 runs into one file for each run of each model"""

import numpy as np
from netCDF4 import Dataset as open_ncfile
from modelsDef import defModels
import glob
import os
# from nco import Nco
# nco = Nco()

indir_hist = '/data/ericglod/Density_binning/Prod_density_april15/historical/'
indir_rcp85 = '/home/ericglod/data/Density_binning/Prod_density_april15/rcp85/'
outdir = '/home/ysilvy/Density_bining/Yona_analysis/data/hist_rcp85/'

models = defModels()

for i, model in enumerate(models):

    print('- Working on ' + model['name'])

    # Read historical files
    listruns2D = glob.glob(indir_hist + 'cmip5.' + model['name'] +'.'+ '*zon2D.nc')
    listruns1D = glob.glob(indir_hist + 'cmip5.' + model['name'] +'.'+ '*zon1D.nc')
    print('    Historical runs :\n'+'\n'.join(listruns2D))
    nruns = len(listruns2D)

    listruns2D_rcp85 = glob.glob(indir_rcp85 + 'cmip5.' + model['name'] +'.'+ '*zon2D.nc')
    listruns1D_rcp85 = glob.glob(indir_rcp85 + 'cmip5.' + model['name'] +'.'+ '*zon1D.nc')
    print('    RCP85 runs:\n'+'\n'.join(listruns2D_rcp85))

    for j in range(nruns):
        fname = os.path.basename(listruns2D[j])
        nb = fname.split('.')[3] # Run number, e.g. r1i1p1
        print('    '+nb)

        for k in range(len(listruns2D_rcp85)):
            fname2 = os.path.basename(listruns2D_rcp85[k])
            nb2 = fname2.split('.')[3]
            if nb == nb2: # Check if run numbers match
                print('    '+nb2)
                # Create new file names, respecting the name convention
                fragments = fname2.split('.')
                fragments = np.delete(fragments,2)
                fragments = np.insert(fragments,2,'historical-rcp85')
                s = '.'
                outfile2D = s.join(fragments)
                print('        '+outfile2D)
                fragments = np.delete(fragments,-2)
                fragments = np.insert(fragments, -1, os.path.basename(listruns1D_rcp85[k]).split('.')[-2])
                outfile1D = s.join(fragments)
                print('        '+outfile1D)
    print('\n')

                # Concatenate
                #nco.ncrcat(listruns2D[j],listruns2D_rcp85[k],outdir+outfile1D)
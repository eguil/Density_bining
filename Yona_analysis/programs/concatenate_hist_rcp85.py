#!/bin/env python
# -*- coding: utf-8 -*-

""" Concatenate historical runs with RCP85 runs into one file for each run of each model,
or duplicate historialNat ensemble mean for each model"""

import numpy as np
from netCDF4 import Dataset as open_ncfile
from modelsDef import defModels
import glob
import os
from subprocess import call

# ----- Work -----

work = 'hist-rcp85' # Concatenate historical runs with rcp8.5 runs for each model
# work = 'histNat' # Duplicate historicalNat ensemble mean for each model

# ----------------

# == Historical + rcp8.5 ==
if work == 'hist-rcp85':

    indir_hist = '/data/ericglod/Density_binning/Prod_density_april15/historical/'
    indir_rcp85 = '/home/ericglod/data/Density_binning/Prod_density_april15/rcp85/'
    outdir = '/home/ysilvy/Density_bining/Yona_analysis/data/hist_rcp85/'

    models = defModels()

    for i, model in enumerate(models):

        print('- Working on ' + model['name'])

        # Read historical files
        listruns2D = glob.glob(indir_hist + 'cmip5.' + model['name'] + '.' + '*zon2D.nc')
        listruns1D = glob.glob(indir_hist + 'cmip5.' + model['name'] + '.' + '*zon1D.nc')
        # print('    Historical runs :\n'+'\n'.join(listruns2D))
        nruns = len(listruns2D)

        listruns2D_rcp85 = glob.glob(indir_rcp85 + 'cmip5.' + model['name'] + '.' + '*zon2D.nc')
        listruns1D_rcp85 = glob.glob(indir_rcp85 + 'cmip5.' + model['name'] + '.' + '*zon1D.nc')
        # print('    RCP85 runs:\n'+'\n'.join(listruns2D_rcp85))

        for j in range(nruns):
            fname = os.path.basename(listruns2D[j])
            nb = fname.split('.')[3]  # Run number, e.g. r1i1p1
            # print('    '+nb)

            for k in range(len(listruns2D_rcp85)):
                fname2 = os.path.basename(listruns2D_rcp85[k])
                nb2 = fname2.split('.')[3]
                if nb == nb2:  # Check if run numbers match
                    # print('    '+nb2)
                    # Create new file names, respecting the name convention
                    fragments = fname2.split('.')
                    fragments = (np.delete(fragments,2)).astype('S30') # Extend string length otherwise historical-rcp85
                    # gets cut
                    fragments = np.insert(fragments,2,'historical-rcp85')
                    s = '.'
                    outfile2D = outdir + s.join(fragments)
                    # print('        '+outfile2D)
                    file_end = fragments[-2].split('_')
                    fragments[-2] = file_end[0]+'_zon1D'
                    outfile1D = outdir + s.join(fragments)
                    # print('        '+outfile1D)
        # print('\n')
                    # Rename for concatenation
                    infile2D_hist = listruns2D[j]
                    infile2D_rcp85 = listruns2D_rcp85[k]
                    infile1D_hist = listruns1D[j]
                    infile1D_rcp85 = listruns1D_rcp85[k]
                    # Concatenate
                    call('ncrcat {0} {1} {2}'.format(infile2D_hist,infile2D_rcp85,outfile2D), shell=True)
                    call('ncrcat {0} {1} {2}'.format(infile1D_hist,infile1D_rcp85,outfile1D), shell=True)
                    print(outfile2D)
                    print(outfile1D)


# == HistoricalNat ==
if work == 'histNat':

    indir_histNat = '/data/ericglod/Density_binning/Prod_density_april15/mme_histNat/'
    outdir = '/home/ysilvy/Density_bining/Yona_analysis/data/histNat_histNat/'

    models = defModels()

    for i, model in enumerate(models):

        print('- Working on ' + model['name'])

        # Read historicalNat files
        filehn_2D = 'cmip5.' + model['name'] + '.historicalNat.ensm.an.ocn.Omon.density.ver-' + \
                    model['file_end_histNat'] + '_zon2D.nc'
        filehn_1D = 'cmip5.' + model['name'] + '.historicalNat.ensm.an.ocn.Omon.density.ver-' + \
                    model['file_end_histNat'] + '_zon1D.nc'

        fname = os.path.basename(filehn_2D)
        fragments = fname.split('.')
        fragments = (np.delete(fragments,2)).astype('S30')
        fragments = np.insert(fragments,2,'historicalNat-historicalNat')
        s = '.'
        outfile2D = outdir + s.join(fragments)
        print(' '+os.path.basename(outfile2D))
        # Now file 1D
        fragments = np.delete(fragments,-2)
        fragments = np.insert(fragments, -1, os.path.basename(filehn_1D).split('.')[-2])
        outfile1D = outdir + s.join(fragments)
        print(' '+os.path.basename(outfile1D))

        # Concatenate
        call('ncrcat {0} {1} {2}'.format(indir_histNat + filehn_2D, indir_histNat + filehn_2D,outfile2D), shell=True)
        call('ncrcat {0} {1} {2}'.format(indir_histNat + filehn_1D,indir_histNat + filehn_1D,outfile1D), shell=True)

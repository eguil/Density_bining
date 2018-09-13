#!/bin/env python
# -*- coding: utf-8 -*-
"""
Python matplotlib
Make surface transformation plots

(c) Eric Guilyardi Sept. 2018

"""
from netCDF4 import Dataset as open_ncfile
import matplotlib.pyplot as plt
import numpy as npy

# -------------------------------------------------------------------------------
#                                Define work
# -------------------------------------------------------------------------------

#inDir = '/Users/ericg/Projets/Density_bining/Nudge_CMIP6_IPSL/'
inDir = '/data/ericglod/Runs_nudges_sept2018/'

run1 = 'historical_r1i1p1f1'
file1  = 'IPSL-CM6A-LR_historical_r1i1p1f1_gn_1950_2009_seasmean_transf_north40.nc'
run2 = 'CM6-pace-TSTr7fgT'
file2  = run2+'_1950_2009_seasmean_transf_north40.nc'
run3 = 'CM6-pace-TSTr7vg'
file3  = run3+'_1950_2009_seasmean_transf_north40.nc'

#
# -- Open netcdf files
nc1 = open_ncfile(inDir + run1 + '/' + file1)
nc2 = open_ncfile(inDir + run2 + '/' + file2)
nc3 = open_ncfile(inDir + run3 + '/' + file3)

# -- Read variables for North Atl (anual mean from cdo yearmean on monthly data)

trfatltot1 = nc1.variables['trsftotAtl'][0,:].squeeze()
trfatlhef1 = nc1.variables['trsfhefAtl'][0,:].squeeze()
trfatlwfo1 = nc1.variables['trsfwfoAtl'][0,:].squeeze()

trfatltot2 = nc2.variables['trsftotAtl'][0,:].squeeze()
trfatlhef2 = nc2.variables['trsfhefAtl'][0,:].squeeze()
trfatlwfo2 = nc2.variables['trsfwfoAtl'][0,:].squeeze()

trfatltot3 = nc3.variables['trsftotAtl'][0,:].squeeze()
trfatlhef3 = nc3.variables['trsfhefAtl'][0,:].squeeze()
trfatlwfo3 = nc3.variables['trsfwfoAtl'][0,:].squeeze()


# axis
levr = nc1.variables['rhon'][:]
sigmin = 22
sigmax = 29
trfmin = -10
trfmax = 40

plt.axis([sigmin, sigmax, trfmin, trfmax])

plt.plot(levr, trfatltot1, c = 'black', label = run1)
plt.plot(levr, trfatlhef1, c = 'black', linestyle ='--')
plt.plot(levr, trfatlwfo1, c = 'black', linestyle =':')

plt.plot(levr, trfatltot2, c = 'b', label = run2)
plt.plot(levr, trfatlhef2, c = 'b', linestyle ='--')
plt.plot(levr, trfatlwfo2, c = 'b', linestyle =':')

plt.plot(levr, trfatltot3, c = 'r', label = run3)
plt.plot(levr, trfatlhef3, c = 'r', linestyle ='--')
plt.plot(levr, trfatlwfo3, c = 'r', linestyle =':')

plt.xlabel('sigma_n', fontsize=14)
plt.ylabel('Tranformation(Sv)', fontsize=14)
plt.hlines(0.,sigmin, sigmax)

plt.legend(loc='upper left', title='', fontsize=10)

plt.text(23, trfmax+3, 'IPSL-CM6A-LR 1950-2009 DJF Surface transformation North Atl. > 40N', fontsize=14, fontweight='bold')

plt.show()


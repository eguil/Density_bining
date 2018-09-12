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
inDir = '/data/ericglod/Runs_nudges_sept2018'

run1 = 'CM6-pace-TSTr7fgT'
file1  = run1+'_1950_2009_mean_transf_north.nc'

#
# -- Open netcdf files
nc1 = open_ncfile(inDir + '/' + file1)

# -- Read variables for NorthAtl
trfatltot1 = npy.ma.average(nc1.variables['trsftotAtl'][:,:].squeeze(),axis=0)
trfatlhef1 = npy.ma.average(nc1.variables['trsfhefAtl'][:,:].squeeze(),axis=0)
trfatlwfo1 = npy.ma.average(nc1.variables['trsfwfoAtl'][:,:].squeeze(),axis=0)

# axis
levr = nc1.variables['rhon'][:]
sigmin = 22
sigmax = 30
trfmin = -15
trfmax = 30

plt.axis([sigmin, sigmax, trfmin, trfmax])

plt.plot(levr, trfatltot1, c = 'b', label = run1)
plt.plot(levr, trfatlhef1, c = 'b', linestyle ='--')
plt.plot(levr, trfatlwfo1, c = 'b', linestyle =':')

plt.xlabel('sigma_n', fontsize=14)
plt.ylabel('Tranformation(Sv)', fontsize=14)
plt.hlines(0.,sigmin, sigmax)

plt.legend(loc='upper left', title='', fontsize=10)

plt.text(23, 32, 'Surface transformation North Atl.', fontsize=14, fontweight='bold')

plt.show()


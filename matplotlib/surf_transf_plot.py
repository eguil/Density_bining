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
file1  = run1+'_1950_2009_mean_transf_north40.nc'
run2 = 'CM6-pace-TSTr7fgT'
file2  = run1+'_1950_2009_seasmean_transf_north40.nc'

#
# -- Open netcdf files
nc1 = open_ncfile(inDir + '/' + file1)
nc2 = open_ncfile(inDir + '/' + file2)

# -- Read variables for North Atl (anual mean from cdo yearmean on monthly data)
trfatltot1 = npy.ma.average(nc1.variables['trsftotAtl'][:,:].squeeze(),axis=0)
trfatlhef1 = npy.ma.average(nc1.variables['trsfhefAtl'][:,:].squeeze(),axis=0)
trfatlwfo1 = npy.ma.average(nc1.variables['trsfwfoAtl'][:,:].squeeze(),axis=0)

#trfatltot2 = npy.ma.average(nc2.variables['trsftotAtl'][:,:].squeeze(),axis=0)
#trfatlhef2 = npy.ma.average(nc2.variables['trsfhefAtl'][:,:].squeeze(),axis=0)
#trfatlwfo2 = npy.ma.average(nc2.variables['trsfwfoAtl'][:,:].squeeze(),axis=0)

trfatltot2 = nc2.variables['trsftotAtl'][0,:].squeeze()
trfatlhef2 = nc2.variables['trsfhefAtl'][0,:].squeeze()
trfatlwfo2 = nc2.variables['trsfwfoAtl'][0,:].squeeze()


# axis
levr = nc1.variables['rhon'][:]
sigmin = 22
sigmax = 30
trfmin = -15
trfmax = 30

plt.axis([sigmin, sigmax, trfmin, trfmax])

plt.plot(levr, trfatltot1, c = 'b', label = run1 + " >40 annual")
plt.plot(levr, trfatlhef1, c = 'b', linestyle ='--')
plt.plot(levr, trfatlwfo1, c = 'b', linestyle =':')

plt.plot(levr, trfatltot2, c = 'r', label = run2 + " >40 DJF")
plt.plot(levr, trfatlhef2, c = 'r', linestyle ='--')
plt.plot(levr, trfatlwfo2, c = 'r', linestyle =':')

plt.xlabel('sigma_n', fontsize=14)
plt.ylabel('Tranformation(Sv)', fontsize=14)
plt.hlines(0.,sigmin, sigmax)

plt.legend(loc='upper left', title='', fontsize=10)

plt.text(23, 32, 'Surface transformation North Atl.', fontsize=14, fontweight='bold')

plt.show()


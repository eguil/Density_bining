#!/bin/env python
# -*- coding: utf-8 -*-
"""
Python matplotlib
Make surface transformation plots

(c) Eric Guilyardi Sept. 2018

"""
from netCDF4 import Dataset as open_ncfile
import matplotlib.pyplot as plt
import numpy as np

# -------------------------------------------------------------------------------
#                                Define work
# -------------------------------------------------------------------------------

#inDir = '/Users/ericg/Projets/Density_bining/Nudge_CMIP6_IPSL/'
inDir = '/data/ericglod/Runs_nudges_sept2018'

file1  = 'CM6-pace-TSTr7fgT_1950_2009_mean_transf_north.nc'

#
# -- Open netcdf files
nc1 = open_ncfile(inDir + '/' + file2dh)

# -- Read variables for NorthAtl
trfatltot1 = nc1.variables['trsftotAtl'][:]
trfatlhef1 = nc1.variables['trsfhefAtl'][:]
trfatlwfo1 = nc1.variables['trsfwfoAtl'][:]

# axis
levr = nc1.variables['rhon'][:]


plt.plot(levr, trfatltot1)

plt.show()


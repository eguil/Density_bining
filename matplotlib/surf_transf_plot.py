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

# use cdo yseasmean to compute seasonal climatology
run1 = 'CM61-LR-hist-03.2110'
file1  = run1+'_1950_2009_seasmean_transf_north40.nc'
run2 = 'CM6-pace-TSTr8fgT'
file2  = run2+'_1950_2009_seasmean_transf_north40.nc'
file2c  = run2+'_1950_2009_seasmean_transf_north40_corr.nc'
run3 = 'CM6-pace-TSTr8vg'
file3  = run3+'_1950_1999_seasmean_transf_north40.nc'
file3c  = run3+'_1950_1999_seasmean_transf_north40_corr.nc'
run4 = 'CM6-pace-TSTr8vgS0'
file4  = run4+'_1950_2009_seasmean_transf_north40.nc'
file4c  = run4+'_1950_2009_seasmean_transf_north40_corr.nc'

#
# -- Open netcdf files
nc1 = open_ncfile(inDir + run1 + '/' + file1)
nc2 = open_ncfile(inDir + run2 + '/' + file2)
nc3 = open_ncfile(inDir + run3 + '/' + file3)
nc4 = open_ncfile(inDir + run4 + '/' + file4)
nc2c = open_ncfile(inDir + run2 + '/' + file2c)
nc3c = open_ncfile(inDir + run3 + '/' + file3c)
nc4c = open_ncfile(inDir + run4 + '/' + file4c)

# -- Read variables for North Atl (anual mean from cdo yearmean on monthly data)

trfatltot1 = nc1.variables['trsftotAtl'][0,:].squeeze()
trfatlhef1 = nc1.variables['trsfhefAtl'][0,:].squeeze()
trfatlwfo1 = nc1.variables['trsfwfoAtl'][0,:].squeeze()

trfatltot2 = nc2.variables['trsftotAtl'][0,:].squeeze()
trfatlhef2 = nc2.variables['trsfhefAtl'][0,:].squeeze()
trfatlwfo2 = nc2.variables['trsfwfoAtl'][0,:].squeeze()

trfatltot2c = nc2c.variables['trsftotAtl'][0,:].squeeze()
trfatlhef2c = nc2c.variables['trsfhefAtl'][0,:].squeeze()
trfatlwfo2c = nc2c.variables['trsfwfoAtl'][0,:].squeeze()

trfatltot3 = nc3.variables['trsftotAtl'][0,:].squeeze()
trfatlhef3 = nc3.variables['trsfhefAtl'][0,:].squeeze()
trfatlwfo3 = nc3.variables['trsfwfoAtl'][0,:].squeeze()

trfatltot3c = nc3c.variables['trsftotAtl'][0,:].squeeze()
trfatlhef3c = nc3c.variables['trsfhefAtl'][0,:].squeeze()
trfatlwfo3c = nc3c.variables['trsfwfoAtl'][0,:].squeeze()

trfatltot4 = nc4.variables['trsftotAtl'][0,:].squeeze()
trfatlhef4 = nc4.variables['trsfhefAtl'][0,:].squeeze()
trfatlwfo4 = nc4.variables['trsfwfoAtl'][0,:].squeeze()

trfatltot4c = nc4c.variables['trsftotAtl'][0,:].squeeze()
trfatlhef4c = nc4c.variables['trsfhefAtl'][0,:].squeeze()
trfatlwfo4c = nc4c.variables['trsfwfoAtl'][0,:].squeeze()


# axis
levr = nc1.variables['rhon'][:]
sigmin = 22
sigmax = 29
trfmin = -10
trfmax = 40
trfdiff = 15


# -- Create figure and axes instances
fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(22, 8))
ax0 = axes[0]
ax1 = axes[1]
ax1 = axes[2]


plt.xlabel('sigma_n', fontsize=14)
plt.ylabel('Tranformation(Sv)', fontsize=14)


ax0.axis([sigmin, sigmax, trfmin, trfmax])

ax0.plot(levr, trfatltot1, c = 'orange', label = run1)
ax0.plot(levr, trfatlhef1, c = 'orange', linestyle ='--')
ax0.plot(levr, trfatlwfo1, c = 'orange', linestyle ='-.')

ax0.plot(levr, trfatltot2, c = 'r', label = run2)
ax0.plot(levr, trfatlhef2, c = 'r', linestyle ='--')
ax0.plot(levr, trfatlwfo2, c = 'r', linestyle ='-.')

ax0.plot(levr, trfatltot3, c = 'black', label = run3)
ax0.plot(levr, trfatlhef3, c = 'black', linestyle ='--')
ax0.plot(levr, trfatlwfo3, c = 'black', linestyle ='-.')

ax0.plot(levr, trfatltot4, c = 'g', label = run4)
ax0.plot(levr, trfatlhef4, c = 'g', linestyle ='--')
ax0.plot(levr, trfatlwfo4, c = 'g', linestyle ='-.')



plt.xlabel('sigma_n', fontsize=14)
plt.ylabel('Tranformation(Sv)', fontsize=14)
ax0.hlines(0.,sigmin, sigmax)

ax0.legend(loc='upper left', title='', fontsize=10)

# difference
ax1.axis([sigmin, sigmax, -trfdiff, trfdiff])

ax1.plot(levr, trfatltot2-trfatltot1, c = 'b', label = run2)
ax1.plot(levr, trfatlhef2-trfatlhef1, c = 'b', linestyle ='--')
ax1.plot(levr, trfatlwfo2-trfatlwfo1, c = 'b', linestyle ='-.')

ax1.plot(levr, trfatltot3-trfatltot1, c = 'r', label = run3)
ax1.plot(levr, trfatlhef3-trfatlhef1, c = 'r', linestyle ='--')
ax1.plot(levr, trfatlwfo3-trfatlwfo1, c = 'r', linestyle ='-.')

ax1.plot(levr, trfatltot4-trfatltot1, c = 'g', label = run4)
ax1.plot(levr, trfatlhef4-trfatlhef1, c = 'g', linestyle ='--')
ax1.plot(levr, trfatlwfo4-trfatlwfo1, c = 'g', linestyle ='-.')

ax1.hlines(0.,sigmin, sigmax)

ax1.legend(loc='upper left', title='', fontsize=10)

ax1.set_title("Differences")

# Correction
ax2.axis([sigmin, sigmax, -trfdiff, trfdiff])

ax2.plot(levr, trfatltot2c, c = 'b', label = run2)
ax2.plot(levr, trfatlhef2c, c = 'b', linestyle ='--')
ax2.plot(levr, trfatlwfo2c, c = 'b', linestyle ='-.')

ax2.plot(levr, trfatltot3c, c = 'r', label = run3)
ax2.plot(levr, trfatlhef3c, c = 'r', linestyle ='--')
ax2.plot(levr, trfatlwfo3c, c = 'r', linestyle ='-.')

ax2.plot(levr, trfatltot4c, c = 'g', label = run4)
ax2.plot(levr, trfatlhef4c, c = 'g', linestyle ='--')
ax2.plot(levr, trfatlwfo4c, c = 'g', linestyle ='-.')

ax2.hlines(0.,sigmin, sigmax)

ax2.legend(loc='upper left', title='', fontsize=10)

ax2.set_title("Correction (hfcorr/wfcorr)")


plt.xlabel('sigma_n', fontsize=14)
plt.ylabel('Tranformation(Sv)', fontsize=14)


ttxt = fig.suptitle('IPSL-CM6A-LR 1950-2009 DJF Surface transformation North Atl. > 40N', fontsize=14, fontweight='bold')
plt.show()


#!/bin/env python
# -*- coding: utf-8 -*-

"""
Python matplotlib
Make time of emergence PDF bar plots from historical vs. historicalNat
Choose which variable and domain to work on

"""

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset as open_ncfile
from maps_matplot_lib import defVarmme, averageDom
from modelsDef import defModelsCO2piC
from libToE import findToE, ToEdomainhistvshistNat

# ----- Workspace ------

indir_1pctCO2 = '/data/ericglod/Density_binning/Prod_density_april15/mme_1pctCO2/'
indir_piC = '/data/ericglod/Density_binning/Prod_density_april15/mme_piControl/'

models = defModelsCO2piC()

# ----- Work ------

varname = defVarmme('salinity'); v = 'S'
#varname = defVarmme('temp'); v = 'T'
#varname = defVarmme('depth'); v = 'Z'

# -- Choose method for computing ToE
method = 'average_ToE' # Determine 2D lat/rho ToE then average in the box
#method = 'average_signal' # Average signal and noise in the box, then compute ToE (much faster)


domain_name = 'Northern ST'
# 'Southern ST', 'SO', 'Northern ST', 'North Atlantic, 'North Pacific'
print domain_name

multStd = 2. # detect ToE at multStd std dev of piControl

labBowl = ['histNat', 'hist']

valmask = 1.e20

iniyear = 1860
finalyear = 2005
deltay = 10.

# density domain
rhomin = 21
rhomid = 26
rhomax = 28
domrho = [rhomin, rhomid, rhomax]

# ----- Variables ------

# Choose random file to read only the basic variables and properties common to all files
file = 'cmip5.' + models[0]['name'] + '.1pctCO2.ensm.an.ocn.Omon.density.ver-' + models[0]['file_end_CO2'] + '_zon2D.nc'
f = open_ncfile(indir_1pctCO2 + file,'r')

lat = f.variables['latitude'][:]; latN = lat.size
density = f.variables['lev'][:]; levN = density.size
time = f.variables['time'][:]; timN = time.size
var = varname['var_zonal']

# Define variable properties
minmax = varname['minmax_zonal']
legVar = varname['legVar']
unit = varname['unit']

#!/bin/env python
# -*- coding: utf-8 -*-

"""
Python matplotlib
Plot the evolution in time of the Southern Subtropics signal in 1pctCO2 vs PiControl runs
"""

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset as open_ncfile
from maps_matplot_lib import defVarmme, custom_div_cmap, averageDom
from modelsDef import defModelsCO2piC
from libToE import ToEdomain1pctCO2vsPiC

# ----- mme -----

# mme 1pct CO2
indir_1pctCO2 = '/data/ericglod/Density_binning/Prod_density_april15/mme_1pctCO2/'
file = 'cmip5.multimodel_piCtl.1pctCO2.ensm.an.ocn.Omon.density_zon2D.nc'
data = indir_1pctCO2 + file
fCO2mme = open_ncfile(data,'r')

# mme PiControl
indir_piC = '/data/ericglod/Density_binning/Prod_density_april15/mme_piControl/'
file = 'cmip5.multimodel_1pct.piControl.ensm.an.ocn.Omon.density_zon2D.nc'
data = indir_piC + file
fpiCmme = open_ncfile(data,'r')

# Read main variables
lat = fCO2mme.variables['latitude'][:]
time = fCO2mme.variables['time'][:]
varname = defVarmme('salinity'); v = 'S'
density = fCO2mme.variables['lev'][:]
var = varname['var_zonal_w/bowl']

models = defModelsCO2piC()
var_change = np.ma.masked_all((len(models)+1, len(time)))

domain_name = 'Southern ST'
box = [-38,-18,25,26.4]
varCO2mme = fCO2mme.variables[var][:,2,:,:].squeeze()
varpiCmme = fpiCmme.variables[var][-20:,2,:,:].squeeze()
varpiCmme = np.ma.average(varpiCmme, axis=0)
var_change[-1,:] = averageDom(varCO2mme, 3, box, lat, density) - averageDom(varpiCmme, 2, box, lat, density)


# ----- Individual runs -----
indir = '/data/ericglod/Density_binning/Prod_density_april15/'


for imodel in range(len(models)) :
    model = models[imodel]
    file_CO2 = 'mme_1pctCO2/cmip5.' + model['name'] + '.1pctCO2.ensm.an.ocn.Omon.density.ver-' + model['file_end_CO2'] + '_zon2D.nc'
    file_piC = 'mme_piControl/cmip5.' + model['name'] + '.piControl.ensm.an.ocn.Omon.density.ver-' + model['file_end_piC'] + '_zon2D.nc'
    fCO2 = open_ncfile(indir + file_CO2, 'r')
    fpiC = open_ncfile(indir + file_piC, 'r')

    # ----- Variables -----
    domain = ToEdomain1pctCO2vsPiC(model['name'], domain_name)[0]
    box = domain['Pacific']
    if box != None:
        varCO2 = fCO2.variables[var][:,2,:,:].squeeze()
        varpiC = fpiC.variables[var][-20:,2,:,:].squeeze()
        varpiC = np.ma.average(varpiC, axis=0)
        var_change[imodel,:] = averageDom(varCO2, 3, box, lat, density) - averageDom(varpiC, 2, box, lat, density)


# ----- Plot -----

fig = plt.figure()

for imodel in range(len(models)):
    plt.plot(time,var_change[imodel,:])
plt.plot(time,var_change[-1,:], color='red', linewidth=3, label='mme')
plt.axhline(y=0, color='black', ls='--')
plt.xlabel('Years', fontweight='bold')
plt.ylabel('Salinity change (PSU)', fontweight='bold')

plt.title('Evolution of Salinity change in the Southern ST for 1pctCO2 vs. PiControl',
          fontweight='bold', va='bottom', fontsize=13)

# plt.show()
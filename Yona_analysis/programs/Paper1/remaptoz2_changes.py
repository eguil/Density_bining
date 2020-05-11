#!/bin/env python
# -*- coding: utf-8 -*-

""" Choose a signal to compute, remap it from sigma back to z, using a climatological reconstructed pseudo-depth, and plot in pseudo-depth/lat coordinates
"""
import sys
sys.path.append('/home/ysilvy/Density_bining/Yona_analysis/programs/')
import os, glob
from netCDF4 import Dataset as open_ncfile
import matplotlib.pyplot as plt
from maps_matplot_lib import defVarmme, defVarDurack, zon_2Dz, custom_div_cmap, modelagree
from lib_remapping import remaptoz
import numpy as np
import datetime
import pickle

# -------------------------------------------------------------------------------
#                                Define work
# -------------------------------------------------------------------------------

# -- Choose what to compute
# name = 'mme_hist_histNat'
# name = 'mme_hist'
name = 'mme_1pctCO2vsPiC'
# name = 'mme_rcp85_histNat'
# name = 'ens_mean_hist'

# -- Choose which model to compute for ensemble means
model_name = 'MIROC-ESM-CHEM'

# -- Choose where to stop for 1%CO2 simulations : 2*CO2 (70 years) or 4*CO2 (140 years) or 1.4*CO2 (34 years)
focus_1pctCO2 = '2*CO2'  # 1.4 or 2*CO2 or 4*CO2

# output format
# outfmt = 'view'
outfmt = 'save'

# Model agreement level
agreelev = 0.5
modelAgree = True

# Show isopycnals
show_isopyc = True

valmask = 1.e20
basinN = 4

# -- Choose work files

if name == 'mme_hist_histNat':
    indirh = '/data/ysilvy/Density_binning/mme_hist/'
    fileh_2d = 'cmip5.multimodel_Nat_rcp85.historical.ensm.an.ocn.Omon.density_zon2D.nc'
    fileh_1d = 'cmip5.multimodel_Nat_rcp85.historical.ensm.an.ocn.Omon.density_zon1D.nc'
    datah_2d = indirh + fileh_2d; datah_1d = indirh + fileh_1d
    indirhn = '/data/ysilvy/Density_binning/mme_histNat/'
    filehn_2d = 'cmip5.multimodel_Nat_rcp85.historicalNat.ensm.an.ocn.Omon.density_zon2D.nc'
    filehn_1d = 'cmip5.multimodel_Nat_rcp85.historicalNat.ensm.an.ocn.Omon.density_zon1D.nc'
    datahn_2d = indirhn + filehn_2d; datahn_1d = indirhn + filehn_1d
    fh2d = open_ncfile(datah_2d,'r')
    fh1d = open_ncfile(datah_1d,'r')
    fhn2d = open_ncfile(datahn_2d,'r')
    fhn1d = open_ncfile(datahn_1d,'r')

if name == 'mme_hist':
    indir = '/data/ysilvy/Density_binning/mme_hist/'
    file_2d = 'cmip5.multimodel_Nat_rcp85.historical.ensm.an.ocn.Omon.density_zon2D.nc'
    file_1d = 'cmip5.multimodel_Nat_rcp85.historical.ensm.an.ocn.Omon.density_zon1D.nc'
    data_2d = indir + file_2d
    data_1d = indir + file_1d
    fh2d = open_ncfile(data_2d, 'r')
    fh1d = open_ncfile(data_1d, 'r')

if name == 'ens_mean_hist':
    indir = '/data/ericglod/Density_binning/Prod_density_april15/mme_hist/'
    if model_name == 'FGOALS-g2':
        indir = '/data/ysilvy/Density_binning/mme_hist/'
    file_2d = glob.glob(indir + '*.'+model_name+'*_zon2D.nc')[0]
    file_1d = glob.glob(indir + '*.'+model_name+'*_zon1D.nc')[0]
    fh2d = open_ncfile(file_2d, 'r')
    fh1d = open_ncfile(file_1d, 'r')

if name == 'mme_1pctCO2vsPiC':
    indir_1pctCO2 = '/data/ericglod/Density_binning/Prod_density_april15/mme_1pctCO2/'
    file_2d = 'cmip5.multimodel_piCtl.1pctCO2.ensm.an.ocn.Omon.density_zon2D.nc'
    file_1d = 'cmip5.multimodel_piCtl.1pctCO2.ensm.an.ocn.Omon.density_zon1D.nc'
    data_2d = indir_1pctCO2 + file_2d
    data_1d = indir_1pctCO2 + file_1d
    fh2d = open_ncfile(data_2d,'r')
    fh1d = open_ncfile(data_1d,'r')

    indir_piC = '/data/ericglod/Density_binning/Prod_density_april15/mme_piControl/'
    file_2d = 'cmip5.multimodel_1pct.piControl.ensm.an.ocn.Omon.density_zon2D.nc'
    file_1d = 'cmip5.multimodel_1pct.piControl.ensm.an.ocn.Omon.density_zon1D.nc'
    data_2d = indir_piC + file_2d
    data_1d = indir_piC + file_1d
    fhn2d = open_ncfile(data_2d,'r')
    fhn1d = open_ncfile(data_1d,'r')

if name == 'mme_rcp85_histNat':
    indir_rcp85 = '/data/ericglod/Density_binning/Prod_density_april15/mme_rcp85/'
    filercp85_2d = 'cmip5.multimodel_Nat.rcp85.ensm.an.ocn.Omon.density_zon2D.nc'
    filercp85_1d = 'cmip5.multimodel_Nat.rcp85.ensm.an.ocn.Omon.density_zon1D.nc'
    indirhn = '/data/ysilvy/Density_binning/mme_histNat/'
    filehn_2d = 'cmip5.multimodel_Nat_rcp85.historicalNat.ensm.an.ocn.Omon.density_zon2D.nc'
    filehn_1d = 'cmip5.multimodel_Nat_rcp85.historicalNat.ensm.an.ocn.Omon.density_zon1D.nc'
    fh2d = open_ncfile(indir_rcp85 + filercp85_2d, 'r')
    fh1d = open_ncfile(indir_rcp85 + filercp85_1d, 'r')
    fhn2d = open_ncfile(indirhn + filehn_2d, 'r')
    fhn1d = open_ncfile(indirhn + filehn_1d, 'r')

# -------------------------------------------------------------------------------
#                                Build variables
# -------------------------------------------------------------------------------

# == Define variables ==

# -- Salinity or temperature
varname = defVarmme('salinity'); v = 'S'
# varname = defVarmme('temp'); v = 'T'

var = varname['var_zonal_w/bowl'] #var_zonal_w/bowl Take all then mask because fields without the bowl are wrong

if name == 'mme_rcp85_histNat':
    minmax = varname['minmax_zonal_rcp85']
else:
    minmax = varname['minmax_zonal']
clevsm = varname['clevsm_zonal']
legVar = varname['legVar']
unit = varname['unit']

# == Read variables ==
lat = fh2d.variables['latitude'][:]
density = fh2d.variables['lev'][:]
basin = fh2d.variables['basin'][:]

# Repeat density into (basin,density,latitude) to remap isopycnal surfaces
lat2d,density2d = np.meshgrid(lat,density)
density3d = np.repeat(density2d[np.newaxis,:,:],4,axis=0)

if name == 'mme_hist_histNat' or name == 'mme_rcp85_histNat':
    field2r = fh2d.variables[var][-20:,:,:,:] # historical or RCP8.5
    bowl2z = fh1d.variables['ptopdepth'][-20:,:,:]
    bowl2r = fh1d.variables['ptopsigma'][-20:,:,:]
    if name == 'mme_hist_histNat':
        field1r = fhn2d.variables[var][-20:,:,:] # historicalNat (last 20 years to average)
        bowl1z = fhn1d.variables['ptopdepth'][-20:,:,:]
        bowl1r = fhn1d.variables['ptopsigma'][-20:,:,:]
    else :
        field1r = fhn2d.variables[var][:,:,:] # historicalNat (entire time series to average)
        bowl1z = fhn1d.variables['ptopdepth'][:,:,:]
        bowl1r = fhn1d.variables['ptopsigma'][:,:,:]
    if name == 'mme_rcp85_histNat':
        labBowl = ['histNat', 'RCP8.5']
    else:
        labBowl = ['histNat', 'hist']

if name == 'mme_hist' or name == 'ens_mean_hist':
    field2r = fh2d.variables[var][-5:,:,:,:]
    field1r = fh2d.variables[var][88:93,:,:,:] # 1950
    bowl2z = fh1d.variables['ptopdepth'][-5:,:,:]
    bowl2r = fh1d.variables['ptopsigma'][-5:,:,:]
    bowl1z = fh1d.variables['ptopdepth'][88:93,:,:]
    bowl1r = fh1d.variables['ptopsigma'][88:93,:,:]
    labBowl = ['1950','2000']
    if modelAgree:
        var_name = varname['var_zonal_w/bowl']
        var_agreer = fh2d.variables[var_name + 'Agree'][-5:,:,:,:]

if name == 'mme_1pctCO2vsPiC':
    if focus_1pctCO2 == '4*CO2':
        y1 = 134
        y2 = 140
    if focus_1pctCO2 == '2*CO2':
        y1 = 69
        y2 = 74
    if focus_1pctCO2 == '1.4*CO2':
        y1 = 33
        y2 = 38
    field2r = fh2d.variables[var][y1:y2,:,:,:] # 1pctCO2
    var_piC = varname['var_zonal_w/bowl'] # Problem with varBowl fields: no data, so read full field
    field1r = fhn2d.variables[var_piC][-20:,:,:,:] # PiControl
    bowl2z = fh1d.variables['ptopdepth'][y1:y2,:,:]
    bowl2r = fh1d.variables['ptopsigma'][y1:y2,:,:]
    bowl1z = fhn1d.variables['ptopdepth'][-10:,:,:]
    bowl1r = fhn1d.variables['ptopsigma'][-10:,:,:]
    labBowl = ['PiControl', focus_1pctCO2]

# == Compute signal hist - histNat or 1pctCO2 - PiControl or rcp8.5 - histNat ==
vardiffr = np.ma.average(field2r, axis=0) - np.ma.average(field1r, axis=0)

# == Average other variables ==
bowl2z = np.ma.average(bowl2z, axis=0)
bowl1z = np.ma.average(bowl1z, axis=0)
bowl2r = np.ma.average(bowl2r, axis=0)
bowl1r = np.ma.average(bowl1r, axis=0)
if name=='mme_hist' and modelAgree :
    var_agreer = np.ma.average(var_agreer, axis=0)

# == Mask above bowl ==
for ilat in range(len(lat)):
    if np.ma.is_masked(bowl2r[1,ilat]) == False :
        inda = np.ma.nonzero(bowl2r[1,ilat]>=density)
        vardiffr[1,inda,ilat] = np.ma.masked
        if modelAgree and name =='mme_hist':
            var_agreer[1,inda,ilat] = np.ma.masked
    if np.ma.is_masked(bowl2r[2,ilat]) == False :
        indp = np.ma.nonzero(bowl2r[2,ilat]>=density)
        vardiffr[2,indp,ilat] = np.ma.masked
        if modelAgree and name =='mme_hist':
            var_agreer[2,indp,ilat] = np.ma.masked
    if np.ma.is_masked(bowl2r[3,ilat]) == False :
        indi = np.ma.nonzero(bowl2r[3,ilat]>=density)
        vardiffr[3,indi,ilat] = np.ma.masked
        if modelAgree and name =='mme_hist':
            var_agreer[3,indi,ilat] = np.ma.masked

    
# == Read reference pseudo-depth used for remapping ==
indir_z = '/home/ysilvy/Density_bining/Yona_analysis/data/remaptoz/'
file_z = 'EN4.pseudo_depth.zonal.pkl'
pseudo_depth = pickle.load( open( indir_z+file_z, "rb" ))

# == Target grid for remapping ==
# WOA13 grid
targetz = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80,
   85, 90, 95, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 375,
   400, 425, 450, 475, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950,
   1000, 1050, 1100, 1150, 1200, 1250, 1300, 1350, 1400, 1450, 1500, 1550,
   1600, 1650, 1700, 1750, 1800, 1850, 1900, 1950, 2000, 2100, 2200, 2300,
   2400, 2500, 2600, 2700, 2800, 2900, 3000, 3100, 3200, 3300, 3400, 3500,
   3600, 3700, 3800, 3900, 4000, 4100, 4200, 4300, 4400, 4500, 4600, 4700,
   4800, 4900, 5000, 5100, 5200, 5300, 5400, 5500]

# == Remap ==
fieldz, zbowl = remaptoz(vardiffr,pseudo_depth,targetz)
print('field remapped')
density_z, zbowl2 = remaptoz(density3d,pseudo_depth,targetz)
print('density remapped')
if name=='mme_hist' and modelAgree:
    var_agreez, zbowl3 = remaptoz(var_agreer,pseudo_depth,targetz)

# import xarray as xr
# density_z = xr.DataArray(density_z,dims=['basin','pseudo-depth','latitude'],
#                           coords=[[0,1,2,3],targetz,lat])
# density_z = density_z.where(xr.DataArray(targetz,dims=['pseudo-depth'],coords=[targetz])>=zbowl)

# -- Make variable bundles for each basin
varAtl = {'name': 'Atlantic', 'var_change': fieldz[1,:,:], 'bowl1': None, 'bowl2': None,
          'labBowl': None, 'density':density_z[1,:,:]}
varPac = {'name': 'Pacific', 'var_change': fieldz[2,:,:], 'bowl1': None, 'bowl2': None,
          'labBowl': None, 'density':density_z[2,:,:]}
varInd = {'name': 'Indian', 'var_change': fieldz[3,:,:], 'bowl1': None, 'bowl2': None,
          'labBowl': None, 'density':density_z[3,:,:]}
    
# # -- Read bathymetry
# # Read masks
# fmask = open_ncfile('/home/ysilvy/Density_bining/Yona_analysis/data/170224_WOD13_masks.nc','r')
# basinmask = fmask.variables['basinmask3'][:] # (latitude, longitude)
# depthmask = fmask.variables['depthmask'][:] # (latitude, longitude)
# longitude = fmask.variables['longitude'][:]
# # Create basin masks
# mask_a = basinmask != 1
# mask_p = basinmask != 2
# mask_i = basinmask != 3
# # Read bathy
# depthmask_a = np.ma.array(depthmask, mask=mask_a) # Mask every basin except Atlantic
# depthmask_p = np.ma.array(depthmask, mask=mask_p)
# depthmask_i = np.ma.array(depthmask, mask=mask_i)
# # Zonal bathy
# bathy_a = np.ma.max(depthmask_a, axis=1)
# bathy_p = np.ma.max(depthmask_p, axis=1)
# bathy_i = np.ma.max(depthmask_i, axis=1)
# # Regroup
# bathy = np.ma.masked_all((basinN,len(lat)))
# bathy[1,:] = bathy_a
# bathy[2,:] = bathy_p
# bathy[3,:] = bathy_i

# -------------------------------------------------------------------------------
#                                Plot
# -------------------------------------------------------------------------------

domzed = [0,500,5000]

# -- Create figure and axes instances
fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(17, 5))

# -- color map
cmap = custom_div_cmap()

# -- levels
levels = np.linspace(minmax[0],minmax[1],minmax[2])

ext_cmap = 'both'
contourDict = {'cmap':cmap, 'levels':levels, 'levels2':levels, 'ext_cmap':ext_cmap, 'isopyc':show_isopyc}

# -- Contourf of signal
cnplot = zon_2Dz(plt, axes[0,0], axes[1,0], 'left', lat, targetz, varAtl,
                 contourDict, domzed)
cnplot = zon_2Dz(plt, axes[0,1], axes[1,1], 'mid', lat, targetz, varPac,
                 contourDict, domzed)
cnplot = zon_2Dz(plt, axes[0,2], axes[1,2], 'right', lat, targetz, varInd,
                 contourDict, domzed)

if name=='mme_hist' and modelAgree:
    modelagree(axes[0,0],axes[1,0],agreelev,lat,targetz,var_agreez[1,:,:])
    modelagree(axes[0,1],axes[1,1],agreelev,lat,targetz,var_agreez[2,:,:])
    modelagree(axes[0,2],axes[1,2],agreelev,lat,targetz,var_agreez[3,:,:])

# Bathymetry
#for i in range(3):
#    axes[0,i].fill_between(lat,[targetz[-1]]*len(lat),bathy[i+1,:],facecolor='0.2')
#    axes[1,i].fill_between(lat,[targetz[-1]]*len(lat),bathy[i+1,:],facecolor='0.2')

for i in range(3):
    # -- Draw bowl
    axes[0,i].fill_between(lat,y1=0,y2=zbowl[i+1],color='0.3',zorder=10)
    
plt.subplots_adjust(hspace=.012, wspace=0.05, left=0.05, right=0.86)

# -- Add colorbar
cb = plt.colorbar(cnplot[0], ax=axes.ravel().tolist(), ticks=levels[::3], fraction=0.015, shrink=2.0, pad=0.05)
cb.set_label('%s (%s)' % (legVar, unit), fontweight='bold',fontsize=14)
cb.ax.set_yticklabels(cb.ax.get_yticklabels(), fontweight='bold')
cb.ax.yaxis.set_tick_params(which='major',width=2)

# -- Add Title text
if name == 'mme_hist_histNat':
    plotTitle = 'Multi-model mean %s changes historical-historicalNat[1985-2005]' %(legVar,)
    plotName = 'remapping2_' + name + '_' + legVar+ '_last20years_paper'
    figureDir = 'models/zonal_remaptoz/'
    plt.figtext(.006,.96,'a',fontweight='bold',fontsize=16)
    #plt.figtext(.2,.01,'Last 20 years',fontsize=8,ha='center')
if name == 'mme_hist':
    plotTitle = 'Multi-model mean %s changes [2000-2005]-[1950-1955]' %(legVar,)
    plotName = 'remapping2_' + name + '_' + legVar + '_paper'
    figureDir = 'models/zonal_remaptoz/'
    plt.figtext(.006,.96,'b',fontweight='bold',fontsize=16)
if name == 'ens_mean_hist':
    plotTitle = '%s changes [2000-2005]-[1950-1955] for %s' %(legVar,model_name)
    plotName = 'remapping2_' + model_name + '_' + name+ '_' + legVar
    figureDir = 'models/zonal_remaptoz/'
if name == 'mme_1pctCO2vsPiC':
    plotTitle  = 'Multi-model mean %s changes %s-piControl' %(legVar,focus_1pctCO2)
    plotName = 'remapping2_' + name + '_' + focus_1pctCO2 + '_' + legVar
    figureDir = 'models/zonal_remaptoz/'
if name == 'mme_rcp85_histNat':
    plotTitle = 'Multi-model mean %s changes RCP8.5[2080-2100]-historicalNat' %(legVar,)
    plotName = 'remapping2_' + name + '_' + legVar + '_last20years_paper'
    figureDir = 'models/zonal_remaptoz/'
    plt.figtext(.006,.96,'b',fontweight='bold',fontsize=16)
    #plt.figtext(.2,.01,'Last 20 years RCP8.5 - mean(histNat)',fontsize=8,ha='center')

# axes[0,1].set_title(plotTitle, y=1.25, fontweight='bold', fontsize=15, verticalalignment='top')
#fig.suptitle(plotTitle, fontsize=14, fontweight='bold')
plt.figtext(.002,.35,'Pseudo-depth (m)',rotation='vertical',fontweight='bold',fontsize=14)



# Date
now = datetime.datetime.now()
date = now.strftime("%Y-%m-%d")

#plt.figtext(.5,.01,'Computed by : remaptoz2_changes.py,  '+date,fontsize=8,ha='center')
#  if name == 'mme_hist' and modelAgree:
#     plt.figtext(.2,.01,'Model agreement level : ' + str(agreelev),fontsize=8,ha='center')

if outfmt == 'view':
    plt.show()
else:
    plt.savefig(plotName+'.png', bbox_inches='tight',dpi=150) #'/home/ysilvy/figures/'+figureDir+

# -*- coding: utf-8 -*-

"""
Use EN4 climatological volume of isopycnals to construct a density/pressure relationship used for remapping data from density to a pseudo-depth coordinate
"""

import numpy as np
import matplotlib.pyplot as plt
import os, glob
from netCDF4 import Dataset as open_ncfile
import datetime

# ===
# === WORKSPACE AND PRE-REQUISITES ===
# ===

# Read EN4 file : will be our climatology
indir = '/home/ericglod/Density_bining/test/'
file = 'EN4.mon.ocean.Omon.1900_2017.density.nc'#'obs.EN4.historical.r0i0p0.mo.ocn.Omon.density.ver-1.latestX_zon2D.nc'
f = open_ncfile(indir+file,'r')

# Read variables
lat = f.variables['latitude'][:]
density = f.variables['lev'][:]
isonvol = np.ma.average(f.variables['isonvol'][:,:,:,:],axis=0)*1.e03 # Volume of isopycnals in km3, make climatology

# Dimensions
basinN=4
latN = len(lat)
densityN = len(density)

# Z grid for calculating ocean volume per depth level
#gridz = np.arange(0,5501,5)
gridz2 = np.concatenate([np.arange(0,21,1),np.arange(25,5501,5)])

# Read bathymetry
# Read masks
fmask = open_ncfile('/home/ysilvy/Density_bining/Yona_analysis/data/170224_WOD13_masks.nc','r')
basinmask = fmask.variables['basinmask3'][:] # (latitude, longitude)
depthmask = fmask.variables['depthmask'][:] # (latitude, longitude)
longitude = fmask.variables['longitude'][:]
# Create basin masks
mask_a = basinmask != 1
mask_p = basinmask != 2
mask_i = basinmask != 3
# Read bathy
depthmask_a = np.ma.array(depthmask, mask=mask_a) # Mask every basin except Atlantic
depthmask_p = np.ma.array(depthmask, mask=mask_p)
depthmask_i = np.ma.array(depthmask, mask=mask_i)
# Read area
area = fmask.variables['basinmask3_area'][:] # (latitude,longitude)

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

# ===
# === OCEAN VOLUME CALCULATION PER DEPTH LEVEL FROM BATHYMETRY ===
# ===

# Compute zonal volume on gridz2 levels (every 5 meters, every 1 meter at the surface)
V = np.ma.masked_all((basinN,len(gridz2),latN))

# Loop on Atlantic, Pacific and Indian basins
for ibasin in range(1,4): 
    print(ibasin)
    if ibasin==1:
        depthmask = depthmask_a
    elif ibasin==2:
        depthmask = depthmask_p
    else:
        depthmask = depthmask_i
    for ilat in range(latN):
        area_lat = np.ma.average(area[ilat,:]) # Same area for all longitudes along ilat: take the average to make sure we don't choose a masked gridpoint
        # Loop on depths
        for iz in range(len(gridz2)):
            cells_above = np.argwhere(depthmask[ilat,:]>gridz2[iz]) #.squeeze() # list of indices where bathymetry is below the current z level
            if len(cells_above)!=0: #Exclude depths where there is no water
                S = len(cells_above)*area_lat # total area in km2 at depth gridz[iz] for the current band of latitude
                V[ibasin,iz,ilat] = S*0.001*(gridz2[iz+1]-gridz2[iz]) # Corresponding volume in km3 (5m interval for gridz, 1m at surface)


# Now compute volume for global zonal mean (ibasin=0)
V[0,:,:] = np.sum(V[1:,:,:],axis=0)

# Computed volume of the ocean, based on WOA13 bathymetry
print('Total ocean volume from bathymetry: ', np.sum(V[0,:,:])/1.e09) # unit = *10^9km3

# EN4 volume of the ocean
print('EN4 total ocean volume: ', np.ma.sum(isonvol[0,:,:])/1.e09) #*10^9km3


# ========
# PSEUDO-Z CONSTRUCTION
# ========

# Constructing remapping relationship pseudo_depth[basin,density,latitude] by re-filling the ocean from the surface down with horizontal layers of constant density

# Initialize
pseudo_depth = np.ma.masked_all((basinN,densityN,latN))
# Start loops
for ibasin in range(basinN):
    for ilat in range(latN):
        if not np.ma.is_masked(np.all(isonvol[ibasin,:,ilat])):
            idx_range = np.ma.flatnotmasked_edges(isonvol[ibasin,:,ilat]) # Indices edges of unmasked densities in the column
            bathy=np.ma.max(depthmask_a[ilat,:]) # Find max depth of the water column
            ibathy=np.argwhere(gridz2==bathy).squeeze() # Index of bathymetry = first masked value of the water column
            cum_V = np.cumsum(V[ibasin,:,ilat])
            cum_isonvol = np.cumsum(isonvol[ibasin,:,ilat])
            # First level
            pseudo_depth[ibasin,idx_range[0],ilat] = 0
            # Loop on all other density levels (stop one before the end)
            for ilev in range(idx_range[0],idx_range[1]):
                isonvol_lev = isonvol[ibasin,ilev,ilat]
                iz = find_nearest(cum_V,cum_isonvol[ilev])# Find index where the cumulated volume at ilev is closest to cumulated volume on gridz2
                pseudo_depth[ibasin,ilev+1,ilat] = gridz2[iz+1] # Save corresponding z level on gridz (+1 because volume of one level corresponds to the water between that level and the level below)
                #print(ilev,density[ilev],iz+1,gridz2[iz+1],cum_isonvol[ilev],cum_V[iz])
                
# Rename so as not to create confusion with netcdf variable
pseudo_z = pseudo_depth
lev = density

# ========
# SAVE PSEUDO-Z TO FILE
# ========

import pickle

# write to pickle
pickle.dump( pseudo_depth, open( "EN4.pseudo_depth.zonal.pkl", "wb" ) )


## Date
#now = datetime.datetime.now()
#date = now.strftime("%Y-%m-%d")
#
#fileName = 'EN4.pseudo_depth.zonal.nc'
#dir = '/home/ysilvy/Density_bining/Yona_analysis/data/remaptoz/'
#description = 'Pseudo_depth as a function of density created from EN4 climatological volume of isopycnals, by refilling the ocean from surface to bottom with horizontal layers of increasing density classes, per latitude band and per ocean basin. Created from ocean_refill_pseudoz_from_EN4.py on ' +date
#fout = open_ncfile(dir+fileName,'w', format='NETCDF4')
#fout.description = description

## dimensions
#fout.createDimension('basin', basinN)
#fout.createDimension('density', densityN)
#fout.createDimension('latitude', latN)

## variables
#basin = fout.createVariable('basin', 'i4', ('basin',))
#density = fout.createVariable('density', 'd', ('density',))
#latitude = fout.createVariable('latitude', 'f4', ('latitude',))
#pseudo_depth = fout.createVariable('pseudo_depth', 'f4', ('basin','density','latitude',))

## data
#basin[:] =  np.arange(0,basinN)
#density[:] = lev
#latitude[:] = lat
#pseudo_depth = pseudo_z

## units
#basin.units = 'basin index'
#density.units = 'kg.m-3'
#latitude.units = ''
#pseudo_depth.units = 'm'

#fout.close()
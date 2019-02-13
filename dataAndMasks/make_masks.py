#!/usr/local/uvcdat/2014-07-21/bin/python
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 22 19:22:06 2014

Paul J. Durack 22nd June 2014

This script interrogates ascii mask files for WOA13 and writes these to netcdf

PJD 22 Jun 2014     - Started
PJD 23 Jun 2014     - Hit issue with ubyte type in netcdf4, byte appears to be written not ubyte
                      kludgey work around implemented below
PJD  6 Aug 2014     - Problems with unsigned types and cdms2, converting to int16
PJD  6 Aug 2014     - Added mask with 3 basins and masked marginal seas
                    http://data.nodc.noaa.gov/woa/WOA13/DOC/woa13documentation.pdf
                    http://www-pcmdi.llnl.gov/publications/pdf/34.pdf
PJD  7 Aug 2014     - Resolved write problems by using int16 rather than uint8 (cdms2 unsupported type)
PJD  7 Aug 2014     - Added area weights
PJD  7 Aug 2014     - Updated area weights to float32 (was int16)
PJD 24 Feb 2017     - Updated to provide depth mask
                    - TODO:
                    - Convert matrix to be indexed from 20E (cut through Africa)
                    - Generate 0.25 degree (*04.msk) files

@author: durack1
"""

import datetime,gc,os,sys
import cdms2 as cdm
import cdutil as cdu
import numpy as np
import MV2 as mv
sys.path.append('/export/durack1/git/durolib/lib')
from durolib import globalAttWrite

# netCDF compression (use 0 for netCDF3)
cdm.setNetcdfShuffleFlag(1)
cdm.setNetcdfDeflateFlag(1)
cdm.setNetcdfDeflateLevelFlag(9) ; # 9(shuf=1) 466.6KB; 9(shuf=0) 504.1KB; 4(shuf=0) 822.93KB; 
cdm.setAutoBounds(1)

#%% Change directory
os.chdir('/work/durack1/Shared/obs_data/WOD13/')

#%%
#del(asc,count,depths,e,lat_ind,latitude,line,lon_ind,longitude,pi,tmp)
#del(asc,count,depths,lat_ind,latitude,line,lon_ind,longitude,tmp)

#%% Deal with 1-degree grid
latitude    = np.arange(-89.5,90.5,1,dtype='float32')
longitude   = np.arange(-179.5,180.5,1,dtype='float32')
depth       = np.float32([0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,125,150,175,200,225,250,275,300,325,350,375,400,\
               425,450,475,500,550,600,650,700,750,800,850,900,950,1000,1050,1100,1150,1200,1250,1300,1350,1400,1450,1500,1550,1600,\
               1650,1700,1750,1800,1850,1900,1950,2000,2100,2200,2300,2400,2500,2600,2700,2800,2900,3000,3100,3200,3300,3400,3500,\
               3600,3700,3800,3900,4000,4100,4200,4300,4400,4500,4600,4700,4800,4900,5000,5100,5200,5300,5400,5500,5600,5700,5800,\
               5900,6000,6100,6200,6300,6400,6500,6600,6700,6800,6900,7000,7100,7200,7300,7400,7500,7600,7700,7800,7900,8000,8100,\
               8200,8300,8400,8500,8600,8700,8800,8900,9000]) ; # WOA13 standard 137 levels
# 0-5500 metres 102 levels
# 0-1500 metres 57 levels
# 0-500  metres 37 levels

# Deal with basin through depth
#%% basinmask
asc         = np.genfromtxt('140622/basinmask_01.msk',dtype=None)
basinmask   = np.ma.zeros([len(depth),len(latitude),len(longitude)],dtype='int16')
for count,line in enumerate(asc):
    if count == 0:
        continue
    if not np.mod(count,10000):
        print count
    # Deal with basinmask
    tmp = line.split(',')
    tmp = filter(None,tmp)
    tmp_len = len(tmp)-2
    lat_ind = np.where(np.equal(latitude,np.float32(tmp[0])))[0][0]
    lon_ind = np.where(np.equal(longitude,np.float32(tmp[1])))[0][0]
    # Assign to variables
    basinmask[0:tmp_len,lat_ind,lon_ind] = np.int16(tmp[2:])
del(asc,count,line,tmp,tmp_len,lat_ind,lon_ind) ; gc.collect()
# Convert to masked array
basinmask   = np.ma.masked_equal(basinmask,0)

#%% mixmask
asc         = np.genfromtxt('140622/mixnumber_01.msk',dtype=None)
mixmask     = basinmask.copy()
for count,line in enumerate(asc):
    if count == 0:
        continue
    if not np.mod(count,10000):
        print count
    # Deal with basinmask
    tmp = line.split(',')
    tmp = filter(None,tmp)
    tmp_len = len(tmp)-2
    lat_ind = np.where(np.equal(latitude,np.float32(tmp[0])))[0][0]
    lon_ind = np.where(np.equal(longitude,np.float32(tmp[1])))[0][0]
    # Deal with mixmask
    line = asc[count]
    tmp = line.split(',')
    tmp = filter(None,tmp)
    tmp_len = len(tmp)-2
    # Assign to variables
    mixmask[0:tmp_len,lat_ind,lon_ind] = np.int16(tmp[2:])
del(asc,count,line,tmp,tmp_len,lat_ind,lon_ind) ; gc.collect()
# Convert to masked array
mixmask     = np.ma.masked_equal(mixmask,0)

#%% Deal with basin at surface
asc         = np.genfromtxt('140622/landsea_01.msk',dtype=None)
landsea     = np.ma.zeros([len(latitude),len(longitude)],dtype='int16')
for count,line in enumerate(asc):
    if count == 0:
        continue
    if not np.mod(count,10000):
        print count
    tmp = line.split(',')
    tmp = filter(None,tmp)
    tmp_len = len(tmp)-2
    lat_ind = np.where(np.equal(latitude,np.float32(tmp[0])))[0][0]
    lon_ind = np.where(np.equal(longitude,np.float32(tmp[1])))[0][0]
    #if np.uint8(tmp[2:]) < 0:
    #    print tmp,np.uint8(tmp[2:])
    landsea[lat_ind,lon_ind] = np.int16(tmp[2:])
del(asc,count,line,tmp,tmp_len,lat_ind,lon_ind) ; gc.collect()
# Convert to masked array
landsea = np.ma.masked_equal(landsea,0)

# Write numpy to cdms2 objects
latitude    = cdm.createAxis(latitude,id='latitude')
longitude   = cdm.createAxis(longitude,id='longitude')
depth       = cdm.createAxis(depth,id='depth')
basinmask   = cdm.createVariable(basinmask,id='basinmask',axes=[depth,latitude,longitude])
mixmask     = cdm.createVariable(mixmask,id='mixmask',axes=[depth,latitude,longitude])
landsea     = cdm.createVariable(landsea,id='landsea',axes=[latitude,longitude])

# Variable attributes
basinmask.index = ''.join(['1: Atlantic Ocean; 2: Pacific Ocean; 3: Indian Ocean; 4: Mediterranean Sea; ',
                           '5: Baltic Sea; 6: Black Sea; 7: Red Sea; 8: Persian Gulf; 9: Hudson Bay; ',
                           '10: Southern Ocean; 11: Arctic Ocean; 12: Sea of Japan; 13: Kara Sea; ',
                           '14: Sulu Sea; 15: Baffin Bay; 16: East Mediterranean; ',
                           '17: West Mediterranean; 18: Sea of Okhotsk; 19: Banda Sea; 20: Caribbean Sea; ',
                           '21: Andaman Basin; 22: North Caribbean; 23: Gulf of Mexico; 24: Beaufort Sea; ',
                           '25: South China Sea; 26: Barents Sea; 27: Celebes Sea; 28: Aleutian Basin; ',
                           '29: Fiji Basin; 30: North American Basin; 31: West European Basin; ',
                           '32: Southeast Indian Basin; 33: Coral Sea; 34: East Indian Basin; ',
                           '35: Central Indian Basin; 36: Southwest Atlantic Basin; 37: Southeast Atlantic Basin; ',
                           '38: Southeast Pacific Basin; 39: Guatemala Basin; 40: East Caroline Basin; ',
                           '41: Marianas Basin; 42: Philippine Sea; 43: Arabian Sea; 44: Chile Basin; ',
                           '45: Somali Basin; 46: Mascarene Basin; 47: Crozet Basin; 48: Guinea Basin; ',
                           '49: Brazil Basin; 50: Argentine Basin; 51: Tasman Sea; 52: Atlantic Indian Basin; ',
                           '53: Caspian Sea; 54: Sulu Sea II; 55: Venezuela Basin; 56: Bay of Bengal; ',
                           '57: Java Sea; 58: East Indian Atlantic Basin;'])
landsea.info  = '0: land; 1-137: standard depth level of last valid data point (ocean bottom for cell)'
mixmask.info  = 'Shows basin interactions used in the objective analysis for each standard depth level'

# Create 3 ocean mask
basinmask3 = basinmask[0,] ; # Trim off top layer
# Try sftbyregion
#mask = cdu.generateLandSeaMask(basinmask3)
#sftb = cdu.sftbyrgn.generateSurfaceTypeByRegionMask(mask)
# Southern Ocean
# Case Atlantic
lonmat = np.ma.zeros([len(latitude),len(longitude)],dtype='int16')
lonlims = [-70.5,20.5]
loninds = [np.where(longitude==np.float32(lonlims[0])),np.where(longitude==np.float32(lonlims[1]))]
loninds = [loninds[0][0][0],loninds[1][0][0]]
lonmat[:,loninds[0]:loninds[1]] = 1
basinmask3 = mv.where(mv.logical_and(basinmask3==10,lonmat),1,basinmask3) ; # Case Atlantic
# Case Indian
lonmat = np.ma.zeros([len(latitude),len(longitude)],dtype='int16')
lonlims = [20.5,147.5]
loninds = [np.where(longitude==np.float32(lonlims[0])),np.where(longitude==np.float32(lonlims[1]))]
loninds = [loninds[0][0][0],loninds[1][0][0]]
lonmat[:,loninds[0]:loninds[1]] = 1
basinmask3 = mv.where(mv.logical_and(basinmask3==10,lonmat),3,basinmask3) ; # Case Indian
# Case Pacific
lonmat = np.ma.zeros([len(latitude),len(longitude)],dtype='int16')
lonlims = [147.5,179.5]
loninds = [np.where(longitude==np.float32(lonlims[0])),np.where(longitude==np.float32(lonlims[1]))]
loninds = [loninds[0][0][0],loninds[1][0][0]]
lonmat[:,loninds[0]:loninds[1]] = 1
lonlims = [-179.5,-69.5]
loninds = [np.where(longitude==np.float32(lonlims[0])),np.where(longitude==np.float32(lonlims[1]))]
loninds = [loninds[0][0][0],loninds[1][0][0]]
lonmat[:,loninds[0]:loninds[1]] = 1
basinmask3 = mv.where(mv.logical_and(basinmask3==10,lonmat),2,basinmask3) ; # Case Pacific
# Marginal seas
basinmask3 = mv.where(basinmask3==11,1,basinmask3) ; # Convert Arctic to Atlantic
basinmask3 = mv.where(basinmask3==56,3,basinmask3) ; # Convert Bay of Bengal to Indian
basinmask3 = mv.masked_where(basinmask3>3,basinmask3) ; # Convert all seas to missing
basinmask3 = cdm.createVariable(np.int16(basinmask3),id='basinmask3',axes=[latitude,longitude])
basinmask3.index = '1: Atlantic Ocean; 2: Pacific Ocean; 3: Indian Ocean;'

# Add area weights
earthSurfaceAreaKm2 = 510.1
earthWaterAreaKm2   = 361.132
earthLandAreaKm2    = 148.94
#grid = cdm.createGenericGrid(latitude,longitude)
grid = cdm.createVariable(np.ma.zeros([len(latitude),len(longitude)],),id='mask',axes=[latitude,longitude])
frac = cdu.area_weights(grid)
area = frac*(earthSurfaceAreaKm2*1e6) # ; Area m2
area = mv.masked_where(basinmask3.mask,area) # ; Masked area 344.3 km^2
basinmask3_area = cdm.createVariable(np.float32(area),id='basinmask3_area',axes=[latitude,longitude])
basinmask3_area.earthSurfaceAreaM2  = earthSurfaceAreaKm2*1e6
basinmask3_area.earthWaterAreaM2    = earthWaterAreaKm2*1e6
basinmask3_area.earthLandAreaM2     = earthLandAreaKm2*1e6
basinmask3_area.oceanSurfaceAreaM2  = area.sum()
basinmask3_area.units = 'm^2'

# Plot to check
#import vcs as vc
#v1 = vc.init()
#v1.plot(area)

# Fix missing value within range - useful in uint8 case
#vl = np.int16(32767) ; # problems with uint8 ; cdms2 is writing byte NOT ubyte
#print vl,type(vl)
#basinmask.set_fill_value(np.uint8(vl))
#basinmask._fill_value=vl
# setMissing fails with - TypeError: 'numpy.uint8' object does not support item assignment
# setMissing fails with - TypeError: 'numpy.int16' object does not support item assignment
# setMissing fails with - TypeError: 'numpy.float32' object does not support item assignment
#basinmask.setMissing=vl

#%% Generate depth mask from lowest value
depthmask = np.ma.zeros([len(latitude),len(longitude)],dtype='int16')
nonmasked = np.transpose(basinmask.nonzero()) ; # Returns depth,lat,lon matrix of nonzero elements
for count in range(0,len(nonmasked)):
    tmpInd = nonmasked[count]
    depthmask[tmpInd[1],tmpInd[2]] = depth[tmpInd[0]] ; # Save depth to indexing, overwriting shallower depths
del(count,tmpInd) ; gc.collect()
# Mask
depthmask = mv.masked_where(depthmask==0,depthmask) ; # Convert all seas to missing
depthmask = np.ma.masked_equal(depthmask,0)
# Convert to cdms variable
depthmask = cdm.createVariable(depthmask,id='depthmask',axes=[latitude,longitude])
depthmask.info = 'Grid point depth with maximum value set at 5500 metres (WOA13 standard grid)'

#%%
# Get time format
time_now = datetime.datetime.now()
time_format = time_now.strftime("%y%m%d")
outfile = '_'.join([time_format,'WOD13_masks.nc'])
# Create output netcdf files
print '** Writing file:',outfile
if os.path.isfile(outfile):
    os.remove(outfile) ; # purge existing file
fOut = cdm.open(outfile,'w')
# Write to outfile
fOut.write(basinmask)
fOut.write(basinmask3)
fOut.write(basinmask3_area)
fOut.write(depthmask)
fOut.write(mixmask)
fOut.write(landsea)
# Write global attributes
fOut.title  = 'NODC World Ocean Atlas 2013 basin masks and areas'
fOut.url    = 'http://www.nodc.noaa.gov/OC5/woa13/masks13.html'
fOut.info   = 'http://data.nodc.noaa.gov/woa/WOA13/DOC/woa13documentation.pdf'
globalAttWrite(fOut,options=None)
fOut.close()
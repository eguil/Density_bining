#!/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 24 11:50:23 2016

Paul J. Durack 24th April 2016

This script generates area weight (areacello) variables from input grids

PJD 24 Apr 2016     - Started
PJD 25 Apr 2016     - Corrected confused logic with IPRC - drive_obs not here
PJD 25 Apr 2016     - Added chdir to local dir
PJD  1 May 2019     - Updated to check areacello values
PJD  1 May 2019     - Update refs to durolib
PJD  1 May 2019     - Corrected issue with axis indexing of gnu.area_weights
                    Check value
                    Check oceanAreaM2:      361132000000000.0
                    IPRC  earthOceanAreaM2: 296726559170288.1
                    UCSD  earthOceanAreaM2: 312227029571367.06
                    SM07  earthOceanAreaM2: 362339513516552.0
                    Jam   earthOceanAreaM2: 306542126017143.9
                    EN4   earthOceanAreaM2: 380806685301225.25
                    - TODO:

@author: durack1
"""
import gc,glob,os,socket,sys
sys.path.append('/export/durack1/git/durolib/durolib')
from durolib import globalAttWrite
from string import replace
import cdms2 as cdm
import genutil as gnu
#import MV2 as mv
import numpy as np

# Set nc classic as outputs
cdm.setNetcdfShuffleFlag(0)
cdm.setNetcdfDeflateFlag(1)
cdm.setNetcdfDeflateLevelFlag(9) ; # Set to auto
cdm.setCompressionWarnings(False)
cdm.setAutoBounds(True)

#%% Set areas
earthAreaKm2 = 510.072e6 ; # Earth area in km^2
earthAreaM2 = earthAreaKm2*1e6 ; # Earth area km^2 -> m^2
earthWaterAreaKm2 = 361.132e6
earthLandAreaKm2 = 148.94e6

#%% Get input grids
os.chdir('/work/durack1/Shared/190213_data_density')
inFiles = glob.glob('*thetao*.xml')

#%% Determine machine and host files
host = socket.gethostname()
if 'ocean' in host:
    fileStr = 'obs'
elif 'crunch' in host:
    fileStr = 'ocerean'
# Subset input file dependent on host
hostFiles = []
for count,inFile in enumerate(inFiles):
    if fileStr in inFile:
        hostFiles += [inFile]
# Reset input to sublist
inFiles = hostFiles
del(hostFiles)

#%% Loop through infiles and create grid
for count,inFile in enumerate(inFiles):
    print count,inFile
    fH = cdm.open(inFile)
    if 'IPRC' in inFile:
        tmp = fH('thetao',time=slice(0,1),level=slice(0,1))(squeeze=1) ; # Load grid only
        tmp1 = fH('thetao',time=slice(0,1),level=slice(0,1))
    elif 'Ishii' in inFile:
        tmp = fH('thetao',time=slice(0,1),lev=slice(0,1))(squeeze=1) ; # Load grid only
        tmp1 = fH('thetao',time=slice(0,1),lev=slice(0,1))
    elif 'JAMSTEC' in inFile or 'UCSD' in inFile:
        tmp = fH('thetao',time=slice(0,1),pressure=slice(0,1))(squeeze=1) ; # Load grid only
        tmp1 = fH('thetao',time=slice(0,1),pressure=slice(0,1))
    elif 'ORAS4' in inFile or 'UCSD' in inFile:
        tmp = fH('thetao',time=slice(0,1),deptht=slice(0,1))(squeeze=1) ; # Load grid only
        tmp1 = fH('thetao',time=slice(0,1),deptht=slice(0,1))
    else:
        tmp = fH('thetao',time=slice(0,1),depth=slice(0,1))(squeeze=1) ; # Load grid only
        tmp1 = fH('thetao',time=slice(0,1),depth=slice(0,1))
    areacello = gnu.area_weights(tmp)
    print('areacelloM2 global:',areacello.sum())
    print('shape(tmp):',tmp.shape)
    print('shape(areacello):',areacello.shape)
    print('earthOceanAreaM2: ',(areacello*earthAreaM2).sum())
    print('Check oceanAreaM2:',earthWaterAreaKm2*1e6)
    #areacello = areacello.squeeze() ; # Loses grid info
    areacello = areacello*earthAreaM2 ; # Multiply by Earth area
    areacello.id = 'areacello'
    areacello.units = 'm2'
    areacello.earthSurfaceAreaM2  = earthAreaM2
    areacello.earthWaterAreaM2    = earthWaterAreaKm2*1e6
    areacello.earthLandAreaM2     = earthLandAreaKm2*1e6
    areacello.oceanSurfaceAreaM2  = areacello.sum()
    # Redress variable
    print('shape(areacello) 0:',areacello.shape)
    #areacello = areacello[:,mv.newaxis]
    areacello = np.expand_dims(areacello,axis=0)
    print('shape(areacello) 1:',areacello.shape)
    areacello = np.expand_dims(areacello,axis=0)
    print('shape(areacello) 2:',areacello.shape)
    areacello.setAxis(0,tmp1.getAxis(0))
    areacello.setAxis(1,tmp1.getAxis(1))
    del(tmp,tmp1) ; gc.collect()
    #continue
    # Create outfile handle
    outFile = replace(replace(inFile,'thetao','areacello'),'.xml','.nc')
    if os.path.isfile(outFile):
        print 'purge'
        os.remove(outFile)
    fHo = cdm.open(outFile,'w')
    fHo.write(areacello)
    globalAttWrite(fHo,options=None) ; # Use function to write standard global atts
    fHo.close()
    fH.close()

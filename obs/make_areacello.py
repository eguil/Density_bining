#!/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 24 11:50:23 2016

Paul J. Durack 24th April 2016

This script generates area weight (areacello) variables from input grids

PJD 24 Apr 2016     - Started
PJD 25 Apr 2016     - Corrected confused logic with IPRC - drive_obs not here
                    - TODO:

@author: durack1
"""
from string import replace
from durolib import globalAttWrite
import glob,os,socket
import cdms2 as cdm
import genutil as gnu

# Set nc classic as outputs
cdm.setNetcdfShuffleFlag(0)
cdm.setNetcdfDeflateFlag(1)
cdm.setNetcdfDeflateLevelFlag(9) ; # Set to auto
cdm.setCompressionWarnings(False)
cdm.setAutoBounds(True)

#%% Set variables
earthArea = 510.072e6 ; # Earth area in km^2
earthAreaM2 = earthArea*1e6 ; # Earth area in m^2

#%% Get input grids
inFiles = glob.glob('*thetao*')

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
        tmp = fH('thetao',time=slice(0,1),level=slice(0,1)) ; # Load grid only
    elif 'Ishii' in inFile:
        tmp = fH('thetao',time=slice(0,1),lev=slice(0,1)) ; # Load grid only
    elif 'JAMSTEC' in inFile or 'UCSD' in inFile:
        tmp = fH('thetao',time=slice(0,1),pressure=slice(0,1)) ; # Load grid only
    elif 'ORAS4' in inFile or 'UCSD' in inFile:
        tmp = fH('thetao',time=slice(0,1),deptht=slice(0,1)) ; # Load grid only
    else:
        tmp = fH('thetao',time=slice(0,1),depth=slice(0,1)) ; # Load grid only
    areacello = gnu.area_weights(tmp)
    #areacello = areacello.squeeze() ; # Loses grid info
    areacello = areacello*earthAreaM2 ; # Multiply by Earth area
    areacello.id = 'areacello'
    areacello.units = 'm2'
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

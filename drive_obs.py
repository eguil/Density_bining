#!/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 20 15:13:29 2016

Paul J. Durack 20th April 2016

This script generates input lists of obs ocean fields and drives densityBin

PJD 20 Apr 2016     - Started file
PJD 24 Apr 2016     - Updated to rewrite inputs rather than function
PJD 24 Apr 2016     - Added host check as inputs are linked and host dependent
                    - TODO:
                    - Resolve temp/to vs thetao inconsistencies
                    - Resolve pressure vs depth inconsistencies

@author: durack1
"""

from binDensity import densityBin
from string import replace
import glob,os,socket

#%% Create input list of files
outDir = '/work/durack1/Shared/160421_data_density'
inFiles = glob.glob(os.path.join(outDir,'*'))
processList = []
for count,inFile in enumerate(inFiles):
    if inFile.lower().endswith(('.nc', '.xml')) and '.density.' not in inFile:
        processList += [inFile]
processList.sort()

#%% Determine machine and host files
host = socket.gethostname()
if 'ocean' in host:
    fileStr = 'obs'
elif 'crunch' in host:
    fileStr = 'ocerean'
# Subset input file dependent on host
hostFiles = []
for count,inFile in enumerate(processList):
    if fileStr in inFile:
        hostFiles += [inFile]
# Reset input to sublist
processList = hostFiles
del(hostFiles)

#%% Call densityBin
iterCount = len(processList)/3 ; # Assumes three input files
listCount = 0
for input in range(0,iterCount):
    print 'Processing: ',input
    areacello = processList[listCount]
    so = processList[listCount+1]
    thetao = processList[listCount+2]
    outfileDensity = os.path.join(outDir,replace(areacello,'areacello','density'))
    if os.path.isfile(outfileDensity):
        os.remove(outfileDensity)
    listCount = listCount+3
    print 'areacello: ',areacello
    print 'so:        ',so
    print 'thetao:    ',thetao
    #densityBin(fileT,fileS,fileFx,'./out.nc',debug=True,timeint='all',mthout=True)
    densityBin(thetao,so,areacello,outfileDensity,debug=True,timeint='all')

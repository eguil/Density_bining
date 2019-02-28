#!/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 20 15:13:29 2016

Paul J. Durack 20th April 2016

This script generates input lists of obs ocean fields and drives densityBin

PJD 20 Apr 2016     - Started file
PJD 24 Apr 2016     - Updated to rewrite inputs rather than function
PJD 24 Apr 2016     - Added host check as inputs are linked and host dependent
PJD 25 Apr 2016     - Updated outFile name from so not areacello (xml -> nc)
PJD 25 Apr 2016     - Added IPRC fudge - needs rewriting
PJD 13 Feb 2019     - Updated to write new bugfixed data to new dir (was 160421)
                    - TODO:
                    - Rewrite: IPRC, ORAS4, SODA224
                    - Resolve temp/to vs thetao inconsistencies (JAMSTEC, UCSD, Ishii, SODA224)
                    - Resolve pressure vs depth inconsistencies (JAMSTEC)

@author: durack1
"""

import glob,os,socket
from binDensity import densityBin
from string import replace

#%% Create input list of files
outDir = '/work/durack1/Shared/190213_data_density'
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

#%% IPRC fudge
noIPRC = []
for count,inFile in enumerate(processList):
    if 'IPRC' not in inFile:
        noIPRC += [inFile]
# Reset input to sublist
processList = noIPRC
del(noIPRC)

#%% Call densityBin
iterCount = len(processList)/3 ; # Assumes three input files
listCount = 0
for input in range(0,iterCount):
    print 'Processing: ',input
    areacello       = processList[listCount]
    so              = processList[listCount+1]
    thetao          = processList[listCount+2]
    outfileDensity  = os.path.join(outDir,replace(replace(so,'so','density'),'.xml','.nc'))
    if os.path.isfile(outfileDensity):
        #print 'found file'
        #print os.path.getsize(outfileDensity),(1024*1024)
        if os.path.getsize(outfileDensity) < (1024*1024):
            # if file exists and is smaller than 1mb:
            print 'purge'
            os.remove(outfileDensity)
        else:
            print 'found big file, skipping'
            listCount = listCount+3
            continue
    listCount = listCount+3
    print 'areacello: ',areacello
    print 'so:        ',so
    print 'thetao:    ',thetao
    #densityBin(fileT,fileS,fileFx,'./out.nc',debug=True,timeint='all',mthout=True)
    #densityBin(thetao,so,areacello,outfileDensity,debug=True,timeint='all')
    densityBin(thetao,so,areacello,'none',outfileDensity,debug=True,timeint='all')

#!/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 28 15:15:04 2014

Paul J. Durack 28th September 2014

This script generates input lists of cmip5 ocean fields and drives surface_transf.py

PJD 28 Oct 2014     - Started file
PJD 22 Nov 2014     - Updated to grab all files
                    TODO:
                    - Need to consider cases where version numbers of input variables don't align
                    - Need to consider adding additional files using so/thetao fields rather than 2D sos/tos

@author: durack1
"""

import argparse,datetime,gc,glob,os,sys ; #re
from surface_transf import surfTransf
from durolib import trimModelList,writeToLog #fixVarUnits,
from string import replace
from socket import gethostname

# Purge spyder variables
if 'e' in locals():
    del(e,pi,sctypeNA,typeNA)
    gc.collect()

#%%
# Set conditional whether files are created or just numbers are calculated
parser = argparse.ArgumentParser()
parser.add_argument('modelSuite',metavar='str',type=str,help='including \'cmip3/5\' as a command line argument will select one model suite to process')
parser.add_argument('experiment',metavar='str',type=str,help='include \'experiment\' as a command line argument')
parser.add_argument('outPath',metavar='str',type=str,help='include \'outPath\' as a command line argument')
args = parser.parse_args()
# Test arguments
if (args.modelSuite in ['cmip3','cmip5']):
    modelSuite  = args.modelSuite ; # 1 = make files
else:
    print "** Invalid arguments - no *.nc files will be written **"
if (args.experiment in ['20c3m','historical','historicalNat','rcp26','rcp45','rcp60','rcp85','sresa1b','sresa2','sresb1']):
    experiment   = args.experiment
else:
    print "** Invalid arguments - no *.nc files will be written **"
if not (os.path.exists(args.outPath)):
    outPath     = os.path.join('/work/guilyardi/Shared/data_density',datetime.datetime.now().strftime("%y%m%d"));
else:
    outPath     = args.outPath
    print "** *.nc files will be written to",args.outPath," **"
   
#%%
'''
## TEST ##
modelSuite = 'cmip5'
experiment = 'historical'
outPath     = '/export/durack1/git/Density_bining/test'
#experiment = 'rcp85'
#outPath     = '/work/guilyardi/git/Density_bining/test_cmip5'
#outPath   = os.path.join('/work/durack1/Shared/data_density',datetime.datetime.now().strftime("%y%m%d"));
#modelSuite = 'cmip3'
#experiment = '20c3m'
#experiment = 'sresa2'
'''

# Create logfile
timeFormat = datetime.datetime.now().strftime("%y%m%d_%H%M%S")
logPath = '/export/durack1/git/Density_bining/'
logfile = os.path.join(logPath,"".join([timeFormat,'_drive_surface-',modelSuite,'-',experiment,'-',gethostname().split('.')[0],'.log'])) ; # WORK MODE
writeToLog(logfile,"".join(['TIME: ',timeFormat]))
writeToLog(logfile,"".join(['HOSTNAME: ',gethostname()]))
print "".join(['** Processing files from ',modelSuite,' for ',experiment,' **'])
writeToLog(logfile,"".join(['** Processing files from ',modelSuite,' for ',experiment,' **']))

#%%
# Set generic paths
soPath      = os.path.join('/work',modelSuite,experiment,'ocn/mo/sos');
thetaoPath  = os.path.join('/work',modelSuite,experiment,'ocn/mo/tos');
hfdsPath    = os.path.join('/work',modelSuite,experiment,'ocn/mo/hfds');
wfoPath     = os.path.join('/work',modelSuite,experiment,'ocn/mo/wfo');
fxPath      = os.path.join('/work',modelSuite,'fx/fx/areacello');

# Validate paths
if not os.path.exists(outPath):
    pass
    #os.makedirs(outPath)
if not os.path.exists(soPath) or not os.path.exists(thetaoPath) or not os.path.exists(hfdsPath) \
   or not os.path.exists(wfoPath) or not os.path.exists(fxPath):
    print "** Invalid source data path - no *.nc files will be written **"
    sys.exit()

# thetao
list_thetao_files = glob.glob(os.path.join(thetaoPath,'*.xml'))
list_thetao_files = trimModelList(list_thetao_files) ; # only need to match model, no version info required
list_thetao_files.sort()
# salinity
list_so_files = glob.glob(os.path.join(soPath,'*.xml'))
list_so_files = trimModelList(list_so_files) ; # only need to match model, no version info required
list_so_files.sort()
# heatflux downward shortwave
list_hfds_files = glob.glob(os.path.join(hfdsPath,'*.xml'))
list_hfds_files = trimModelList(list_hfds_files) ; # only need to match model, no version info required
list_hfds_files.sort()
# freshwater flux
list_wfo_files = glob.glob(os.path.join(wfoPath,'*.xml'))
list_wfo_files = trimModelList(list_wfo_files) ; # only need to match model, no version info required
list_wfo_files.sort()
# areacello
list_fx_files = glob.glob(os.path.join(fxPath,'*.xml'))
list_fx_files = trimModelList(list_fx_files) ; # only need to match model, no version info required
list_fx_files.sort()

# Create comparable model+ver for pairing
list_thetao = [] ; list_thetao_noVer = []
for infile in list_thetao_files:
    tmp = replace(replace(replace(infile.split('/')[-1],'.tos',''),''.join([modelSuite,'.']),''),'.latestX.xml','')
    list_thetao += [tmp];
    tmp = '.'.join(tmp.split('.')[0:-1]) ; # truncate version info
    list_thetao_noVer += [tmp]
del(infile,tmp)
list_so = [] ; list_so_noVer = []
for infile in list_so_files:
    tmp = replace(replace(replace(infile.split('/')[-1],'.sos',''),''.join([modelSuite,'.']),''),'.latestX.xml','')
    list_so += [tmp]
    tmp = '.'.join(tmp.split('.')[0:-1]) ; # truncate version info
    list_so_noVer += [tmp]
del(infile,tmp)
list_hfds = [] ; list_hfds_noVer = []
for infile in list_hfds_files:
    tmp = replace(replace(replace(infile.split('/')[-1],'.hfds',''),''.join([modelSuite,'.']),''),'.latestX.xml','')
    list_hfds += [tmp]
    tmp = '.'.join(tmp.split('.')[0:-1]) ; # truncate version info
    list_hfds_noVer += [tmp]
del(infile,tmp)
list_wfo = [] ; list_wfo_noVer = []
for infile in list_wfo_files:
    tmp = replace(replace(replace(infile.split('/')[-1],'.wfo',''),''.join([modelSuite,'.']),''),'.latestX.xml','')
    list_wfo += [tmp]
    tmp = '.'.join(tmp.split('.')[0:-1]) ; # truncate version info
    list_wfo_noVer += [tmp]
del(infile,tmp)
list_fx_model = []
for infile in list_fx_files:
    tmp = infile.split('/')[-1].split('.')[1]
    list_fx_model += [tmp]
del(infile,tmp) ; gc.collect()    

#%%
# Match hfds with wfo
list_inFiles = [[None] * 6 for i in range(len(list_hfds_files))]
del(i) ; gc.collect()
for x,model in enumerate(list_hfds):
    modelNoVer = list_hfds_noVer[x]
    list_inFiles[x][0] = list_hfds_files[x]
    # Pair hfds with wfo
    try:
        index = list_wfo_noVer.index(modelNoVer)
        list_inFiles[x][1] = list_wfo_files[index]
    except:
        print format(x,'03d'),''.join(['No wfo match for hfds: ',model])
    # Pair hfds with so
    try:
        index = list_so_noVer.index(modelNoVer)
        list_inFiles[x][2] = list_so_files[index]
    except:
        print format(x,'03d'),''.join(['No so match for hfds: ',model])
    # Pair hfds with thetao
    try:
        index = list_thetao_noVer.index(modelNoVer)
        list_inFiles[x][3] = list_thetao_files[index]
    except:
        print format(x,'03d'),''.join(['No thetao match for hfds: ',model])
    # Pair hfds with areacello
    model = modelNoVer.split('.')[0]
    try:
        index = list_fx_model.index(model)
        list_inFiles[x][4] = list_fx_files[index]
    except:
        print format(x,'03d'),''.join(['No fx match for hfds: ',model])
    # Create output fileName
    list_inFiles[x][5] = replace(replace(list_inFiles[x][0].split('/')[-1],'.hfds.','.surfTrans.'),'.latestX.xml','.nc')

#%%
# Remove blank entries
tmp = []
for x in list_inFiles:
    if None not in x and x not in tmp:
        tmp.append(x)
list_inFiles = tmp
del(tmp,x) ; gc.collect()

#%%
# Process model list
# Use for debugging
modelInd = [0] ; # Test suite to capture all errors
list_sht = []
for count,x in enumerate(list_inFiles):
    if count in modelInd:
        list_sht.append(x)
for x,model in enumerate(list_sht):
#for x,model in enumerate(list_inFiles):
    # Get steric outfile name
    outfileTransf = os.path.join(outPath,model[5])
    writeToLog(logfile,''.join(['Processing:   ',outfileTransf.split('/')[-1]]))
    print 'FileCount: ',x    
    print 'outPath:   ','/'.join(outfileTransf.split('/')[0:-1])
    print 'outfile:   ',outfileTransf.split('/')[-1]
    print 'hfds:      ',model[0].split('/')[-1]
    print 'wfo:       ',model[1].split('/')[-1]
    print 'sos:       ',model[2].split('/')[-1]
    print 'tos:       ',model[3].split('/')[-1]
    print 'areacello: ',model[4].split('/')[-1]
    # Call surfTransf
    #surfTransf(fileFx,fileTos,fileSos,fileHef,fileWfo,outFile,debug=True,timeint='all'):
    surfTransf(model[4],model[3],model[2],model[0],model[1],outfileTransf,debug=True,timeint='1,24')

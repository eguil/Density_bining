#!/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 23 13:53:04 2013

Paul J. Durack 23rd July 2013

This script generates input lists of cmip5 ocean fields and drives densit_bin

PJD 14 Sep 2014     - Started file
                    - TODO:

@author: durack1
"""

import argparse,datetime,gc,glob,os,sys ; #re
import cdms2 as cdm
from durolib import fixVarUnits,trimModelList,writeToLog
from makeStericLib import makeSteric
#from socket import gethostname
from string import replace

# Set netcdf file criterion - turned on from default 0s
#cdm.setCompressionWarnings(0) ; # Suppress warnings
#cdm.setNetcdfShuffleFlag(0)
#cdm.setNetcdfDeflateFlag(1)
#cdm.setNetcdfDeflateLevelFlag(9)
# Hi compression: 1.4Gb file ; # Single salt variable
# No compression: 5.6Gb ; Standard (compression/shuffling): 1.5Gb ; Hi compression w/ shuffling: 1.5Gb
#cdm.setAutoBounds(1) ; # Ensure bounds on time and depth axes are generated

# Purge spyder variables
if 'e' in locals():
    del(e,pi,sctypeNA,typeNA)
    gc.collect()

# Get time format
timeFormat = datetime.datetime.now().strftime("%y%m%d")

# Set conditional whether files are created or just numbers are calculated
parser = argparse.ArgumentParser()
parser.add_argument('modelSuite',metavar='str',type=str,help='including \'cmip3/5\' as a command line argument will select one model suite to process')
parser.add_argument('experiment',metavar='str',type=str,help='include \'experiment\' as a command line argument')
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

#%%
# Create logfile
#timeFormat = datetime.datetime.now().strftime("%y%m%d_%H%M%S")
#logfile = os.path.join(replace(outPath,out_path_ext,''),"".join([time_format,'_drive_models-',model_suite,'-',experiment,'-',time_str,'-',drift,'-',gethostname().split('.')[0],'.log'])) ; # WORK MODE
#writeToLog(logfile,"".join(['TIME: ',time_format]))
#writeToLog(logfile,"".join(['HOSTNAME: ',gethostname()]))
#print "".join(['** Processing files from ',model_suite,' for ',experiment,': ',drift,' **'])
#writeToLog(logfile,"".join(['** Processing files from ',model_suite,' for ',experiment,': ',drift,' **']))
#%%

## TEST ## 
'''
modelSuite = 'cmip5'
experiment = 'historical'
experiment = 'rcp85'

modelSuite = 'cmip3'
experiment = '20c3m'
experiment = 'sresa2'
'''
#%%
# Set generic paths
outPath    = os.path.join('/work/durack1/Shared/data_density',timeFormat);
soPath     = os.path.join('/work',modelSuite,experiment,'ocn/mo/so');
thetaoPath = os.path.join('/work',modelSuite,experiment,'ocn/mo/thetao');
fxPath     = os.path.join('/work',modelSuite,'fx/fx/areacello');

# Validate paths
if not os.path.exists(outPath):
    pass
    #os.makedirs(outPath

if not os.path.exists(soPath) or not os.path.exists(thetaoPath) or not os.path.exists(fxPath):
    print "** Invalid source data path - no *.nc files will be written **"
    sys.exit()

# thetao
list_thetao_files = glob.glob(os.path.join(thetaoPath,'*.xml'))
list_thetao_files.sort()
# salinity
list_so_files = glob.glob(os.path.join(soPath,'*.xml'))
list_so_files.sort()
# areacello
list_fx_files = glob.glob(os.path.join(fxPath,'*.xml'))
list_fx_files = trimModelList(list_fx_files) ; # only need to match model, no version info required
list_fx_files.sort()

# Create comparable model+ver for pairing
list_thetao = [] ; list_thetao_noVer = []
for infile in list_thetao_files:
    tmp = replace(replace(replace(infile.split('/')[-1],'.thetao',''),''.join([modelSuite,'.']),''),'.latestX.xml','')
    list_thetao += [tmp];
    tmp = '.'.join(tmp.split('.')[0:-1]) ; # truncate version info
    list_thetao_noVer += [tmp]
del(infile,tmp)

list_so = []
for infile in list_so_files:
    tmp = replace(replace(replace(infile.split('/')[-1],'.so',''),''.join([modelSuite,'.']),''),'.latestX.xml','')
    list_so += [tmp]
del(infile,tmp)

# Match so with thetao
list_soAndthetao = [[None] * 5 for i in range(len(list_so_files))]
del(i) ; gc.collect()
for x,model in enumerate(list_so):
    modelNoVer = '.'.join(model.split('.')[0:-1]) ; # so file
    try:
        index = list_thetao_noVer.index(modelNoVer)
        list_soAndthetao[x][0] = x ; # Write so index
        list_soAndthetao[x][1] = list_so_files[x] ; # Write so file
        list_soAndthetao[x][2] = index ; # Write thetao index
        list_soAndthetao[x][3] = list_thetao_files[index] ; # Write thetao file
        list_soAndthetao[x][4] = ''.join([''.join([modelSuite,'.']),replace(model,'.ocn.Omon','.ocn.Omon.density'),'.nc']) ; # Write output filename
    except:
        print x,''.join(['No so match for thetao: ',model])
    # Test for version inconsistency
    matching = [s for s in list_thetao if model in s]
    if not matching and list_soAndthetao[x] != [None, None, None, None, None]:
        print ''.join(['** Version clash - so: ',model]) ; #,' thetao: ',','.join(matching)])
        #writeToLog(logfile,''.join(['** Version clash - so: ',model]))
        print ''.join(['so    : ',list_soAndthetao[x][1].split('/')[-1]])
        #writeToLog(logfile,''.join(['so    : ',list_soAndthetaoAndfx[x][1].split('/')[-1]]))
        print ''.join(['thetao: ',list_soAndthetao[x][3].split('/')[-1]])
        #writeToLog(logfile,''.join(['so    : ',list_soAndthetaoAndfx[x][3].split('/')[-1]]))
del(model,x,index,list_so,list_so_files,list_thetao,list_thetao_files,list_thetao_noVer,matching,modelNoVer,s) ; gc.collect()

# Remove blank entries - First remove duplicate None*5 fields then remove final
tmp = []
for x in list_soAndthetao:
    if x not in tmp:
        tmp.append(x)
try:
    ind = tmp.index([None, None, None, None, None])
    del tmp[ind]
    del(ind)
except Exception,err:
    print 'Exception thrown: ',err
list_soAndthetao = tmp
del(tmp,x) ; gc.collect()

# Trim out dupes using versioning info
if 'cmip5' in modelSuite:
    list_so = [tmp[1] for tmp in list_soAndthetao]
    list_so = trimModelList(list_so) ; # Trim out dupes using versioning info
    # Use trim list to remove dupes
    list_soAndthetao_noDupes = []
    for count,pairs in enumerate(list_soAndthetao):
        if pairs[1] in list_so:
            list_soAndthetao_noDupes += [pairs]
    del(list_so,count,pairs,tmp); gc.collect()
    list_soAndthetao = list_soAndthetao_noDupes
    del(list_soAndthetao_noDupes) ; gc.collect()

# Add in fx file
list_soAndthetaoAndfx = [[None] * 6 for i in range(len(list_soAndthetao))]
for count,pairs in enumerate(list_soAndthetao):
    for x,val in enumerate(pairs):
        list_soAndthetaoAndfx[count][x] = val
    del(x,val)
    # Match with fx and fill
    model = list_soAndthetaoAndfx[0][1].split('.')[1]
    for count2,fx in enumerate(list_fx_files):
        model_fx = fx.split('.')[1]
        if model == model_fx:
            list_soAndthetaoAndfx[count][5] = list_fx_files[count2]
            continue
del(i,fx,list_fx_files,list_soAndthetao,count,count2,pairs,model,model_fx) ; gc.collect()
#%%

# Process model list
for x,model in enumerate(list_soAndthetaoAndfx):
    # Get steric outfile name
    outfileDensity = os.path.join(outPath,list_soAndthetaoAndfx[x][4])
    print "".join(['Processing:   ',outfileDensity.split('/')[-1]])
    #writeToLog(logfile,"".join(['Processing:   ',outfileDensity.split('/')[-1]]))

    # Start steric file creation
    #make_steric(salinity,salinityChg,temp,tempChg,outFileName,thetao,pressure)   
    makeSteric(so,so_chg,thetao,thetao_chg,outfile_steric,True,pressure)

'''
# Check code for input variables
    # Load z-axis from salt file to check validity
    f_h = cdm.open(os.path.join(soPath,list_soAndthetaoAndfx[x][1]))
    # Check for pressure or depth
    if 'cmip5' in modelSuite:
        z_coord = f_h.getAxis('lev')
    elif 'cmip3' in modelSuite:
        z_coord = f_h.getAxis('depth')
    if z_coord.units == 'dbar':
        pressure = True
    elif z_coord.units == 'm':
        pressure = False
    elif z_coord.units == '':
        
        # DEAL WITH SIGMA LEVEL MODELS        
        print "".join(['** infile: ',list_soAndthetaoAndfx[x][1],' has non-recognised depth coord, skipping.. **','\n'])
        #writeToLog(logfile,"".join(['** infile: ',list_soAndthetaoAndfx[x][1],' has non-recognised depth coord, skipping.. **','\n']))
        continue

    # Load salt variables
    so      = f_h('so',squeeze=1)
    f_h.close()
    # Correct so to PSS-78
    [so,so_fixed] = fixVarUnits(so,'so',True) #,logfile)
    # Validate variable ranges
    if so.min() < -1 or so.max() > 150:
        print "".join(['** infile: ',list_soAndthetaoAndfx[x][1],' has invalid data - max: ','{:08.2f}'.format(float(so.max())),' min: ','{:08.2f}'.format(float(so.min())),' skipping.. **','\n'])
        #writeToLog(logfile,"".join(['** infile: ',list_soAndthetaoAndfx[x][1],' has invalid data - max: ','{:08.2f}'.format(float(so.max())),' min: ','{:08.2f}'.format(float(so.min())),' skipping.. **','\n']))
        continue
    
    # Load thetao variables
    f_h         = cdm.open(os.path.join(thetaoPath,list_soAndthetaoAndfx[x][3]))
    thetao      = f_h('thetao',squeeze=1)
    f_h.close()
    # Correct K to degrees_C
    [thetao,thetao_fixed] = fixVarUnits(thetao,'thetao',True) #,logfile)
    # Validate variable ranges
    if thetao.min() < -10 or thetao.max() > 50:
        print "".join(['** infile: ',list_soAndthetaoAndfx[x][3],' has invalid data - max: ','{:08.2f}'.format(float(thetao.max())),' min: ','{:08.2f}'.format(float(thetao.min())),' skipping.. **','\n'])
        writeToLog(logfile,"".join(['** infile: ',list_soAndthetaoAndfx[x][3],' has invalid data - max: ','{:08.2f}'.format(float(thetao.max())),' min: ','{:08.2f}'.format(float(thetao.min())),' skipping.. **','\n']))
        continue
    
    # HadGEM2-AO.historical.r1i1p1 reverse axis
    if 'HadGEM2-AO' in list_soAndthetaoAndfx[x][1]:
        print "Dealing with HadGEM2-AO inverted z-axis.."
        #writeToLog(logfile,"Dealing with HadGEM2-AO inverted z-axis..")
        f_h             = cdm.open(os.path.join(soPath,list_soAndthetaoAndfx[x][1]))
        z_coord         = f_h.getAxis('lev')
        z_coord2        = cdm.createAxis(z_coord[::-1])
        f_h.close()
        z_coord2.id     = 'lev'
        for k in z_coord.attributes.keys():
            setattr(z_coord2,k,z_coord.attributes[k])
        z_coord         = z_coord2
        del(z_coord2,k) ; gc.collect()
        so              = so[::-1,...] ; # Syntax so[39::-1] & so[::-1] also works
        so.setAxis(0,z_coord)
        so_chg          = so_chg[::-1,...]
        so_chg.setAxis(0,z_coord)
        thetao          = thetao[::-1,...]
        thetao.setAxis(0,z_coord)
        thetao_chg      = thetao_chg[::-1,...]
        thetao_chg.setAxis(0,z_coord)
'''
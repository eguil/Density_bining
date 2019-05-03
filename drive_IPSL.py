#!/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 17 17:40:42 2014

Paul J. Durack 17th September 2014

This script drives densityBin with a single model input

PJD 17 Sep 2014     - Started file
PJD 18 Sep 2014     - Updated output directory
PJD  1 Oct 2014     - Updated using timeint
                    - TODO:

@author: durack1
"""

from binDensity import densityBin

modelThetao = '/prodigfs/project/CMIP5/main/IPSL/IPSL-CM5B-LR/piControl/mon/ocean/Omon/r1i1p1/latest/thetao/thetao_Omon_IPSL-CM5B-LR_piControl_r1i1p1_183001-187912.nc'
modelSo = '/prodigfs/project/CMIP5/main/IPSL/IPSL-CM5B-LR/piControl/mon/ocean/Omon/r1i1p1/latest/so/so_Omon_IPSL-CM5B-LR_piControl_r1i1p1_183001-187912.nc'
modelVo = '/prodigfs/project/CMIP5/main/IPSL/IPSL-CM5B-LR/piControl/mon/ocean/Omon/r1i1p1/latest/vo/vo_Omon_IPSL-CM5B-LR_piControl_r1i1p1_183001-187912.nc'

modelAreacello = '/prodigfs/project/CMIP5/main/IPSL/IPSL-CM5A-LR/piControl/fx/ocean/fx/r0i0p0/latest/areacello/areacello_fx_IPSL-CM5A-LR_piControl_r0i0p0.nc'
outfileDensity = '/home/ysilvy/Density_bining/test/cmip5.IPSL-VLR0.piControl.rip.mon.ocean.Omon.density.nc'

grid_T_file = '/prodigfs/project/CMIP5/main/IPSL/IPSL-CM5B-LR/piControl/mon/ocean/Omon/r1i1p1/latest/thetao/thetao_Omon_IPSL-CM5B-LR_piControl_r1i1p1_183001-187912.nc'
grid_S_file = '/prodigfs/project/CMIP5/main/IPSL/IPSL-CM5B-LR/piControl/mon/ocean/Omon/r1i1p1/latest/so/so_Omon_IPSL-CM5B-LR_piControl_r1i1p1_183001-187912.nc'
grid_V_file = '/prodigfs/project/CMIP5/main/IPSL/IPSL-CM5B-LR/piControl/mon/ocean/Omon/r1i1p1/latest/vo/vo_Omon_IPSL-CM5B-LR_piControl_r1i1p1_183001-187912.nc'
targetGridFile = '/home/ericglod/Density_bining/170224_WOD13_masks.nc'

print 'outfile:   ',outfileDensity
print 'so:        ',modelSo
print 'thetao:    ',modelThetao
print 'areacello: ',modelAreacello
# Call densityBin
#densityBin(modelThetao,modelSo,modelAreacello,outfileDensity, debug=False)
#densityBin(modelThetao,modelSo,modelAreacello,outfileDensity,timeint='1,24')
print grid_T_file

#densityBin(modelThetao,modelSo,modelAreacello,fileV=modelVo,outFile=outfileDensity,timeint='121,12', gridfT=grid_T_file, gridfS=grid_S_file, gridfV=grid_V_file)
densityBin(modelThetao,modelSo,modelAreacello,targetGrid=targetGridFile,fileV='none',outFile=outfileDensity,timeint='1,12', gridfT=grid_T_file, gridfS=grid_S_file, gridfV=grid_V_file)

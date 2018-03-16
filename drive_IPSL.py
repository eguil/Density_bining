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

modelSo = '/work/cmip5/historical/ocn/mo/so/cmip5.IPSL-CM5A-LR.historical.r1i1p1.mo.ocn.Omon.so.ver-v20111119.latestX.xml'
modelThetao = '/work/cmip5/historical/ocn/mo/thetao/cmip5.IPSL-CM5A-LR.historical.r1i1p1.mo.ocn.Omon.thetao.ver-v20111119.latestX.xml'
modelAreacello = '/work/cmip5/fx/fx/areacello/cmip5.IPSL-CM5A-LR.historical.r0i0p0.fx.ocn.fx.areacello.ver-1.latestX.xml'
outfileDensity = 'test/cmip5.IPSL-CM5A-LR.historical.r1i1p1.mo.ocn.Omon.density.ver-v20111119-compressed.nc'

grid_T_file = '/prodigfs/project/CMIP5/main/IPSL/IPSL-CM5B-LR/piControl/mon/ocean/Omon/r1i1p1/latest/thetao/thetao_Omon_IPSL-CM5B-LR_piControl_r1i1p1_183001-187912.nc'
grid_S_file = '/prodigfs/project/CMIP5/main/IPSL/IPSL-CM5B-LR/piControl/mon/ocean/Omon/r1i1p1/latest/so/so_Omon_IPSL-CM5B-LR_piControl_r1i1p1_183001-187912.nc'

print 'outfile:   ',outfileDensity
print 'so:        ',modelSo
print 'thetao:    ',modelThetao
print 'areacello: ',modelAreacello
# Call densityBin
#densityBin(modelThetao,modelSo,modelAreacello,outfileDensity, debug=False)
#densityBin(modelThetao,modelSo,modelAreacello,outfileDensity,timeint='1,24')
densityBin(modelThetao,modelSo,modelAreacello,outfileDensity,timeint='1,12', gridfT = grid_T_file, gridfS = grid_S_file)

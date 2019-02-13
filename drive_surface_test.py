#!/bin/env python
# -*- coding: utf-8 -*-
#
# Test surface transformation on ciclad
#
# Oct 2017 EG
#
# ----------------------------------------
#
from surface_transf import surfTransf
indir ='/prodigfs/project/CMIP5/main/CNRM-CERFACS/CNRM-CM5-2/piControl/mon/ocean/Omon/r1i1p1/latest'
# we need: tos, sos, hfds, wfo, fixed
fileFx  ='/prodigfs/project/CMIP5/main/CNRM-CERFACS/CNRM-CM5-2/piControl/fx/ocean/fx/r0i0p0/latest/areacello/areacello_fx_CNRM-CM5-2_piControl_r0i0p0.nc'
fileTos = indir+'/tos/tos_Omon_CNRM-CM5-2_piControl_r1i1p1_185001-185912.nc'
fileSos = indir+'/sos/sos_Omon_CNRM-CM5-2_piControl_r1i1p1_185001-185912.nc'
fileHef = indir+'/hfds/hfds_Omon_CNRM-CM5-2_piControl_r1i1p1_185001-185912.nc'
fileWfo = indir+'/wfo/wfo_Omon_CNRM-CM5-2_piControl_r1i1p1_185001-185912.nc'
varnames =['tos','sos','hfds','wfo']
outFile = 'surf_transf_test.nc'

surfTransf(fileFx,fileTos,fileSos,fileHef,fileWfo,varnames,outFile,debug=True,timeint='all')

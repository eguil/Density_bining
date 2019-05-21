# -*- coding: utf-8 -*-

"""
Drive readzfile function
"""

from functions_z_analysis import readzfile
import glob

work = 'rcp85' # 'historical' 'historicalNat' 'rcp85'

if work == 'historical':
    # Historical
    indir = '/bdd/CMIP5/main/IPSL/IPSL-CM5A-LR/historical/mon/ocean/Omon/r2i1p1/latest/'
    list_so = sorted(glob.glob(indir+'so/*.nc'))
    list_thetao = sorted(glob.glob(indir+'thetao/*.nc'))
    print('Historical files :', list_so, list_thetao)
    compute_gamma = True

if work == 'historicalNat':
    # HistoricalNat
    indir = '/bdd/CMIP5/main/IPSL/IPSL-CM5A-LR/historicalNat/mon/ocean/Omon/r2i1p1/latest/'
    list_so = sorted(glob.glob(indir+'so/*.nc'))
    list_thetao = sorted(glob.glob(indir+'thetao/*.nc'))
    print('HistoricalNat files :', list_so, list_thetao)
    compute_gamma = False

if work == 'rcp85':
    # RCP8.5
    indir = '/bdd/CMIP5/main/IPSL/IPSL-CM5A-LR/rcp85/mon/ocean/Omon/r2i1p1/latest/'
    list_so = sorted(glob.glob(indir+'so/*.nc'))
    list_thetao = sorted(glob.glob(indir+'thetao/*.nc'))
    print('HistoricalNat files :', list_so, list_thetao)
    compute_gamma = True
    
# Input arguments
targetGridFile = '/home/ericglod/Density_bining/170224_WOD13_masks.nc'
outdir = '/data/ysilvy/CMIP5_annual/'
# fileT = indir+'thetao/thetao_Omon_IPSL-CM5A-LR_historical_r2i1p1_185001-189912.nc'
# fileS = indir+'so/so_Omon_IPSL-CM5A-LR_historical_r2i1p1_185001-189912.nc'
for i in range(len(list_so)):
    fileT = list_thetao[i]
    fileS = list_so[i]

    readzfile(fileT,fileS,compute_gamma,targetGridFile,outdir)

#!/usr/local/uvcdat/latest/bin/cdat
#
#
# -------------------------------------------------------------------------------------------------
# Procedure to compute water mass transformation from surface buoyancy fluxes in density space
#
#  Input fields:
#    - sst, sss, E-P, Qnet (3D, time,j,i)
#    - density grid sigrid, and delta_s (1D)
#    - target grid, including ocean basins
#  Output fields (on target grid):
#    - density flux (total, heat, fresh water) (2D rho,time, per basin or specific region)
#    - transformation (2D rho, time, per basin or specific region)
#
# Following Walin (1982) and Speer and Tziperman (1992)
# -------------------------------------------------------------------------------------------------
#   E. Guilyardi Sept 2014
#

import cdms2 as cdm
import MV2 as mv
import os, sys
import argparse
import string
import numpy as npy
import numpy.ma as ma
import cdutil as cdu
from genutil import statistics
#import support_density as sd
from binDensity import maskVal
from binDensity import eosNeutral
from binDensity import rhonGrid
from binDensity import computeArea
import time as timc
import timeit
import resource
import ESMP
from cdms2 import CdmsRegrid
from durolib import fixVarUnits
from string import replace
#
# inits
# -----
#
def alpha (t, s):
    # compute alpha=-1/rho (d rho / d T)
    dt = 0.05
    siga = eso_neutral(t, s)
    sigb = eso_neutral(t+0.05, s)
    alpha = -0.001*(sigb-siga)/dt/(1.+1.e-3*siga)
    return alpha
def betar (t, s):
    # compute beta= 1/rho (d rho / d S)
    ds = 0.01
    siga = eos_neutral(t, s)-1000.
    sigb = eos_neutral(t, s+ds)-1000.    
    beta = 0.001*(sigb-siga)/ds/(1.+1.e-3*siga)
    return beta
def cpsw (t, s, p):
    # Specific heat of sea water (J/KG C)
    CP1 = 0.
    CP2 = 0.
    SR=SQRT(ABS(S))
    # SPECIFIC HEAT CP0 FOR P=0 (MILLERO ET AL. 1973)
    A = (-1.38E-3*T+0.10727)*T-7.644
    B = (5.35E-5*T-4.08E-3)*T+0.177
    C = (((2.093236E-5*T-2.654387E-3)*T+0.1412855)*T-3.720283)*T+4217.4
    CP0 = (B*SR + A) * S + C
    # CP1 PRESSURE AND TEMPERATURE TERMS FOR S = 0
    A = (((1.7168E-8*T+2.0357E-6)*T-3.13885E-4)*T+1.45747E-2)*T-0.49592
    B = (((2.2956E-11*T-4.0027E-9)*T+2.87533E-7)*T-1.08645E-5)*T+2.4931E-4
    C = ((6.136E-13*T-6.5637E-11)*T+2.6380E-9)*T-5.422E-8
    CP1 = ((C*P+B)*P+A)*P
    # CP2 PRESSURE AND TEMPERATURE TERMS FOR S > 0
    A = (((-2.9179E-10*T+2.5941E-8)*T+9.802E-7)*T-1.28315E-4)*T+4.9247E-3
    B = (3.122E-8*T-1.517E-6)*T-1.2331E-4
    A = (A+B*SR)*S
    B = ((1.8448E-11*T-2.3905E-9)*T+1.17054E-7)*T-2.9558E-6
    B = (B+9.971E-8*SR)*S
    C = (3.513E-13*T-1.7682E-11)*T+5.540E-10
    C = (C-1.4300E-12*T*SR)*S
    CP2 = ((C*P+B)*P+A)*P
    cp = CP0 + CP1 + CP2
    return cp    

def surface_transf(fileFx, fileSst, fileSss, fileHef, fileWfo, outFile, debug=True,timeint='all'):
    # Keep track of time (CPU and elapsed)
    cpu0 = timc.clock()
    #
    # netCDF compression (use 0 for netCDF3)
    comp = 1
    cdm.setNetcdfShuffleFlag(comp)
    cdm.setNetcdfDeflateFlag(comp)
    cdm.setNetcdfDeflateLevelFlag(comp)
    cdm.setAutoBounds('on')
    # 
    # == Inits
    #
    npy.set_printoptions(precision = 2)
    # Determine file name from inputs
    modeln = fileSst.split('/')[-1].split('.')[1]
    #
    if debug >= '1':
        print ' Debug - File names:'
        print '    ', file_sst
        print '    ', file_sss
        debugp = True
    else:
        debugp = False
    #
    # Open files
    fsst  = cdm.open(fileSst)
    fsss  = cdm.open(fileSss)
    fhef  = cdm.open(fileHef)
    fwfo  = cdm.open(fileWfo)
    timeax = fsst.getAxis('time')
    #
    # Dates to read
    if timeint == 'all':
        tmin = 0
        tmax = timeax.shape[0]
    else:
        tmin = int(timeint.split(',')[0]) - 1
        tmax = tmin + int(timeint.split(',')[1])

    if debugp:
        print; print ' Debug mode'
 
    # Read file attributes to carry on to output files
    list_file   = fsst.attributes.keys()
    file_dic    = {}
    for i in range(0,len(list_file)):
        file_dic[i]=list_file[i],fsst.attributes[list_file[i] ]
    #
    # Read data
    if debugp:
        print' Read sst, sss',tmin,tmax
    sst = fsst('tos' , time = slice(tmin,tmax))
    sss = fsss('sos' , time = slice(tmin,tmax))
    if debugp:
        print' Read hfds, wfo'
    qnet = fhef('hfds', time = slice(tmin,tmax))
    emp  = fwfo('wfo' , time = slice(tmin,tmax))
    #
    # Read masking value
    valmask = sst._FillValue
    #
    # Read time and grid
    lon  = sst.getLongitude()
    lat  = sst.getLatitude()
    ingrid = sst.getGrid()
    #
    # Read cell area
    ff = cdm.open(fileFx)
    area = ff('areacello')
    ff.close()
    #
    # Define dimensions
    N_i = int(lon.shape[1])
    N_j = int(lon.shape[0])
    #
    # Define sigma grid 
    rho_min = 19
    rho_int = 26
    rho_max = 28.5
    del_s1  = 0.2
    del_s2  = 0.1
    s_s, s_sax, del_s, N_s = rhonGrid(rho_min, rho_int, rho_max, del_s1, del_s2)
    print
    print ' ==> model:', modeln
    #
    # File output inits
    #
    s_axis = cdm.createAxis(s_sax, id = 'rhon')
    s_axis.long_name = 'Neutral density'
    s_axis.units = ''
    s_axis.designateLevel()
    #
    # Monthly transformation
    outFile = replace(outFile,'.mo.','.an.')
    if os.path.exists(outFile):
        os.remove(outFile)
    outFile_f = cdm.open(outFile,'w')
    # Define dimensions
    N_i = int(sst.shape[2])
    N_j = int(sst.shape[1])
    N_t = int(sst.shape[0])
    N_s = len(sigrid)
    # Read masking value
    valmask = sst._FillValue
    # reorganise i,j dims in single dimension data
    sst  = npy.reshape(sst, (N_t, N_i*N_j))
    sss  = npy.reshape(sss, (N_t, N_i*N_j))
    emp  = npy.reshape(emp, (N_t, N_i*N_j))
    qnet = npy.reshape(qnet, (N_t, N_i*N_j))
    area = npy.reshape(area, (N_i*N_j))
    # Test variable units
    [sss,sssFixed] = fixVarUnits(sss,'sss',True)#,'logfile.txt')
    if sosFixed:
        print '     sss: units corrected'
    [sst,sstFixed] = fixVarUnits(sst,'sst',True)#,'logfile.txt')
    if sstFixed:
        print '     sst: units corrected'        
    # Physical inits
    P = 0          # surface pressure
    conwf = 1.e-3  # kg/m2/s=mm/s -> m/s
    # find non-masked points
    maskin = mv.masked_values(sst.data[0], valmask).mask 
    nomask = npy.equal(maskin,0)
    # init arrays
    areabin = npy.ones((N_t,N_s+1))*valmask # surface of bin
    denflxh = npy.ones((N_t,N_s+1))*valmask # heat flux contrib
    denflxw = npy.ones((N_t,N_s+1))*valmask # E-P contrib
    t_heat  = npy.ones((N_t))*valmask # integral heat flux
    t_wafl  = npy.ones((N_t))*valmask # integral E-P
    #
    # target horizonal grid for interp 
    fileg = '/work/guilyardi/Density_bining/WOD13_masks.nc'
    gt = cdm.open(fileg)
    maskg = gt('basinmask3')
    outgrid = maskg.getGrid()
    # global mask
    maski = maskg.mask
    # regional masks
    maskAtl = maski*1 ; maskAtl[...] = True
    idxa = npy.argwhere(maskg == 1).transpose()
    maskAtl[idxa[0],idxa[1]] = False
    maskPac = maski*1 ; maskPac[...] = True
    idxp = npy.argwhere(maskg == 2).transpose()
    maskPac[idxp[0],idxp[1]] = False
    maskInd = maski*1 ; maskInd[...] = True
    idxi = npy.argwhere(maskg == 3).transpose()
    maskInd[idxi[0],idxi[1]] = False
    masks = [maski, maskAtl, maskPac, maskInd]
    loni    = maskg.getLongitude()
    lati    = maskg.getLatitude()
    Nii     = int(loni.shape[0])
    Nji     = int(lati.shape[0])
    #
    # Compute area of target grid and zonal sums
    areai = computeArea(loni[:], lati[:])
    gt.close()
    # Interpolation init (regrid)
    ESMP.ESMP_Initialize()
    regridObj = CdmsRegrid(ingrid, outgrid, denflxh.dtype, missing = valmask, regridMethod = 'linear', regridTool = 'esmf')
    # init integration intervals
    dt   = 1./float(N_t) 

    # Bin on density grid
    for t in range(N_t):
        sstt = sst[t,:]
        ssst = sss[t,:]
        # Compute density
        rhon = eosNeutral(sstt, ssst)
        # Compute buoyancy flux as mass fluxes in kg/m2/s (SI unts)
        fheat = (-alpha(sstt,ssst)/cpsw(sstt,ssst,P))*qnet[t,:]
        fwafl = (rhon+1000.)*beta(sstt,ssst)*ssst*emp[t,:]*convwf
        # bining loop
        for ks in range(N_s):
            # find indices of points in density bin
            idxbin = npy.argwhere( (rhon[t,:] >= sigrid[ks]) & (rhon[t,:] < sigrid[ks+1]) )
            denflxh[t,ks] = cdu.averager(fheat[t, idxbin] * area[idxbin], axis=1, action='sum')
            denflxw[t,ks] = cdu.averager(fwafl[t, idxbin] * area[idxbin], axis=1, action='sum')
            areabin[t,ks] = cdu.averager(area[idxbin], axis=1, action='sum')
        # last bin
        idxbin = npy.argwhere( (rhon[t,:] >= sigrid[N_s]))
        denflxh[t,N_s] = cdu.averager(fheat[t, idxbin] * area[idxbin], axis=1, action='sum')
        denflxw[t,N_s] = cdu.averager(fwafl[t, idxbin] * area[idxbin], axis=1, action='sum')
        areabin[t,N_s] = cdu.averager(area[idxbin], axis=1, action='sum')
        # Total density flux
        denflx[t,:] = denflxh[y,:] + denflxw[t,:]
        # Transformation
        #transfh[t,:] = ...
        # domain integrals
        # heat flux (conv W -> PW)
        convt  = 1.e-15
        t_heat[t] = cdu.averager(fheat*area, action='sum')*dt*convt
        # fw flux (conv mm -> m and m3/s to Sv)
        convw = 1.e-3*1.e-6
        t_wafl[t] = cdu.averager(fwafl*area, action='sum')*dt*convw
      
        
    #+ create a basins variables (loop on n masks)
    # CPU use
    print
    print ' CPU use', timc.clock() - cpu0

    print ' -> Calculated  denflx, denflxh, denflxw, t_heat, t_wafl'
    if debugp:
        print ' t_heat',t_heat
        print ' t_wafl',t_wafl

#
# <------------------------------------------------>
#       Test driver for surface transformation
# <------------------------------------------------>
#
#
# IPSL-CM5A-LR
#
file_fx = '/work/cmip5/fx/fx/areacello/cmip5.MPI-ESM-LR.historical.r0i0p0.fx.ocn.fx.areacello.ver-v20111006.latestX.xml'
file_sst = '/work/cmip5/historical/ocn/mo/tos/cmip5.MPI-ESM-LR.historical.r1i1p1.mo.ocn.Omon.tos.ver-1.latestX.xml'
file_sss = '/work/cmip5/historical/ocn/mo/sos/cmip5.MPI-ESM-LR.historical.r1i1p1.mo.ocn.Omon.sos.ver-1.latestX.xml'
file_hef = '/work/cmip5/historical/ocn/mo/hfds/cmip5.MPI-ESM-LR.historical.r1i1p1.mo.ocn.Omon.hfds.ver-1.latestX.xml'
file_wfo = '/work/cmip5/historical/ocn/mo/wfo/cmip5.MPI-ESM-LR.historical.r1i1p1.mo.ocn.Omon.wfo.ver-1.latestX.xml'
outfileSurfDen = 'test/cmip5.MPI-ESM-LR.historical.r1i1p1.mo.ocn.Omon.surfden.ver-ver-1.latest.nc'

#
# -----------------------------------------
#  Compute density flux and transformation
# -----------------------------------------
#
surface_transf(file_fx, file_sst, file_sss, file_hef, file_wfo, outfileSurfDen, debug=True,timeint='1,1')


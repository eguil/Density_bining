#!/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  8 09:11:11 CEST 2014

Eric Guilyardi 8 October 2014

This script computes water mass transformation from surface buoyancy fluxes in density space

  Input fields:
    - sst, sss, E-P, Qnet (3D, time,j,i)
    - density grid sigrid, and delta_s (1D)
    - target grid, including ocean basins
  Output fields (on target grid):
    - density flux (total, heat, fresh water) (2D rho,time, per basin or specific region)
    - transformation (2D rho, time, per basin or specific region)

 Following Walin (1982) and Speer and Tziperman (1992)computes density bins and indexes T,S values from the vertical z
 Uses McDougall and Jackett 2005 EOS (IDL routine provided by G. Madec)
 Uses python seawater package (https://github.com/ocefpaf/python-seawater)
---------------------------------------------------------------------------------

EG 8 Oct 2014     - Started file
                    - TODO: 
                     - Bug in integral fluxes calculations
                     - North vs. South calculation
                     - add ekman pumping bining (cf wcurl in densit)


@author: eguil
"""

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
import seawater as sw
#
# inits
# -----
#

def surfTransf(fileFx, fileTos, fileSos, fileHef, fileWfo, outFile, debug=True, timeint='all'):
    '''
    The surfTransf() function takes files and variable arguments and creates
    density bined surface transformation fields which are written to a specified outfile
    Author:    Eric Guilyardi : Eric.Guilyardi@locean-ipsl.upmc.fr
    Co-author: Paul J. Durack : pauldurack@llnl.gov : @durack1.
    
    Created on Wed Oct  8 09:15:59 CEST 2014

    Inputs:
    ------
    - fileTos(time,lat,lon)     - 3D SST array
    - fileSos(time,lat,lon)     - 3D SSS array
    - fileHef(time,lat,lon)     - 3D net surface heat flux array
    - fileWfo(time,lat,lon)     - 3D fresh water flux array
    - fileFx(lat,lon)           - 2D array containing the cell area values
    - outFile(str)              - output file with full path specified.
    - debug <optional>          - boolean value
    - timeint <optional>        - specify temporal step for binning <init_idx>,<ncount>

    Outputs:
    --------
    - netCDF file with monthly surface rhon, density fluxes, transformation (global and per basin)
    - use cdo yearmean to compute annual mean

    Usage:
    ------
    >>> from binDensity import surfTransf
    >>> surfTransf(file_fx, file_tos, file_sos, file_hef, file_wfo, ./output.nc, debug=True,timeint='all')

    Notes:
    -----
    - EG 8 Oct 2014   - Initial function write and tests ok

    '''
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
    modeln = fileTos.split('/')[-1].split('.')[1]
    #
    if debug:
        print ' Debug - File names:'
        print '    ', fileTos
        print '    ', fileSos
        debugp = True
    else:
        debugp = False
    #
    # Open files
    ftos  = cdm.open(fileTos)
    fsos  = cdm.open(fileSos)
    fhef  = cdm.open(fileHef)
    fwfo  = cdm.open(fileWfo)
    timeax = ftos.getAxis('time')
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
    list_file   = ftos.attributes.keys()
    file_dic    = {}
    for i in range(0,len(list_file)):
        file_dic[i]=list_file[i],ftos.attributes[list_file[i] ]
    #
    # Read data
    if debugp:
        print' Read tos, sos',tmin,tmax
    tos = ftos('tos' , time = slice(tmin,tmax))
    sos = fsos('sos' , time = slice(tmin,tmax))
    if debugp:
        print' Read hfds, wfo'
    qnet = fhef('hfds', time = slice(tmin,tmax))
    emp  = fwfo('wfo' , time = slice(tmin,tmax))
    tos_h = ftos['tos']
    #
    # Read masking value
    valmask = tos._FillValue
    #
    # Read time and grid
    time = tos.getTime()
    lon  = tos_h.getLongitude()
    lat  = tos_h.getLatitude()
    ingrid = tos_h.getGrid()
    #
    # Read cell area
    ff = cdm.open(fileFx)
    area = ff('areacello')
    ff.close()
    areain = area.data
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
    sigrid, s_sax, del_s, N_s = rhonGrid(rho_min, rho_int, rho_max, del_s1, del_s2)
    print
    print ' ==> model:', modeln
    #
    # File output inits
    #
    s_axis = cdm.createAxis(s_sax, id = 'rhon')
    s_axis.long_name = 'Neutral density'
    s_axis.units = 'kg m-3 (anomaly, minus 1000)'
    s_axis.designateLevel()
    #
    # Monthly transformation
    #outFile = replace(outFile,'.mo.','.an.')
    if os.path.exists(outFile):
        os.remove(outFile)
    outFile_f = cdm.open(outFile,'w')
    # Define dimensions
    N_i = int(tos.shape[2])
    N_j = int(tos.shape[1])
    N_t = int(tos.shape[0])
    print ' ==> dimensions N_t, N_j, N_i:', N_t, N_j, N_i
    # Read masking value
    valmask = tos._FillValue
    # Test variable units
    [sos,sosFixed] = fixVarUnits(sos,'sos',True)#,'logfile.txt')
    if sosFixed:
        print '     sos: units corrected'
    [tos,tosFixed] = fixVarUnits(tos,'thetao',True)#,'logfile.txt')
    if tosFixed:
        print '     tos: units corrected'        
    # Physical inits
    P = 0          # surface pressure
    # find non-masked points
    maskin = mv.masked_values(tos.data[0], valmask).mask 
    nomask = npy.equal(maskin,0)
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
    # init arrays
    tmp     = npy.ma.ones([N_t, Nji, Nii], dtype='float32')*valmask 
    denflx  = tmp.copy() # Total density flux
    denflxh = tmp.copy() # heat flux contrib
    denflxw = tmp.copy() # E-P contrib
    rhon    = tmp.copy() # surface density
    tmpi    = npy.ma.ones([Nji, Nii], dtype='float32')*valmask 
    tost    = tmpi.copy()
    sost    = tmpi.copy() 
    heft    = tmpi.copy()
    empt    = tmpi.copy()
    # Global
    atmp    = npy.ma.ones([N_t, N_s+1], dtype='float32')*valmask 
    transf  = atmp.copy() # Total tranformation
    transfh = atmp.copy() # Heat flux tranformation
    transfw = atmp.copy() # Water flux tranformation
    areabin = atmp.copy() # surface of bin
    transfh = maskVal(transfh, valmask)
    transfw = maskVal(transfw, valmask)
    transf  = maskVal(transf , valmask)
    areabin = maskVal(areabin, valmask)
    # Basin
    transfa  = atmp.copy() # Total tranformation
    transfha = atmp.copy() # Heat flux tranformation
    transfwa = atmp.copy() # Water flux tranformation
    areabina = atmp.copy() # surface of bin
    transfha = maskVal(transfha, valmask)
    transfwa = maskVal(transfwa, valmask)
    transfa  = maskVal(transfa , valmask)
    areabina = maskVal(areabina, valmask)
    #
    transfp  = atmp.copy() # Total tranformation
    transfhp = atmp.copy() # Heat flux tranformation
    transfwp = atmp.copy() # Water flux tranformation
    areabinp = atmp.copy() # surface of bin
    transfhp = maskVal(transfhp, valmask)
    transfwp = maskVal(transfwp, valmask)
    transfp  = maskVal(transfp , valmask)
    areabinp = maskVal(areabinp, valmask)
    #
    transfi  = atmp.copy() # Total tranformation
    transfhi = atmp.copy() # Heat flux tranformation
    transfwi = atmp.copy() # Water flux tranformation
    areabini = atmp.copy() # surface of bin
    transfhi = maskVal(transfhi, valmask)
    transfwi = maskVal(transfwi, valmask)
    transfi  = maskVal(transfi , valmask)
    areabini = maskVal(areabini, valmask)
    #
    tmp = npy.ma.ones((N_t))*valmask
    intHeatFlx  = tmp.copy() # integral heat flux
    intWatFlx   = tmp.copy() # integral E-P
    intHeatFlxa = tmp.copy() # integral heat flux Atl
    intWatFlxa  = tmp.copy() # integral E-P Atl
    intHeatFlxp = tmp.copy() # integral heat flux Pac
    intWatFlxp  = tmp.copy() # integral E-P Pac
    intHeatFlxi = tmp.copy() # integral heat flux Ind
    intWatFlxi  = tmp.copy() # integral E-P Ind
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
        tost = regridObj(tos [t,:,:])
        sost = regridObj(sos [t,:,:])
        heft = regridObj(qnet[t,:,:])
        empt = regridObj(emp [t,:,:])
        tost.mask = maski
        sost.mask = maski
        heft.mask = maski
        empt.mask = maski
        tost = maskVal(tost, valmask)
        sost = maskVal(sost, valmask)
        heft = maskVal(heft, valmask)
        empt = maskVal(empt, valmask)
        # define basin heat and water fluxes
        hefta = heft*1.
        heftp = heft*1.
        hefti = heft*1.
        hefta.mask = maskAtl
        heftp.mask = maskPac
        hefti.mask = maskInd
        hefta = maskVal(hefta, valmask)
        heftp = maskVal(heftp, valmask)
        hefti = maskVal(hefti, valmask)
        #
        empta = empt*1.
        emptp = empt*1.
        empti = empt*1.
        empta.mask = maskAtl
        emptp.mask = maskPac
        empti.mask = maskInd
        empta = maskVal(empta, valmask)
        emptp = maskVal(emptp, valmask)
        empti = maskVal(empti, valmask)
        #
        # Compute density
        rhon[t,...] = eosNeutral(tost.data, sost.data) - 1000.
        rhonl = rhon.data[t,...]
         # Compute buoyancy/density flux as mass fluxes in kg/m2/s (SI unts)
        #  convwf : kg/m2/s = mm/s -> m/s
        convwf = 1.e-3
        pres = tost.data*0.
        denflxh[t,...] = (-sw.alpha(sost.data,tost.data,pres)/sw.cp(sost.data,tost.data,pres))*heft.data
        denflxw[t,...] = (rhonl+1000.)*sw.beta(sost.data,tost.data,pres)*sost.data*empt.data*convwf
        denflx [t,...] = denflxh[t,...] + denflxw[t,...]
        denflx [t,...].mask  = maski
        denflxh[t,...].mask  = maski
        denflxw[t,...].mask  = maski
        denflx [t,...] = maskVal(denflx [t,...], valmask)
        denflxh[t,...] = maskVal(denflxh[t,...], valmask)
        denflxw[t,...] = maskVal(denflxw[t,...], valmask)
        dflxh = denflxh.data[t,:,:]
        dflxw = denflxw.data[t,:,:]
        #
        # Transformation (integral of density flux on density outcrops)
        for ks in range(N_s-1):
            # find indices of points in density bin
            # Global
            idxbin = npy.argwhere( (rhonl >= sigrid[ks]) & (rhonl < sigrid[ks+1]) ).transpose()
            idj = idxbin[0] ; idi = idxbin[1]
            transfh[t,ks] = cdu.averager(dflxh[idj,idi] * areai[idj,idi], axis=0, action='sum')/del_s[ks]
            transfw[t,ks] = cdu.averager(dflxw[idj,idi] * areai[idj,idi], axis=0, action='sum')/del_s[ks]
            areabin[t,ks] = cdu.averager(areai[idxbin[0],idxbin[1]], axis=0, action='sum')
            # Basin
            idxbina = npy.argwhere( (rhonl*maskAtl >= sigrid[ks]) & (rhonl*maskAtl < sigrid[ks+1]) ).transpose()
            idj = idxbina[0] ; idi = idxbina[1]
            transfha[t,ks] = cdu.averager(dflxh[idj,idi] * areai[idj,idi], axis=0, action='sum')/del_s[ks]
            transfwa[t,ks] = cdu.averager(dflxw[idj,idi] * areai[idj,idi], axis=0, action='sum')/del_s[ks]
            areabina[t,ks] = cdu.averager(areai[idj,idi], axis=0, action='sum')
            #
            idxbinp = npy.argwhere( (rhonl*maskPac >= sigrid[ks]) & (rhonl*maskPac < sigrid[ks+1]) ).transpose()
            idj = idxbinp[0] ; idi = idxbinp[1]
            transfhp[t,ks] = cdu.averager(dflxh[idj,idi] * areai[idj,idi], axis=0, action='sum')/del_s[ks]
            transfwp[t,ks] = cdu.averager(dflxw[idj,idi] * areai[idj,idi], axis=0, action='sum')/del_s[ks]
            areabinp[t,ks] = cdu.averager(areai[idj,idi], axis=0, action='sum')
            #
            idxbini = npy.argwhere( (rhonl*maskInd >= sigrid[ks]) & (rhonl*maskInd < sigrid[ks+1]) ).transpose()
            idj = idxbini[0] ; idi = idxbini[1]
            transfhi[t,ks] = cdu.averager(dflxh[idj,idi] * areai[idj,idi], axis=0, action='sum')/del_s[ks]
            transfwi[t,ks] = cdu.averager(dflxw[idj,idi] * areai[idj,idi], axis=0, action='sum')/del_s[ks]
            areabini[t,ks] = cdu.averager(areai[idj,idi], axis=0, action='sum')
            
        # last bin
        # Global
        idxbin = npy.argwhere( (rhonl >= sigrid[N_s-1])).transpose()
        idj = idxbin[0] ; idi = idxbin[1]
        transfh[t,N_s] = cdu.averager(dflxh[idj,idi] * areai[idj,idi], axis=0, action='sum')/del_s[ks]
        transfw[t,N_s] = cdu.averager(dflxw[idj,idi] * areai[idj,idi], axis=0, action='sum')/del_s[ks]
        areabin[t,N_s] = cdu.averager(areai[idj,idi], axis=0, action='sum')
        # Basins
        idxbina = npy.argwhere( (rhonl*maskAtl >= sigrid[N_s-1])).transpose()
        idj = idxbina[0] ; idi = idxbina[1]
        transfha[t,N_s] = cdu.averager(dflxh[idj,idi] * areai[idj,idi], axis=0, action='sum')/del_s[ks]
        transfwa[t,N_s] = cdu.averager(dflxw[idj,idi] * areai[idj,idi], axis=0, action='sum')/del_s[ks]
        areabina[t,N_s] = cdu.averager(areai[idj,idi], axis=0, action='sum')
        # 
        idxbinp = npy.argwhere( (rhonl*maskPac >= sigrid[N_s-1])).transpose()
        idj = idxbinp[0] ; idi = idxbinp[1]
        transfhp[t,N_s] = cdu.averager(dflxh[idj,idi] * areai[idj,idi], axis=0, action='sum')/del_s[ks]
        transfwp[t,N_s] = cdu.averager(dflxw[idj,idi] * areai[idj,idi], axis=0, action='sum')/del_s[ks]
        areabinp[t,N_s] = cdu.averager(areai[idj,idi], axis=0, action='sum')
        # 
        idxbini = npy.argwhere( (rhonl*maskInd >= sigrid[N_s-1])).transpose()
        idj = idxbini[0] ; idi = idxbini[1]
        transfhi[t,N_s] = cdu.averager(dflxh[idj,idi] * areai[idj,idi], axis=0, action='sum')/del_s[ks]
        transfwi[t,N_s] = cdu.averager(dflxw[idj,idi] * areai[idj,idi], axis=0, action='sum')/del_s[ks]
        areabini[t,N_s] = cdu.averager(areai[idj,idi], axis=0, action='sum')
        # Total transformation
        transf[t,:] = transfh[t,:] + transfw[t,:]        
        transfa[t,:] = transfha[t,:] + transfwa[t,:]        
        transfp[t,:] = transfhp[t,:] + transfwp[t,:]        
        transfi[t,:] = transfhi[t,:] + transfwi[t,:] 
        # Formation = divergence of transformation in density space
        # done in postpro:  form [t,ks] = -(transf [t,ks+1] - transf [t,ks])
        #
        # domain integrals
        # heat flux (conv W -> PW)
        convt  = 1.e-15
        intHeatFlx [t] = cdu.averager(npy.reshape(heft*areai , (Nji*Nii)), action='sum')*dt*convt
        intHeatFlxa[t] = cdu.averager(npy.reshape(hefta*areai, (Nji*Nii)), action='sum')*dt*convt
        intHeatFlxp[t] = cdu.averager(npy.reshape(heftp*areai, (Nji*Nii)), action='sum')*dt*convt
        intHeatFlxi[t] = cdu.averager(npy.reshape(hefti*areai, (Nji*Nii)), action='sum')*dt*convt
        # fw flux (conv mm -> m and m3/s to Sv)
        convw = 1.e-3*1.e-6
        intWatFlx [t]  = cdu.averager(npy.reshape(empt*areai , (Nji*Nii)), action='sum')*dt*convw
        intWatFlxa[t]  = cdu.averager(npy.reshape(empta*areai, (Nji*Nii)), action='sum')*dt*convw
        intWatFlxp[t]  = cdu.averager(npy.reshape(emptp*areai, (Nji*Nii)), action='sum')*dt*convw
        intWatFlxi[t]  = cdu.averager(npy.reshape(empti*areai, (Nji*Nii)), action='sum')*dt*convw

#        if debugp:
#            print '    integral Q flux ',t,intHeatFlx [t], intHeatFlxa[t], intHeatFlxp[t], intHeatFlxi[t]
#            print '    integral W flux ',t,intWatFlx [t], intWatFlxa[t], intWatFlxp[t], intWatFlxi[t]
      
    # Wash mask over variables
    maskt        = mv.masked_values(rhon, valmask).mask
    denflx.mask  = maskt
    denflxh.mask = maskt
    denflxw.mask = maskt
    denflx       = maskVal(denflx , valmask)
    denflxh      = maskVal(denflxh, valmask)
    denflxw      = maskVal(denflxw, valmask)

    maskin       = mv.masked_values(transf, valmask).mask
    transfh._FillValue = valmask
    transfw._FillValue = valmask
    transf._FillValue  = valmask
    transfh      = maskVal(transfh, valmask)
    transfw      = maskVal(transfw, valmask)
    transf       = maskVal(transf , valmask)
    #
    maskin       = mv.masked_values(transfa, valmask).mask
    transfha._FillValue = valmask
    transfwa._FillValue = valmask
    transfa._FillValue  = valmask
    transfha      = maskVal(transfha, valmask)
    transfwa      = maskVal(transfwa, valmask)
    transfa       = maskVal(transfa , valmask)
    #
    maskin       = mv.masked_values(transfp, valmask).mask
    transfhp._FillValue = valmask
    transfwp._FillValue = valmask
    transfp._FillValue  = valmask
    transfhp      = maskVal(transfhp, valmask)
    transfwp      = maskVal(transfwp, valmask)
    transfp       = maskVal(transfp , valmask)
    #
    maskin       = mv.masked_values(transfi, valmask).mask
    transfhi._FillValue = valmask
    transfwi._FillValue = valmask
    transfi._FillValue  = valmask
    transfhi      = maskVal(transfhi, valmask)
    transfwi      = maskVal(transfwi, valmask)
    transfi       = maskVal(transfi , valmask)
       
    #+ create a basins variables (loop on n masks)

    #
    # Output files as netCDF
    # Density flux (3D: time, lon, lat)
    convw = 1.e6
    rhon    = cdm.createVariable(rhon          , axes = [time, lati, loni], id = 'densurf')
    denFlx  = cdm.createVariable(denflx*convw  , axes = [time, lati, loni], id = 'denflux')
    denFlxh = cdm.createVariable(denflxh*convw , axes = [time, lati, loni], id = 'hdenflx')
    denFlxw = cdm.createVariable(denflxw*convw , axes = [time, lati, loni], id = 'wdenflx')
    denFlx.long_name   = 'Surface density'
    denFlx.units       = 'kg.m-3 (anomaly, minus 1000)'
    denFlx.long_name   = 'Total density flux'
    denFlx.units       = '1.e-6 kg/m2/s'
    denFlxh.long_name  = 'Heat density flux'
    denFlxh.units      = '1.e-6 kg/m2/s'
    denFlxw.long_name  = 'Water density flux'
    denFlxw.units      = '1.e-6 kg/m2/s'
    #
    # Transformation (2D: time, sigma)
    convw = 1.e-6
    totTransf   = cdm.createVariable(transf*convw  , axes = [time, s_axis], id = 'trsftot')
    hefTransf   = cdm.createVariable(transfh*convw , axes = [time, s_axis], id = 'trsfhef')
    wfoTransf   = cdm.createVariable(transfw*convw , axes = [time, s_axis], id = 'trsfwfo')
    totTransfa  = cdm.createVariable(transfa*convw , axes = [time, s_axis], id = 'trsftotAtl')
    hefTransfa  = cdm.createVariable(transfha*convw, axes = [time, s_axis], id = 'trsfhefAtl')
    wfoTransfa  = cdm.createVariable(transfwa*convw, axes = [time, s_axis], id = 'trsfwfoAtl')
    totTransfp  = cdm.createVariable(transfp*convw , axes = [time, s_axis], id = 'trsftotPac')
    hefTransfp  = cdm.createVariable(transfhp*convw, axes = [time, s_axis], id = 'trsfhefPac')
    wfoTransfp  = cdm.createVariable(transfwp*convw, axes = [time, s_axis], id = 'trsfwfoPac')
    totTransfi  = cdm.createVariable(transfi*convw , axes = [time, s_axis], id = 'trsftotInd')
    hefTransfi  = cdm.createVariable(transfhi*convw, axes = [time, s_axis], id = 'trsfhefInd')
    wfoTransfi  = cdm.createVariable(transfwi*convw, axes = [time, s_axis], id = 'trsfwfoInd')
    totTransf.long_name   = 'Total transformation'
    totTransf.units       = 'Sv'
    hefTransf.long_name   = 'Heat flux transformation'
    hefTransf.units       = 'Sv'
    wfoTransf.long_name   = 'Water flux transformation'
    wfoTransf.units       = 'Sv'
    totTransfa.long_name  = 'Atl. Total transformation'
    totTransfa.units      = 'Sv'
    hefTransfa.long_name  = 'Atl. Heat flux transformation'
    hefTransfa.units      = 'Sv'
    wfoTransfa.long_name  = 'Atl. Water flux transformation'
    wfoTransfa.units      = 'Sv'
    totTransfp.long_name  = 'Pac. Total transformation'
    totTransfp.units      = 'Sv'
    hefTransfp.long_name  = 'Pac. Heat flux transformation'
    hefTransfp.units      = 'Sv'
    wfoTransfp.long_name  = 'Pac. Water flux transformation'
    wfoTransfp.units      = 'Sv'
    totTransfi.long_name  = 'Ind. Total transformation'
    totTransfi.units      = 'Sv'
    hefTransfi.long_name  = 'Ind. Heat flux transformation'
    hefTransfi.units      = 'Sv'
    wfoTransfi.long_name  = 'Ind. Water flux transformation'
    wfoTransfi.units      = 'Sv'
    #
    # Integral heat and emp fux (1D: time)
    intQFlx   = cdm.createVariable(intHeatFlx  , axes = [time], id = 'intQflx')
    intWFlx   = cdm.createVariable(intWatFlx   , axes = [time], id = 'intWflx')
    intQFlxa  = cdm.createVariable(intHeatFlxa , axes = [time], id = 'intQflxAtl')
    intWFlxa  = cdm.createVariable(intWatFlxa  , axes = [time], id = 'intWflxAtl')
    intQFlxp  = cdm.createVariable(intHeatFlxp , axes = [time], id = 'intQflxPac')
    intWFlxp  = cdm.createVariable(intWatFlxp  , axes = [time], id = 'intWflxPac')
    intQFlxi  = cdm.createVariable(intHeatFlxi , axes = [time], id = 'intQflxInd')
    intWFlxi  = cdm.createVariable(intWatFlxi  , axes = [time], id = 'intWflxInd')
    intQFlx.long_name   = 'Integral Surface Heat Flux'
    intQFlx.units       = 'PW'
    intWFlx.long_name   = 'Integral Surface E minus P'
    intWFlx.units       = 'Sv'
    intQFlxa.long_name  = 'Atl. Integral Surface Heat Flux'
    intQFlxa.units      = 'PW'
    intWFlxa.long_name  = 'Atl. Integral Surface E minus P'
    intWFlxa.units      = 'Sv'
    intQFlxp.long_name  = 'Pac. Integral Surface Heat Flux'
    intQFlxp.units      = 'PW'
    intWFlxp.long_name  = 'Pac. Integral Surface E minus P'
    intWFlxp.units      = 'Sv'
    intQFlxi.long_name  = 'Ind. Integral Surface Heat Flux'
    intQFlxi.units      = 'PW'
    intWFlxi.long_name  = 'Ind. Integral Surface E minus P'
    intWFlxi.units      = 'Sv'

    outFile_f.write(rhon)
    outFile_f.write(denFlx)
    outFile_f.write(denFlxh)
    outFile_f.write(denFlxw)
    outFile_f.write(totTransf)
    outFile_f.write(hefTransf)
    outFile_f.write(wfoTransf)
    outFile_f.write(totTransfa)
    outFile_f.write(hefTransfa)
    outFile_f.write(wfoTransfa)
    outFile_f.write(totTransfp)
    outFile_f.write(hefTransfp)
    outFile_f.write(wfoTransfp)
    outFile_f.write(totTransfi)
    outFile_f.write(hefTransfi)
    outFile_f.write(wfoTransfi)
    outFile_f.write(intQFlx)
    outFile_f.write(intWFlx)
    outFile_f.write(intQFlxa)
    outFile_f.write(intWFlxa)
    outFile_f.write(intQFlxp)
    outFile_f.write(intWFlxp)
    outFile_f.write(intQFlxi)
    outFile_f.write(intWFlxi)
    #
    # File global attributes
    for i in range(0,len(file_dic)):
        dm = file_dic[i]
        setattr(outFile_f,dm[0],dm[1])
        post_txt = 'Density flux via surfTransf using delta_sigma = '+str(del_s1)+' and '+str(del_s2)
        setattr(outFile_f, 'Post_processing_history', post_txt)

    ftos.close()
    fsos.close()
    fhef.close()
    fwfo.close()
    outFile_f.close()
    print ' Wrote file: ',outFile




    # CPU use
    print
    print ' CPU use', timc.clock() - cpu0

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

EG   8 Oct 2014     - Started file
PJD 22 Nov 2014     - Updated to comment out unused statements and imports

- TODO:
    - add ekman pumping bining (cf wcurl in densit)
    - add time chunks (as in binDensity.py) to deal with CPU/memory balance
        - read variables in time loop
        - add extension in outputs files


@author: eguil
"""

import cdms2 as cdm
import MV2 as mv
import os, resource
import numpy as npy
import cdutil as cdu

from binDensity import maskVal
from binDensity import eosNeutral
from binDensity import rhonGrid
from binDensity import computeAreaScale
import time as timc
import ESMP
from cdms2 import CdmsRegrid
from durolib import fixVarUnits
import seawater as sw

#
# inits
# -----
#

def surfTransf(fileFx, fileTos, fileSos, fileHef, fileWfo, varNames, outFile, debug=True, timeint='all',noInterp=False, domain='global'):
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
    - varNames[4]               - 1D array containing the names of the variables
    - outFile(str)              - output file with full path specified.
    - debug <optional>          - boolean value
    - timeint <optional>        - specify temporal step for binning <init_idx>,<ncount>
    - noInterp <optional>       - if true no interpolation to target grid
    - domain <optional>         - specify domain for averaging when interpolated to WOA grid ('global','north',
                                  'north40', 'south' for now)

    Outputs:
    --------
    - netCDF file with monthly surface rhon, density fluxes, transformation (global and per basin)
    - use cdo yearmean to compute annual mean

    Usage:
    ------
    '>>> from binDensity import surfTransf
    '>>> surfTransf(file_fx, file_tos, file_sos, file_hef, file_wfo, [var1,var2,var3,var4]./output.nc, debug=True,timeint='all')

    Notes:
    -----
    - EG   8 Oct 2014   - Initial function write and tests ok
    - PJD 22 Nov 2014   - Code cleanup
    - EG   4 Oct 2017   - code on ciclad, more cleanup and options
    - EG  12 Sep 2018   - Add North vs. South calculation

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
    #timeax = ftos.getAxis('time')
    timeax = ftos.getAxis('time_counter')
    #print 'timeax'
    #print timeax
    #
    # Dates to read
    if timeint == 'all':
        tmin = 0
        tmax = timeax.shape[0]
        timeaxis = timeax
    else:
        tmin = int(timeint.split(',')[0]) - 1
        tmax = tmin + int(timeint.split(',')[1])
        # update time axis
        timeaxis   = cdm.createAxis(timeax[tmin:tmax])
        timeaxis.id       = 'time'
        timeaxis.units    = timeax.units
        timeaxis.designateTime()
        #print timeaxis

    if debugp:
        print; print ' Debug mode'
 
    # Read file attributes to carry on to output files
    list_file   = ftos.attributes.keys()
    file_dic    = {}
    for i in range(0,len(list_file)):
        file_dic[i]=list_file[i],ftos.attributes[list_file[i] ]
    #
    # Read data
        
    # varnames
    tos_name = varNames[0]
    sos_name = varNames[1]
    hef_name = varNames[2]
    wfo_name = varNames[3]

    if debugp:
        print' Read ',tos_name, sos_name,tmin,tmax
    tos = ftos(tos_name , time = slice(tmin,tmax))
    sos = fsos(sos_name , time = slice(tmin,tmax))
    if debugp:
        print' Read ',hef_name, wfo_name
    qnet = fhef(hef_name, time = slice(tmin,tmax))
    try:
        emp  = fwfo(wfo_name , time = slice(tmin,tmax))
        empsw = 0
    except Exception,err:
        emp  = fwfo('wfos' , time = slice(tmin,tmax))
        print ' Reading concentration dillution fresh water flux'
        empsw = 0
    tos_h = ftos[tos_name]
    if debugp:
        print tos_h
    #
    # Read input grid
    ingrid = tos.getGrid()
    #
    # Read cell area
    #ff = cdm.open(fileFx)
    #area = ff('areacello')
    #ff.close()
    #areain = area.data
    #
    # Define dimensions
    N_i = int(tos.shape[1])
    N_j = int(tos.shape[0])
    #
    # Define sigma grid 
    rho_min = 22
    rho_int = 25
    rho_max = 29
    del_s1  = 0.1
    del_s2  = 0.05
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
    # Write 3D density flux ?
    #
    writedenflx = False
    #
    # Monthly transformation
    if os.path.exists(outFile):
        os.remove(outFile)
    outFile_f = cdm.open(outFile,'w')
    # Define dimensions
    N_i = int(tos.shape[2])
    N_j = int(tos.shape[1])
    N_t = int(tos.shape[0])
    print ' ==> dimensions N_t, N_j, N_i:', N_t, N_j, N_i
    try:
        valmask = tos.missing_value
        if valmask == None:
            print 'EC-EARTH missing_value fix'
            valmask = 1.e20            
            sos = mv.masked_equal(sos,0.)
            print sos.count()
            sos.data[:] = sos.filled(valmask)
            tos.mask = sos.mask
            tos.data[:] = tos.filled(valmask)
            qnet.mask = sos.mask
            qnet.data[:] = qnet.filled(valmask)
            emp.mask = sos.mask
            emp.data[:] = emp.filled(valmask)
    except Exception,err:
        print 'Exception: ',err
        if 'EC-EARTH' == modeln:
            print 'EC-EARTH missing_value fix'
            valmask = 1.e20
            sos = mv.masked_equal(sos,0.)
            print sos.count()
            sos.data[:] = sos.filled(valmask)
            tos.mask = sos.mask
            tos.data[:] = tos.filled(valmask)
            qnet.mask = sos.mask
            qnet.data[:] = qnet.filled(valmask)
            emp.mask = sos.mask
            emp.data[:] = emp.filled(valmask)
    # Read masking value
    # added if wfcorr == masked values everywhere
    emp.mask = sos.mask
    emp.data[:] = emp.filled(valmask)

    # Test variable units
    [sos,sosFixed] = fixVarUnits(sos,'sos',True)#,'logfile.txt')
    if sosFixed:
        print '     sos: units corrected'
    [tos,tosFixed] = fixVarUnits(tos,'thetao',True)#,'logfile.txt')
    if tosFixed:
        print '     tos: units corrected'        
    # Physical inits
    #P = 0          # surface pressure
    # find non-masked points
    #maskin = mv.masked_values(tos.data[0], valmask).mask 
    #nomask = npy.equal(maskin,0)
    #
    # target horizonal grid for interp 
    if noInterp:
        outgrid = ingrid
        fileg = '/data/vestella/Masks/Mask_Convect_Atlantic_40N.nc'
        gt = cdm.open(fileg)
        maskr = gt('mask_part')
        maski = maskr[...]
        gt.close
        #maski = tos[0,...].mask
        fileg = '/data/vestella/Masks/Mask_Convect_GS.nc'
        gt = cdm.open(fileg)
        maskr = gt('mask_part')
        maskAtl = maskr[...]
        gt.close
        fileg = '/data/vestella/Masks/Mask_Convect_SI.nc'
        gt = cdm.open(fileg)
        maskr = gt('mask_part')
        maskPac = maskr[...]
        gt.close
        fileg = '/data/vestella/Masks/Mask_Convect_LS.nc'
        gt = cdm.open(fileg)
        maskr = gt('mask_part')
        maskInd = maskr[...]
        gt.close

        # Read area of target grid and zonal sums
        areai = npy.ma.ones([N_j, N_i], dtype='float32')*0.
        fileg = '/data/igcmg/database/grids/ORCA2.3_area.nc'
        gt = cdm.open(fileg)
        arear = gt('area')
        arear_h = gt['area']
        areai = arear.data[:,:]

        loni  = tos_h.getLongitude()
        lati  = tos_h.getLatitude()

        Nii   = int(loni.shape[1])
        Nji   = int(lati.shape[0])

        gt.close
    else:
        # Interpolate on WOA grid
        #fileg = '/work/guilyardi/Density_bining/WOD13_masks.nc'
        #fileg = '/export/durack1/git/Density_bining/140807_WOD13_masks.nc'
        fileg = '170224_WOD13_masks.nc'
        gt = cdm.open(fileg)
        maskg = gt('basinmask3')
        outgrid = maskg.getGrid()
        # global mask
        maski = maskg.mask
        # regional masks
        maskAtl = maski*1 ; maskAtl[...] = False
        idxa = npy.argwhere(maskg == 1).transpose()
        maskAtl[idxa[0],idxa[1]] = True
        maskPac = maski*1 ; maskPac[...] = False
        idxp = npy.argwhere(maskg == 2).transpose()
        maskPac[idxp[0],idxp[1]] = True
        maskInd = maski*1 ; maskInd[...] = False
        idxi = npy.argwhere(maskg == 3).transpose()
        maskInd[idxi[0],idxi[1]] = True
        #masks = [maski, maskAtl, maskPac, maskInd]
        loni    = maskg.getLongitude()
        lati    = maskg.getLatitude()
        Nii     = int(loni.shape[0])
        Nji     = int(lati.shape[0])
        # Compute area of target grid and zonal sums
        areai, scalex, scaley = computeAreaScale(loni[:], lati[:])
        #areai = gt('basinmask3_area')
        #print areai.shape
        #print areai[:,90]
        gt.close()

        # Reduce domain to North/South ?
        print ' Domain : ',domain
        if domain == 'north':
            lati2d = npy.tile(lati, Nii).reshape(Nii, Nji).transpose()
            indn = npy.argwhere(lati2d <= 0).transpose()
            maski  [indn[0],indn[1]] = False
            maskAtl[indn[0],indn[1]] = False
            maskPac[indn[0],indn[1]] = False
            maskInd[indn[0],indn[1]] = False
        elif domain == 'north40':
            lati2d = npy.tile(lati, Nii).reshape(Nii,Nji).transpose()
            indn = npy.argwhere(lati2d <= 40).transpose()
            maski  [indn[0],indn[1]] = False
            maskAtl[indn[0],indn[1]] = False
            maskPac[indn[0],indn[1]] = False
            maskInd[indn[0],indn[1]] = False
        elif domain == 'south':
            lati2d = npy.tile(lati, Nii).reshape(Nii,Nji).transpose()
            indn = npy.argwhere(lati2d >= 0).transpose()
            maski  [indn[0],indn[1]] = False
            maskAtl[indn[0],indn[1]] = False
            maskPac[indn[0],indn[1]] = False
            maskInd[indn[0],indn[1]] = False

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
    # Interpolation init (regrid)
    if noInterp == False:
        ESMP.ESMP_Initialize()
        regridObj = CdmsRegrid(ingrid, outgrid, denflxh.dtype, missing = valmask, regridMethod = 'linear', regridTool = 'esmf')
    # init integration intervals
    dt   = 1./float(N_t) 


    # Bin on density grid
    for t in range(N_t):
        if noInterp:
            tost = tos [t,:,:]
            sost = sos [t,:,:]
            heft = qnet[t,:,:]
            empt = emp [t,:,:]
        else:
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
        #
        empta = empt*1.
        emptp = empt*1.
        empti = empt*1.
        #
        # Compute density
        rhon[t,...] = eosNeutral(tost.data, sost.data) - 1000.
        rhon[t,...].mask  = maski
        rhon[t,...] = maskVal(rhon[t,...], valmask)
        rhonl = rhon.data[t,...]
        # Compute buoyancy/density flux as mass fluxes in kg/m2/s (SI units)
        #  convwf : kg/m2/s = mm/s -> m/s
        convwf = 1.e-3
        pres = tost.data*0.
        denflxh[t,...] = (-sw.alpha(sost.data,tost.data,pres)/sw.cp(sost.data,tost.data,pres))*heft.data
        if empsw == 0:
            denflxw[t,...] = (rhonl+1000.)*sw.beta(sost.data,tost.data,pres)*sost.data*empt.data*convwf
        else:
            denflxw[t,...] = (rhonl+1000.)*sw.beta(sost.data,tost.data,pres)*empt.data*convwf            
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
            areabin[t,ks] = cdu.averager(areai[idj,idi], axis=0, action='sum')

            # Basin
            #print maskAtl.shape
            #print maskAtl
            #print rhonl.shape
            #print rhonl

            idxbina = npy.argwhere( (rhonl*maskAtl >= sigrid[ks]) & (rhonl*maskAtl < sigrid[ks+1]) ).transpose()
            idj = idxbina[0] ; idi = idxbina[1]
            #print t,ks,dflxh.shape,areai.shape,idj,idi
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
        intHeatFlxa[t] = cdu.averager(npy.reshape(hefta*areai*maskAtl, (Nji*Nii)), action='sum')*dt*convt
        intHeatFlxp[t] = cdu.averager(npy.reshape(heftp*areai*maskPac, (Nji*Nii)), action='sum')*dt*convt
        intHeatFlxi[t] = cdu.averager(npy.reshape(hefti*areai*maskInd, (Nji*Nii)), action='sum')*dt*convt
        # fw flux (conv mm -> m and m3/s to Sv)
        convw = 1.e-3*1.e-6
        intWatFlx [t]  = cdu.averager(npy.reshape(empt*areai , (Nji*Nii)), action='sum')*dt*convw
        intWatFlxa[t]  = cdu.averager(npy.reshape(empta*areai*maskAtl, (Nji*Nii)), action='sum')*dt*convw
        intWatFlxp[t]  = cdu.averager(npy.reshape(emptp*areai*maskPac, (Nji*Nii)), action='sum')*dt*convw
        intWatFlxi[t]  = cdu.averager(npy.reshape(empti*areai*maskInd, (Nji*Nii)), action='sum')*dt*convw

        if debugp and t == 0:
            print '    integral Q flux ',t,intHeatFlx [t], intHeatFlxa[t], intHeatFlxp[t], intHeatFlxi[t]
            print '    integral W flux ',t,intWatFlx [t], intWatFlxa[t], intWatFlxp[t], intWatFlxi[t]
      
    # Wash mask over variables
    maskt        = mv.masked_values(rhon, valmask).mask
    denflx.mask  = maskt
    denflxh.mask = maskt
    denflxw.mask = maskt
    denflx       = maskVal(denflx , valmask)
    denflxh      = maskVal(denflxh, valmask)
    denflxw      = maskVal(denflxw, valmask)

    #maskin       = mv.masked_values(transf, valmask).mask
    transfh._FillValue = valmask
    transfw._FillValue = valmask
    transf._FillValue  = valmask
    transfh      = maskVal(transfh, valmask)
    transfw      = maskVal(transfw, valmask)
    transf       = maskVal(transf , valmask)
    #
    #maskin       = mv.masked_values(transfa, valmask).mask
    transfha._FillValue = valmask
    transfwa._FillValue = valmask
    transfa._FillValue  = valmask
    transfha      = maskVal(transfha, valmask)
    transfwa      = maskVal(transfwa, valmask)
    transfa       = maskVal(transfa , valmask)
    #
    #maskin       = mv.masked_values(transfp, valmask).mask
    transfhp._FillValue = valmask
    transfwp._FillValue = valmask
    transfp._FillValue  = valmask
    transfhp      = maskVal(transfhp, valmask)
    transfwp      = maskVal(transfwp, valmask)
    transfp       = maskVal(transfp , valmask)
    #
    #maskin       = mv.masked_values(transfi, valmask).mask
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
    if writedenflx:
        rhon    = cdm.createVariable(rhon          , axes = [timeaxis, lati, loni], id = 'densurf')
        denFlx  = cdm.createVariable(denflx*convw  , axes = [timeaxis, lati, loni], id = 'denflux')
        denFlxh = cdm.createVariable(denflxh*convw , axes = [timeaxis, lati, loni], id = 'hdenflx')
        denFlxw = cdm.createVariable(denflxw*convw , axes = [timeaxis, lati, loni], id = 'wdenflx')
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
    totTransf   = cdm.createVariable(transf*convw  , axes = [timeaxis, s_axis], id = 'trsftot')
    hefTransf   = cdm.createVariable(transfh*convw , axes = [timeaxis, s_axis], id = 'trsfhef')
    wfoTransf   = cdm.createVariable(transfw*convw , axes = [timeaxis, s_axis], id = 'trsfwfo')
    totTransfa  = cdm.createVariable(transfa*convw , axes = [timeaxis, s_axis], id = 'trsftotAtl')
    hefTransfa  = cdm.createVariable(transfha*convw, axes = [timeaxis, s_axis], id = 'trsfhefAtl')
    wfoTransfa  = cdm.createVariable(transfwa*convw, axes = [timeaxis, s_axis], id = 'trsfwfoAtl')
    totTransfp  = cdm.createVariable(transfp*convw , axes = [timeaxis, s_axis], id = 'trsftotPac')
    hefTransfp  = cdm.createVariable(transfhp*convw, axes = [timeaxis, s_axis], id = 'trsfhefPac')
    wfoTransfp  = cdm.createVariable(transfwp*convw, axes = [timeaxis, s_axis], id = 'trsfwfoPac')
    totTransfi  = cdm.createVariable(transfi*convw , axes = [timeaxis, s_axis], id = 'trsftotInd')
    hefTransfi  = cdm.createVariable(transfhi*convw, axes = [timeaxis, s_axis], id = 'trsfhefInd')
    wfoTransfi  = cdm.createVariable(transfwi*convw, axes = [timeaxis, s_axis], id = 'trsfwfoInd')
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
    intQFlx   = cdm.createVariable(intHeatFlx  , axes = [timeaxis], id = 'intQflx')
    intWFlx   = cdm.createVariable(intWatFlx   , axes = [timeaxis], id = 'intWflx')
    intQFlxa  = cdm.createVariable(intHeatFlxa , axes = [timeaxis], id = 'intQflxAtl')
    intWFlxa  = cdm.createVariable(intWatFlxa  , axes = [timeaxis], id = 'intWflxAtl')
    intQFlxp  = cdm.createVariable(intHeatFlxp , axes = [timeaxis], id = 'intQflxPac')
    intWFlxp  = cdm.createVariable(intWatFlxp  , axes = [timeaxis], id = 'intWflxPac')
    intQFlxi  = cdm.createVariable(intHeatFlxi , axes = [timeaxis], id = 'intQflxInd')
    intWFlxi  = cdm.createVariable(intWatFlxi  , axes = [timeaxis], id = 'intWflxInd')
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

    if writedenflx:
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
    print ' [ Time stamp',(timc.strftime("%d/%m/%Y %H:%M:%S")),']'
    print ' CPU use', timc.clock() - cpu0
    print ' Max memory use',resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1.e6,'GB'


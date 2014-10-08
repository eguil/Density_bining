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
---------------------------------------------------------------------------------

EG 8 Oct 2014     - Started file
                    - TODO: 
                     - add ekman pumping bining (cf wcurl in densit)
                     - modify alpha and betar functions to more physical
test

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
#
# inits
# -----
#
def alpha (t, s):
    # compute alpha=-1/rho (d rho / d T)
    dt = 0.05
    siga = eosNeutral(t, s)
    sigb = eosNeutral(t+0.05, s)
    alpha = -0.001*(sigb-siga)/dt/(1.+1.e-3*siga)
    return alpha
def betar (t, s):
    # compute beta= 1/rho (d rho / d S)
    ds = 0.01
    siga = eosNeutral(t, s)-1000.
    sigb = eosNeutral(t, s+ds)-1000.    
    beta = 0.001*(sigb-siga)/ds/(1.+1.e-3*siga)
    return beta
def cpsw (t, s, p):
    # Specific heat of sea water (J/KG C)
    CP1 = 0.
    CP2 = 0.
    S = s
    T = t
    P = p
    SR = npy.ma.sqrt(S)
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

def surfTransf(fileFx, fileTos, fileSos, fileHef, fileWfo, outFile, debug=True,timeint='all'):
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

    Usage:
    ------
    >>> from binDensity import surfTransf
    >>> surfTransf(file_fx, file_tos, file_sos, file_hef, file_wfo, ./output.nc, debug=True,timeint='all')

    Notes:
    -----
    - EG 8 Oct 2014   - Initial function write
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
    #
    # Read masking value
    valmask = tos._FillValue
    #
    # Read time and grid
    lon  = tos.getLongitude()
    lat  = tos.getLatitude()
    ingrid = tos.getGrid()
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
    s_axis.units = ''
    s_axis.designateLevel()
    #
    # Monthly transformation
    outFile = replace(outFile,'.mo.','.an.')
    if os.path.exists(outFile):
        os.remove(outFile)
    outFile_f = cdm.open(outFile,'w')
    # Define dimensions
    N_i = int(tos.shape[2])
    N_j = int(tos.shape[1])
    N_t = int(tos.shape[0])
    # Read masking value
    valmask = tos._FillValue
    # reorganise i,j dims in single dimension data
    tos  = npy.reshape(tos, (N_t, N_i*N_j))
    sos  = npy.reshape(sos, (N_t, N_i*N_j))
    emp  = npy.reshape(emp, (N_t, N_i*N_j))
    qnet = npy.reshape(qnet, (N_t, N_i*N_j))
    areain = npy.reshape(areain, (N_i*N_j))
    print 'tos', tos.data[0,1000:1100]
    print 'sos', sos.data[0,1000:1100]
    # Test variable units
    [sos,sosFixed] = fixVarUnits(sos,'sos',True)#,'logfile.txt')
    if sosFixed:
        print '     sos: units corrected'
    [tos,tosFixed] = fixVarUnits(tos,'thetao',True)#,'logfile.txt')
    if tosFixed:
        print '     tos: units corrected'        
    # Physical inits
    P = 0          # surface pressure
    conwf = 1.e-3  # kg/m2/s=mm/s -> m/s
    # find non-masked points
    maskin = mv.masked_values(tos.data[0], valmask).mask 
    nomask = npy.equal(maskin,0)
    # init arrays
    tmp    = npy.ma.ones([N_t, N_j*N_i], dtype='float32')*valmask 
    denflx  = tmp.copy() # Total density flux
    denflxh = tmp.copy() # heat flux contrib
    denflxw = tmp.copy() # E-P contrib
    #
    atmp    = npy.ma.ones([N_t, N_s+1], dtype='float32')*valmask 
    transf  = atmp.copy() # Total tranformation
    transfh = atmp.copy() # Heat flux tranformation
    transfw = atmp.copy() # Water flux tranformation
    areabin = atmp.copy() # surface of bin
    t_heat  = npy.ma.ones((N_t))*valmask # integral heat flux
    t_wafl  = npy.ma.ones((N_t))*valmask # integral E-P
    transfh = maskVal(transfh, valmask)
    transfw = maskVal(transfw, valmask)
    transf  = maskVal(transf , valmask)
    areabin = maskVal(areabin, valmask)
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
        tost = tos.data[t,:]
        sost = sos.data[t,:]
        # Compute density
        rhon = eosNeutral(tost, sost)-1000.
        # Compute buoyancy/density flux as mass fluxes in kg/m2/s (SI unts)
        #  convwf : kg/m2/s=mm/s -> m/s
        convwf = 1.e-3
        denflxh[t,:] = (-alpha(tost,sost)/cpsw(tost,sost,P))*qnet.data[t,:]
        denflxw[t,:] = (rhon+1000.)*betar(tost,sost)*sost*emp.data[t,:]*convwf
        denflx [t,:] = denflxh[t,:] + denflxw[t,:]
        # Transformation (integral of density flux on density outcrops)
        for ks in range(N_s-1):
            # find indices of points in density bin
            idxbin = npy.argwhere( (rhon >= sigrid[ks]) & (rhon < sigrid[ks+1]) )
            transfh[t,ks] = cdu.averager(denflxh[t,idxbin] * areain[idxbin], axis=0, action='sum')/del_s[ks]
            transfw[t,ks] = cdu.averager(denflxw[t,idxbin] * areain[idxbin], axis=0, action='sum')/del_s[ks]
            areabin[t,ks] = cdu.averager(areain[idxbin], axis=0, action='sum')
        # last bin
        idxbin = npy.argwhere( (rhon >= sigrid[N_s-1]))
        transfh[t,N_s] = cdu.averager(denflxh[t,idxbin] * areain[idxbin], axis=0, action='sum')/del_s[ks]
        transfw[t,N_s] = cdu.averager(denflxw[t,idxbin] * areain[idxbin], axis=0, action='sum')/del_s[ks]
        areabin[t,N_s] = cdu.averager(areain[idxbin], axis=0, action='sum')
        # Total transformation
        transf[t,:] = transfh[t,:] + transfw[t,:]        
        # domain integrals
        # heat flux (conv W -> PW)
        convt  = 1.e-15
        t_heat[t] = cdu.averager(denflxh[t,:]*areain, action='sum')*dt*convt
        # fw flux (conv mm -> m and m3/s to Sv)
        convw = 1.e-3*1.e-6
        t_wafl[t] = cdu.averager(denflxw[t,:]*areain, action='sum')*dt*convw
      
    # Reshape i*j back to i,j for output of density flux
    denflxo  = npy.reshape(denflx , (N_t, N_j, N_i))
    denflxho = npy.reshape(denflxh, (N_t, N_j, N_i))
    denflxwo = npy.reshape(denflxw, (N_t, N_j, N_i))
    # Wash mask over variables for transformation
    maskt         = mv.masked_values(tos, valmask).mask
    denflxo.mask  = maskt
    denflxho.mask = maskt
    denflxwo.mask = maskt
    denflxo       = maskVal(denflxo , valmask)
    denflxho      = maskVal(denflxho, valmask)
    denflxwo      = maskVal(denflxwo, valmask)

    maskin       = mv.masked_values(transf, valmask).mask
    transfh._FillValue = valmask
    transfw._FillValue = valmask
    transf._FillValue  = valmask
    transfh.mask = maskin
    transfw.mask = maskin
    transf.mask  = maskin
    transfh      = maskVal(transfh, valmask)
    transfw      = maskVal(transfw, valmask)
    transf       = maskVal(transf , valmask)
       
    #+ create a basins variables (loop on n masks)
    # CPU use
    print
    print ' CPU use', timc.clock() - cpu0

    print ' -> Calculated  denflx, denflxh, denflxw, t_heat, t_wafl'
    if debugp:
        print ' t_heat',t_heat
        print ' t_wafl',t_wafl
    #
    # Output files as netCDF
    # Def variables 
    denFlx  = cdm.createVariable(denflxo , axes = [time, ingrid], id = 'denflux')
    denFlxh = cdm.createVariable(denflxho, axes = [time, ingrid], id = 'hdenflx')
    denFlxw = cdm.createVariable(denflxwo, axes = [time, ingrid], id = 'wdenflx')
    denFlx.long_name   = 'Total density flux'
    denFlx.units       = 'kg/m2/s'
    denFlxh.long_name  = 'Heat density flux'
    denFlxh.units      = 'kg/m2/s'
    denFlxw.long_name  = 'Water density flux'
    denFlxw.units      = 'kg/m2/s'


    outFile_f.write(denFlx)
    outFile_f.write(denFlxh)
    outFile_f.write(denFlxw)

    for i in range(0,len(file_dic)):
        dm = file_dic[i]
        setattr(outFile_f,dm[0],dm[1])
        post_txt = 'Density flux vi surfTransf using delta_sigma = '+str(del_s1)+' and '+str(del_s2)
        setattr(outFile_f, 'Post_processing_history', post_txt)

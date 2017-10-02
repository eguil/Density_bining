#!/home/nillod/anaconda2/bin/python
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 14 21:32:13 2014

Paul J. Durack 14th September 2014

This script computes density bins and indexes T,S values from the vertical z
coordinate onto an approximate neutral density coordinate

Reads in netCDF T(x,y,z,t) and S(x,y,z,t) files and writes
- T or S(x,y,sigma,t)
- D(x,y,sigma,t) (depth of isopycnal)
- Z(x,y,sigma,t) (thickness of isopycnal)

Uses McDougall and Jackett 2005 EOS (IDL routine provided by G. Madec)
Inspired from IDL density bin routines by G. Roullet and G. Madec (1999)
---------------------------------------------------------------------------------
E. Guilyardi - while at LBNL/LLNL -  March 2014

PJD 14 Sep 2014     - Started file
PJD 16 Sep 2014     - Turned off numpy warnings - should be aware of this for final code version
PJD 18 Sep 2014     - Added fixVarUnits function to densityBin
EG  23 Sep 2014     - Clean up and more comments
PJD 16 Oct 2014     - Added getGitInfo,globalAttWrite for metadata writing to outfiles
EG  03 Feb 2015     - Code optimisation (removing loop in persistence)
                    - TODO:
test

@author: durack1
"""

import numpy as npy
from string import replace
import time as timc
from scipy.interpolate import interp1d
from scipy.interpolate._fitpack import _bspleval
import gc,os,resource,timeit
import cdms2 as cdm
import cdutil as cdu
import MV2 as mv
from cdms2 import CdmsRegrid, mvCdmsRegrid
from durolib import fixVarUnits,getGitInfo,globalAttWrite
import ESMP
from libDensity import *

# Turn off numpy warnings
npy.seterr(all='ignore') ; # Cautious use of this turning all error reporting off - shouldn't be an issue as using masked arrays

# Constant values definition

# temporary fix to read grid from CMIP5 IPSL file
GRIDFILE_THETAO = '/prodigfs/project/CMIP5/main/IPSL/IPSL-CM5B-LR/piControl/mon/ocean/Omon/r1i1p1/latest/thetao/thetao_Omon_IPSL-CM5B-LR_piControl_r1i1p1_183001-187912.nc'
GRIDFILE_SO = '/prodigfs/project/CMIP5/main/IPSL/IPSL-CM5B-LR/piControl/mon/ocean/Omon/r1i1p1/latest/so/so_Omon_IPSL-CM5B-LR_piControl_r1i1p1_183001-187912.nc'
GRIDFILE_V = '/prodigfs/project/CMIP5/main/IPSL/IPSL-CM5B-LR/piControl/mon/ocean/Omon/r1i1p1/latest/vo/vo_Omon_IPSL-CM5B-LR_piControl_r1i1p1_183001-187912.nc'

MAX_DEPTH_OCE = 6000. # Ocean max depth
TARGET_GRID_FILENAME = '140807_WOD13_masks.nc' # Target horizonal grid file name for interp
TARGET_GRID_VARNAME = 'basinmask3' # Target horizonal grid var name to use for interp

# Define rho grid with zoom on higher densities
RHO_MIN = 19.
RHO_INT = 26.
RHO_MAX = 28.5
DEL_S1  = 0.2
DEL_S2  = 0.1



def initDensityBin():
    print ''

def densityBin(fileT,fileS,fileV,fileFx,outFile,debug=True,timeint='all',mthout=False):
    # Keep track of time (CPU and elapsed)
    ti0 = timc.clock()
    te0 = timeit.default_timer()

    # CDMS initialisation - netCDF compression
    comp = 1 ; # 0 for no compression
    cdm.setNetcdfShuffleFlag(comp)
    cdm.setNetcdfDeflateFlag(comp)
    cdm.setNetcdfDeflateLevelFlag(comp)
    cdm.setAutoBounds('on')
    # Numpy initialisation
    npy.set_printoptions(precision=2)

    # Determine file name from inputs
    modeln = fileT.split('/')[-1].split('.')[1]

    # Declare and open files for writing too
    #/work/durack1/Shared/data_density/140915/cmip5.ACCESS1-0.historical.r1i1p1.mo.ocn.Omon.density.ver-1.nc
    outFile = replace(outFile,'.mo.','.an.')
    if os.path.isfile(outFile):
        os.remove(outFile)
    if not os.path.exists(os.path.join(*outFile.split('/')[0:-2])):
        os.makedirs(os.path.join(*outFile.split('/')[0:-2]))
    if not os.path.exists(os.path.join(*outFile.split('/')[0:-1])):
        os.makedirs(os.path.join(*outFile.split('/')[0:-1]))
        # Need to convert to shutil - tree create
    outFile_f = cdm.open(outFile,'w')
    if mthout:
        outFileMon = replace(outFile,'.an.','.mo.')
        if os.path.isfile(outFileMon):
            os.remove(outFileMon)
        outFileMon_f = cdm.open(outFileMon,'w') ; # g
    # Size of uncompressed files:
    #  Monthly mean of T,S, thickness and depth on neutral density bins on source grid - IPSL (182x149x61) ~6GB 20yrs
    #  Annual zonal mean of T,S, thick, depth and volume per basin on WOA grid - IPSL 60MB 275yrs
    #  Annual mean persistence variables on WOA grid - IPSL 200MB 150yrs
    #  Annual mean zonal mean of persistence on WOA grid - IPSL 60MB 150yrs

    if debug:
        debug = True
        # CPU analysis
        cpuan = True
    else:
        debug = False
        cpuan = False

    # Open files to read
    ft = cdm.open(fileT)
    fs = cdm.open(fileS)
    fv = cdm.open(fileV)
    # temporary fix to read grid from CMIP5 IPSL file
    ft2 = cdm.open(GRIDFILE_THETAO)
    fs2 = cdm.open(GRIDFILE_SO)
    fv2 = cdm.open(GRIDFILE_V)
    timeax  = ft.getAxis('time')
    # Define temperature and salinity arrays
    thetao_h    = ft2['thetao'] ; # Create variable handle
    so_h        = fs2['so'] ; # Create variable handle
    vo_h        = fv2['vo'] ; # Create variable handle
    tur = timc.clock()
    # Read time and grid
    lon     = thetao_h.getLongitude()
    lat     = thetao_h.getLatitude()
    depth   = thetao_h.getLevel()
    # depth profiles:
    z_zt = depth[:]
    try:
        bounds  = ft2('lev_bnds')
        z_zw = bounds.data[:,0]
    except Exception,err:
        print 'Exception: ',err
        bounds  = depth.getBounds() ; # Work around for BNU-ESM
        z_zw = bounds[:,0]
    max_depth_ocean = MAX_DEPTH_OCE # maximum depth of ocean
    # Horizontal grid
    ingrid  = thetao_h.getGrid()
    # Get grid objects
    axesList = thetao_h.getAxisList()
    # Define dimensions
    lonN    = so_h.shape[3]
    latN    = so_h.shape[2]
    depthN  = so_h.shape[1]
    # Read masking value
    try:
        valmask = so_h.missing_value
        if valmask == None:
            print 'EC-EARTH missing_value fix'
            valmask = 1.e20
    except Exception,err:
        print 'Exception: ',err
        if 'EC-EARTH' == modeln:
            print 'EC-EARTH missing_value fix'
            valmask = 1.e20
    # Test to ensure thetao and so are equivalent sized (times equal)
    if so_h.shape[3] != thetao_h.shape[3] or so_h.shape[2] != thetao_h.shape[2] \
        or so_h.shape[1] != thetao_h.shape[1] or so_h.shape[0] != thetao_h.shape[0]:
        print '** Input variables have different dimensions, exiting..'
        return
    # Test to ensure thetao and vo are equivalent sized (times equal)
    if vo_h.shape[3] != thetao_h.shape[3] or vo_h.shape[2] != thetao_h.shape[2] \
        or vo_h.shape[1] != thetao_h.shape[1] or vo_h.shape[0] != thetao_h.shape[0]:
        print '** Input variables have different dimensions (problem in vo), exiting..'
        return
    #
    thetaoLongName = thetao_h.long_name
    soLongName = so_h.long_name
    soUnits = so_h.units
    del(thetao_h,so_h); gc.collect()

    # Compute level thickness
    lev_thick     = npy.roll(z_zw,-1)-z_zw
    lev_thick[-1] = lev_thick[-2]*.5
    lev_thickt    = npy.swapaxes(mv.reshape(npy.tile(lev_thick,lonN*latN),(lonN*latN,depthN)),0,1)
    
    # Dates to read
    if timeint == 'all':
        tmin = 0
        tmax = timeax.shape[0]
    else:
        tmin = int(timeint.split(',')[0]) - 1
        tmax = tmin + int(timeint.split(',')[1])

    # Target horizonal grid for interp
    msk = Mask(TARGET_GRID_FILENAME, TARGET_GRID_VARNAME)
    tmsk = timc.clock()
    loni    = msk.lonI
    lati    = msk.latI
    Nii     = len(loni)
    Nji     = len(lati)

    # Read cell area and create associated variables
    area = Area(fileFx, loni, lati, {'glob':msk.maski, 'Atl':msk.maskAtl, 'Pac':msk.maskPac, 'Ind':msk.maskInd})
    tarea = timc.clock()

    # Define rho grid with zoom on higher densities
    s_s, s_sax, del_s, N_s = rhonGrid(RHO_MIN, RHO_INT, RHO_MAX, DEL_S1, DEL_S2)
    s_s = npy.tile(s_s, lonN*latN).reshape(lonN*latN,N_s).transpose() # make 3D for matrix computation
    rhoAxis, basinAxis = createAxisRhoBassin(s_sax)
    del(s_sax) ; gc.collect()
    # Create rho axis list
    rhoAxesList             = [axesList[0],rhoAxis,axesList[2],axesList[3]] ; # time, rho, lat, lon
    # Create basin-zonal axes lists
    basinTimeList           = [axesList[0],basinAxis] ; # time, basin
    basinAxesList           = [axesList[0],basinAxis,axesList[2]] ; # time, basin, lat
    basinRhoAxesList        = [axesList[0],basinAxis,rhoAxis,axesList[2]] ; # time, basin, rho, lat

    tinit     = timc.clock()
    # ---------------------
    #  Init density bining
    # ---------------------
    # test point (for debug only)
    itest = 80
    jtest = 30
    ijtest = jtest*lonN + itest

    # Define time read interval (as function of 3D array size)  #TODO: review to optimize
    grdsize = lonN * latN * depthN
    print 'grdsize:',grdsize

    # define number of months in each chunk
    if grdsize > 5.e7:
        tcdel = min(12,tmax) ; # MIROC4h 24 months ~60Gb/50%
    elif grdsize > 2.5e7:
        tcdel = min(24,tmax)
    else:
        tcdel = min(48,tmax)
    tcdel = min(12,tmax) # min(24, tmax) # faster than higher tcdel ?
    nyrtc = tcdel/12
    tcmax = (tmax-tmin)/tcdel ; # number of time chunks
    print ' ==> model:', modeln,' (grid size:', grdsize,')'
    print ' ==> time interval: ', tmin, tmax - 1
    print ' ==> size of time chunk, number of time chunks (memory optimization) :', tcdel, tcmax

    # Preallocate masked arrays on target grid
    # Global arrays on target grid
    depthBini   = npy.ma.ones([nyrtc, N_s+1, Nji, Nii], dtype='float32')*valmask
    thickBini,x1Bini,x2Bini,x3Bini = [npy.ma.ones(npy.ma.shape(depthBini)) for _ in range(4)]
    # Basin zonal on target grid
    depthBinia,thickBinia,x1Binia,x2Binia,x3Binia,depthBinip,thickBinip,\
    x1Binip,x2Binip,x3Binip,depthBinii,thickBinii,x1Binii,x2Binii,x3Binii = [npy.ma.ones(npy.shape(depthBini)) for _ in range(15)]
    # Persistence arrays on original grid
    persist     = npy.ma.ones([nyrtc, N_s+1, latN, lonN], dtype='float32')*valmask
    persisti,persistia,persistip,persistii,persistv = [npy.ma.ones(npy.shape(depthBini)) for _ in range(5)]
    # Persistence arrays on target grid
    persistm    = npy.ma.ones([nyrtc, Nji, Nii], dtype='float32')*valmask
    ptopdepthi,ptopsigmai,ptoptempi,ptopsalti,ptophvmi = [npy.ma.ones(npy.shape(persistm)) for _ in range(5)]
    # Basin zonal on target grid
    ptopdepthia,ptopsigmaia,ptoptempia,ptopsaltia,ptophvmia,ptopdepthip,ptopsigmaip,\
    ptoptempip,ptopsaltip,ptophvmip,ptopdepthii,ptopsigmaii,ptoptempii,ptopsaltii,ptophvmii = [npy.ma.ones(npy.shape(persistm)) for _ in range(15)]
    # Volume/thetao/so of persistent ocean
    volpersist = npy.ma.ones([nyrtc], dtype='float32')*valmask
    volpersista,volpersistp,volpersisti, \
      tempersist,tempersista,tempersistp,tempersisti, \
      salpersist,salpersista,salpersistp,salpersisti, \
      hvmpersist,hvmpersista,hvmpersistp,hvmpersisti  = [npy.ma.ones(npy.shape(volpersist)) for _ in range(15)]

    ######################################### START CHUNK #########################################
    # Interpolation init (regrid)
    ESMP.ESMP_Initialize()
    regridObj = CdmsRegrid(ingrid,msk.outgrid,depthBini.dtype,missing=valmask,regridMethod='distwgt',regridTool='esmf', coordSys='deg', diag = {},periodicity=1)
    tintrp     = timc.clock()
    # testing
    voltotij0 = npy.ma.ones([latN*lonN], dtype='float32')*0.
    temtotij0 = npy.ma.ones([latN*lonN], dtype='float32')*0.
    saltotij0 = npy.ma.ones([latN*lonN], dtype='float32')*0.
    hvmtotij0 = npy.ma.ones([latN*lonN], dtype='float32')*0.
    tin1     = timc.clock()
    
    if cpuan:
        print ' '
    # -----------------------------------------
    #  Density bining loop (on time chunks tc)
    # -----------------------------------------
    for tc in range(tcmax):
        tuc     = timc.clock()
        # output arrays for each chunk
        tmp         = npy.ma.ones([tcdel, N_s+1, latN*lonN], dtype='float32')*valmask
        depth_bin   = tmp.copy() ; depth_bin   = maskVal(depth_bin, valmask)
        thick_bin   = tmp.copy() ; thick_bin   = maskVal(thick_bin, valmask)
        x1_bin      = tmp.copy() ; x1_bin      = maskVal(x1_bin, valmask)
        x2_bin      = tmp.copy() ; x2_bin      = maskVal(x2_bin, valmask)
        x3_bin      = tmp.copy() ; x3_bin      = maskVal(x3_bin, valmask) #x3 is added for the meridional transport
        del(tmp) ; gc.collect()
        # read tcdel month by tcdel month to optimise memory
        trmin   = tmin + tc*tcdel ; # define as function of tc and tcdel
        trmax   = tmin + (tc+1)*tcdel ; # define as function of tc and tcdel
        print ' --> time chunk (bounds) = ',tc+1, '/',tcmax,' (',trmin,trmax-1,')', modeln
        thetao  = ft('thetao', time = slice(trmin,trmax))
        so      = fs('so'    , time = slice(trmin,trmax))
        vo      = fv('vo'    , time = slice(trmin,trmax))
        time    = thetao.getTime()
        testval = valmask
        # Check for missing_value/mask
        if ( 'missing_value' not in thetao.attributes.keys() and modeln == 'EC-EARTH' ) \
           or (modeln == 'MIROC4h' ):
            print 'trigger mask fix - EC-EARTH/MIROC4h'
            so = mv.masked_equal(so,0.)
            print so.count()
            so.data[:] = so.filled(valmask)
            thetao.mask = so.mask
            thetao.data[:] = thetao.filled(valmask)
        # Define rho output axis
        rhoAxesList[0]  = time ; # replace time axis
        
        # Test variable units
        [so,soFixed] = fixVarUnits(so,'so',True)#,'logfile.txt')
        #if soFixed:
        #    print '     so:     units corrected'
        [thetao,thetaoFixed] = fixVarUnits(thetao,'thetao',True)#,'logfile.txt')
        # 17/08/17: something should be added for vo I guess (which should be in m/s) but I am not sure how yet. 

        turd = timc.clock()
        # Compute neutral density
        rhon = eosNeutral(thetao,so)-1000.
        turr = timc.clock()

        # reorganise i,j dims in single dimension data (speeds up loops)
        thetao  = mv.reshape(thetao,(tcdel, depthN, lonN*latN))
        so      = mv.reshape(so    ,(tcdel, depthN, lonN*latN))
        rhon    = mv.reshape(rhon  ,(tcdel, depthN, lonN*latN))
        vo      = mv.reshape(vo  ,(tcdel, depthN, lonN*latN))
        #print 'thetao.shape:',thetao.shape
        if debug and tc < 0 :
            print ' thetao :',thetao.data[0,:,ijtest]
            print ' thetao :',thetao[0,:,ijtest]
            print ' so     :',so.data    [0,:,ijtest]
            print ' so     :',so[0,:,ijtest]
            print ' vo     :',vo.data    [0,:,ijtest]
            print ' vo     :',vo[0,:,ijtest]
        # Reset output arrays to missing for binned fields
        depth_Bin,thick_bin,x1_bin,x2_bin,x3_bin = [npy.ma.ones([tcdel, N_s+1, latN*lonN])*valmask for _ in range(5)]

        tucz0     = timc.clock()
        if cpuan:
            cpu1 = 0.
            cpu2 = 0.
            cpu3 = 0.
            cpu4 = 0.
            cpu40 = 0.
            cpu5 = 0.
        # Loop on time within chunk tc
        for t in range(trmax-trmin):
            tcpu0 = timc.clock()
            
            # find bottom level at each lat/lon point
            # x1 contents on vertical (not yet implemented - may be done to ensure conservation)
            x1_content  = thetao.data[t]
            x2_content  = so.data[t]
            x3_content  = npy.array(vo.data[t]*lev_thickt) # x3 doit representer la vitesse integree sur l epaisseur de la grille. 

            # Find indexes of masked points
            vmask_3D    = mv.masked_values(so.data[t],testval).mask ; # Returns boolean
            nomask      = npy.equal(vmask_3D[0],0) ; # Returns boolean
            # Check integrals on source z coordinate grid
            if debug and t == 0:
                voltotij0 = npy.sum(lev_thickt*(1-vmask_3D[:,:]), axis=0)
                temtotij0 = npy.sum(lev_thickt*(1-vmask_3D[:,:])*x1_content[:,:], axis=0)
                saltotij0 = npy.sum(lev_thickt*(1-vmask_3D[:,:])*x2_content[:,:], axis=0)
                hvmtotij0 = npy.sum(lev_thickt*(1-vmask_3D[:,:])*x3_content[:,:], axis=0)
                voltot = npy.sum(voltotij0*mv.reshape(area.area,lonN*latN))
                temtot = npy.sum(temtotij0*mv.reshape(area.area,lonN*latN))/voltot
                saltot = npy.sum(saltotij0*mv.reshape(area.area,lonN*latN))/voltot
                hvmtot = npy.sum(hvmtotij0*mv.reshape(area.area,lonN*latN))/voltot #NL: line added!
                print '  Total volume in z coordinates source grid (ref = 1.33 e+18) : ', voltot
                print '  Mean Temp./Salinity in z coordinates source grid            : ', temtot, saltot
                print '  Mean Meridional TRANSPORT in z coordinates source grid       : ', hvmtot
                # 17/08/17 JM: consider checking the total integrated transport rather than mean velocity. 

            # init arrays for this time chunk
            z_s,c1_s,c2_s,c3_s,t_s  = [npy.ma.ones((N_s+1, lonN*latN))*valmask for _ in range(5)]
            szmin,szmax,delta_rho   = [npy.ma.ones(lonN*latN)*valmask for _ in range(3)]
            i_min,i_max             = [npy.ma.zeros(lonN*latN) for _ in range(2)]
            tcpu1 = timc.clock()
            # find bottom level at each lat/lon point
            i_bottom                = vmask_3D.argmax(axis=0)-1
            z_s [N_s, nomask]   = z_zw[i_bottom[nomask]+1] ; # Cell depth limit
            c1_s[N_s, nomask]   = x1_content[depthN-1,nomask] ; # Cell bottom temperature/salinity
            c2_s[N_s, nomask]   = x2_content[depthN-1,nomask] ; # Cell bottom tempi_profilerature/salinity
            c3_s[N_s, nomask]   = x3_content[depthN-1,nomask] ; # Cell bottom meridional transport
            # init arrays as a function of depth = f(z)
            s_z     = rhon.data[t]
            c1_z    = x1_content
            c2_z    = x2_content
            c3_z    = x3_content
            # Extract a strictly increasing sub-profile
            i_min[nomask]           = s_z.argmin(axis=0)[nomask]
            i_max[nomask]           = s_z.argmax(axis=0)[nomask]-1 # why -1 ?
            i_min[i_min > i_max]    = i_max[i_min > i_max]
            # Test on bottom minus surface stratification to check that it is larger than delta_rho
            delta_rho[nomask]           = s_z[i_bottom[nomask],nomask] - s_z[0,nomask]
            i_min[delta_rho < DEL_S1]   = 0
            i_max[delta_rho < DEL_S1]   = i_bottom[delta_rho < DEL_S1]

            # General case
            # find min/max of density for each z profile
            for i in range(lonN*latN):
                if nomask[i]:
                    szmin[i] = s_z[int(i_min[i]),i]
                    szmax[i] = s_z[int(i_max[i]),i]
                else:
                    szmin[i] = 0.
                    szmax[i] = RHO_MAX+10.
            tcpu2 = timc.clock()
            # Find indices between density min and density max
            #
            # Construct arrays of szm/c1m/c2m = s_z[i_min[i]:i_max[i],i] and valmask otherwise
            # same for zzm from z_zt
            szm,zzm,c1m,c2m,c3m  = [npy.ma.ones(s_z.shape)*valmask for _ in range(5)]

            for k in range(depthN):
                k_ind = i_min*1.; k_ind[:] = valmask
                k_ind = npy.argwhere( (k >= i_min) & (k <= i_max))
                szm[k,k_ind] = s_z [k,k_ind]
                c1m[k,k_ind] = c1_z[k,k_ind]
                c2m[k,k_ind] = c2_z[k,k_ind]
                c3m[k,k_ind] = c3_z[k,k_ind]
                zzm[k,k_ind] = z_zt[k]

            # interpolate depth(z) (=z_zt) to depth(s) at s_s densities (=z_s) using density(z) (=s_z)
            # TODO: use ESMF ? # TODO check that interp in linear or/and stabilise column as post-pro
            tcpu3 = timc.clock()
            for i in range(lonN*latN):
                if nomask[i]:
                    z_s [0:N_s,i] = npy.interp(s_s[:,i], szm[:,i], zzm[:,i], right = valmask) ; # depth - consider spline
                    c1_s[0:N_s,i] = npy.interp(z_s[0:N_s,i], zzm[:,i], c1m[:,i], right = valmask) ; # thetao
                    c2_s[0:N_s,i] = npy.interp(z_s[0:N_s,i], zzm[:,i], c2m[:,i], right = valmask) ; # so
                    c3_s[0:N_s,i] = npy.interp(z_s[0:N_s,i], zzm[:,i], c3m[:,i], right = valmask) ; # hvo
            # if level in s_s has lower density than surface, isopycnal is put at surface (z_s = 0)
            tcpu40 = timc.clock()

            # if level of s_s has higher density than bottom density,
            # isopycnal is set to bottom (z_s = z_zw[i_bottom])
            # TODO:  add half level to depth to ensure thickness integral conservation
            inds = npy.argwhere(s_s > szmax).transpose()
            z_s [inds[0],inds[1]] = z_s[N_s-1,inds[1]]
            c1_s[inds[0],inds[1]] = c1_s[N_s-1,inds[1]]
            c2_s[inds[0],inds[1]] = c2_s[N_s-1,inds[1]]
            c3_s[inds[0],inds[1]] = c3_s[N_s-1,inds[1]]
            tcpu4 = timc.clock()
            # Thickness of isopycnal
            t_s [0,:] = 0.
            t_s [1:N_s,:] = z_s[1:N_s,:]-z_s[0:N_s-1,:]
            # Use thickness of isopycnal (less than zero) to create masked point for all bined arrays
            inds = npy.argwhere( (t_s <= 0.) ^ (t_s >= max_depth_ocean)).transpose()
            t_s [inds[0],inds[1]] = valmask
            z_s [inds[0],inds[1]] = valmask
            c1_s[inds[0],inds[1]] = valmask
            c2_s[inds[0],inds[1]] = valmask
            c3_s[inds[0],inds[1]] = valmask
            if debug and t < 0: #t == 0:
                i = ijtest
                print
                print ' density target array s_s[i]'
                print s_s[:,i]
                print ' density profile on Z grid szm[i]'
                print szm[:,i]
                print ' depth profile on Z grid zzm[i]'
                print zzm[:,i]
                print ' depth profile on rhon target grid z_s[i]'
                print z_s[:,i]
                print 'tc = ',tc
                print ' thickness profile on rhon grid t_s[i]'
                print t_s[:,i]
                print ' bined temperature profile on rhon grid c1_s[i]'
                print c1_s[:,i]
            # assign to final arrays
            depth_bin[t,:,:] = z_s
            thick_bin[t,:,:] = t_s
            x1_bin[t,:,:]    = c1_s
            x2_bin[t,:,:]    = c2_s
            x3_bin[t,:,:]    = c3_s
            # CPU analysis
            tcpu5 = timc.clock()
            if cpuan:
                cpu1 = cpu1 + tcpu1 - tcpu0
                cpu2 = cpu2 + tcpu2 - tcpu1
                cpu3 = cpu3 + tcpu3 - tcpu2
                cpu4 = cpu4 + tcpu40 - tcpu3
                cpu40 = cpu40 + tcpu4 - tcpu40
                cpu5 = cpu5 + tcpu5 - tcpu4
        #
        # end of loop on t <===
        #
        # CPU analysis
        if cpuan:
            print ' Bining CPU analysis tc/tcdel = ',tc,tcdel
            print '    average cpu1  = ',cpu1/float(tcdel)
            print '    average cpu2  = ',cpu2/float(tcdel)
            print '    average cpu3  = ',cpu3/float(tcdel)
            print '    average cpu4  = ',cpu4/float(tcdel)
            print '    average cpu40 = ',cpu40/float(tcdel)
            print '    average cpu5  = ',cpu5/float(tcdel)
            print '    CPU read T/S  = ',turd-tuc
            print '    CPU comp. rho = ',turr-turd

        ticz0 = timc.clock()
        # Free memory
        del(rhon, x1_content, x2_content, x3_content, vmask_3D, szm, zzm, c1m, c2m, c3m, z_s, c1_s, c2_s, c3_s, t_s, inds, c1_z, c2_z, c3_z) ; gc.collect()
        # Wash mask (from temp) over variables
        maskb          = mv.masked_values(x1_bin, valmask).mask
        depth_bin.mask = maskb
        x1_bin.mask    = maskb
        x2_bin.mask    = maskb
        x3_bin.mask    = maskb
        maskt          = mv.masked_values(thick_bin, valmask).mask
        thick_bin.mask = maskt
        depth_bin      = maskVal(depth_bin, valmask)
        thick_bin      = maskVal(thick_bin, valmask)
        x1_bin         = maskVal(x1_bin, valmask)
        x2_bin         = maskVal(x2_bin, valmask)
        x3_bin         = maskVal(x3_bin, valmask)

        # 17/08/17: JM:  specificity of vmask to be checked.
        #                1st approx: v on T-points.
        
        # Reshape i*j back to i,j
        depth_bin = npy.ma.reshape(depth_bin, (tcdel, N_s+1, latN, lonN))
        thick_bin = npy.ma.reshape(thick_bin, (tcdel, N_s+1, latN, lonN))
        x1_bin    = npy.ma.reshape(x1_bin,    (tcdel, N_s+1, latN, lonN))
        x2_bin    = npy.ma.reshape(x2_bin,    (tcdel, N_s+1, latN, lonN))
        x3_bin    = npy.ma.reshape(x3_bin,    (tcdel, N_s+1, latN, lonN))

        if debug and (tc < 0):
            # test write
            i = itest
            j = jtest
            print 'test point',i,j, area.area[j,i]
            try:
                print 'lon,lat',lon[j,i],lat[j,i]
            except Exception,err:
                print 'Exception: ',err
                print 'lon,lat',lon[i],lat[j]
            print 'depth_bin', depth_bin[0,:,j,i]
            print 'thick_bin', thick_bin[0,:,j,i]
            print 'x1_bin', x1_bin[0,:,j,i]
            print 'x2_bin', x2_bin[0,:,j,i]
            print 'x3_bin', x3_bin[0,:,j,i]
        if debug and tc == 0:
            # Check integrals/mean on source density grid
            print voltotij0[ijtest], temtotij0[ijtest],saltotij0[ijtest],hvmtotij0[ijtest]
            voltotij0 = npy.sum(npy.ma.reshape(thick_bin,(tcdel, N_s+1, latN*lonN)).data[tc,:,:]*(1-npy.ma.reshape(thick_bin,(tcdel, N_s+1, latN*lonN)).mask[tc,:,:]), axis=0)
            temtotij0 = npy.sum(npy.ma.reshape(thick_bin,(tcdel, N_s+1, latN*lonN)).data[tc,:,:]*npy.ma.reshape(x1_bin,(tcdel, N_s+1, latN*lonN)).data[tc,:,:]*(1-npy.ma.reshape(thick_bin,(tcdel, N_s+1, latN*lonN)).mask[tc,:,:]), axis=0)
            saltotij0 = npy.sum(npy.ma.reshape(thick_bin,(tcdel, N_s+1, latN*lonN)).data[tc,:,:]*npy.ma.reshape(x2_bin,(tcdel, N_s+1, latN*lonN)).data[tc,:,:]*(1-npy.ma.reshape(thick_bin,(tcdel, N_s+1, latN*lonN)).mask[tc,:,:]), axis=0)
            hvmtotij0 = npy.sum(npy.ma.reshape(thick_bin,(tcdel, N_s+1, latN*lonN)).data[tc,:,:]*npy.ma.reshape(x3_bin,(tcdel, N_s+1, latN*lonN)).data[tc,:,:]*(1-npy.ma.reshape(thick_bin,(tcdel, N_s+1, latN*lonN)).mask[tc,:,:]), axis=0)
            voltot = npy.sum(voltotij0*npy.ma.reshape(area.area,lonN*latN))
            temtot = npy.sum(temtotij0*npy.ma.reshape(area.area,lonN*latN))/voltot
            saltot = npy.sum(saltotij0*npy.ma.reshape(area.area,lonN*latN))/voltot
            hvmtot = npy.sum(hvmtotij0*npy.ma.reshape(area.area,lonN*latN))/voltot
            print voltotij0[ijtest], temtotij0[ijtest],saltotij0[ijtest],hvmtotij0[ijtest]
            print '  Total volume in rho coordinates source grid (ref = 1.33 e+18) : ', voltot
            print '  Mean Temp./Salinity in rho coordinates source grid            : ', temtot, saltot
            print '  Mean Meridional velocity in z coordinates source grid       : ', hvmtot
        #
        # Output files as netCDF
        # Def variables
        depthBin = cdm.createVariable(depth_bin, axes = rhoAxesList, id = 'isondepth')
        thickBin = cdm.createVariable(thick_bin, axes = rhoAxesList, id = 'isonthick')
        x1Bin    = cdm.createVariable(x1_bin   , axes = rhoAxesList, id = 'thetao')
        x2Bin    = cdm.createVariable(x2_bin   , axes = rhoAxesList, id = 'so')
        x3Bin    = cdm.createVariable(x3_bin   , axes = rhoAxesList, id = 'hvo')
        #
        del (depth_bin,thick_bin,x1_bin,x2_bin,x3_bin) ; gc.collect()
        if mthout:
            if tc == 0:
                depthBin.long_name  = 'Depth of isopycnal'
                depthBin.units      = 'm'
                thickBin.long_name  = 'Thickness of isopycnal'
                thickBin.units      = 'm'
                x1Bin.long_name     = thetaoLongName
                x1Bin.units         = 'C'
                x2Bin.long_name     = soLongName
                x2Bin.units         = soUnits
                x3Bin.long_name     = 'Oceanic meridional transport h*vo'
                x3Bin.units         = 'm/s * m'
                
                outFileMon_f.write(area.area.astype('float32')) ; # Added area so isonvol can be computed

        # -------------------------------------------------------------
        #  Compute annual mean, persistence, make zonal mean and write
        # -------------------------------------------------------------
        ticz = timc.clock()
        if tcdel >= 12:
            # Annual mean
            dym  = npy.ma.reshape(depthBin, (nyrtc, 12, N_s+1, latN, lonN))
            tym  = npy.ma.reshape(thickBin, (nyrtc, 12, N_s+1, latN, lonN))
            x1ym = npy.ma.reshape(x1Bin,    (nyrtc, 12, N_s+1, latN, lonN))
            x2ym = npy.ma.reshape(x2Bin,    (nyrtc, 12, N_s+1, latN, lonN))
            x3ym = npy.ma.reshape(x3Bin,    (nyrtc, 12, N_s+1, latN, lonN))

            # this divided the CPU by 5 for annual mean
            validPoints = dym/dym
            validMonths = npy.ma.sum(validPoints , axis=1)
            dy  = npy.ma.sum(dym , axis=1)/validMonths
            ty  = npy.ma.sum(tym , axis=1)/validMonths
            x1y = npy.ma.sum(x1ym, axis=1)/validMonths
            x2y = npy.ma.sum(x2ym, axis=1)/validMonths
            x3y = npy.ma.sum(x3ym, axis=1)/validMonths

            del (dym,tym,x1ym,x2ym,x3ym) ; gc.collect()
            # create annual time axis
            timeyr          = cdm.createAxis(dy.getAxis(0))
            timeyr.id       = 'time'
            timeyr.units    = time.units
            timeyr.designateTime()
            rhoAxesList[0]  = timeyr ; # replace time axis

            dy   = cdm.createVariable(dy , axes = rhoAxesList, id = 'isondy')
            ty   = cdm.createVariable(ty , axes = rhoAxesList, id = 'isonty')
            x1y  = cdm.createVariable(x1y, axes = rhoAxesList, id = 'isonx1y')
            x2y  = cdm.createVariable(x2y, axes = rhoAxesList, id = 'isonx2y')
            x3y  = cdm.createVariable(x3y, axes = rhoAxesList, id = 'isonx3y')

            toz = timc.clock()

            # Interpolate onto common grid
            for t in range(nyrtc):
                for ks in range(N_s+1):
                    # Global
                    depthBini[t,ks,:,:]         = regridObj(dy [t,ks,:,:])
                    thickBini[t,ks,:,:]         = regridObj(ty [t,ks,:,:])
                    x1Bini[t,ks,:,:]            = regridObj(x1y[t,ks,:,:])
                    x2Bini[t,ks,:,:]            = regridObj(x2y[t,ks,:,:])
                    x3Bini[t,ks,:,:]            = regridObj(x3y[t,ks,:,:])
                    depthBini[t,ks,:,:].mask    = msk.maski
                    thickBini[t,ks,:,:].mask    = msk.maski
                    x1Bini[t,ks,:,:].mask       = msk.maski
                    x2Bini[t,ks,:,:].mask       = msk.maski
                    x3Bini[t,ks,:,:].mask       = msk.maski
                    # Atl
                    depthBinia[t,ks,:,:]        = depthBini[t,ks,:,:]*1.
                    thickBinia[t,ks,:,:]        = thickBini[t,ks,:,:]*1.
                    x1Binia[t,ks,:,:]           = x1Bini[t,ks,:,:]*1.
                    x2Binia[t,ks,:,:]           = x2Bini[t,ks,:,:]*1.
                    x3Binia[t,ks,:,:]           = x3Bini[t,ks,:,:]*1.
                    depthBinia[t,ks,:,:].mask   = msk.maskAtl
                    thickBinia[t,ks,:,:].mask   = msk.maskAtl
                    x1Binia[t,ks,:,:].mask      = msk.maskAtl
                    x2Binia[t,ks,:,:].mask      = msk.maskAtl
                    x3Binia[t,ks,:,:].mask      = msk.maskAtl
                    # Pac
                    depthBinip[t,ks,:,:]        = depthBini[t,ks,:,:]*1.
                    thickBinip[t,ks,:,:]        = thickBini[t,ks,:,:]*1.
                    x1Binip[t,ks,:,:]           = x1Bini[t,ks,:,:]*1.
                    x2Binip[t,ks,:,:]           = x2Bini[t,ks,:,:]*1.
                    x3Binip[t,ks,:,:]           = x3Bini[t,ks,:,:]*1.
                    depthBinip[t,ks,:,:].mask   = msk.maskPac
                    thickBinip[t,ks,:,:].mask   = msk.maskPac
                    x1Binip[t,ks,:,:].mask      = msk.maskPac
                    x2Binip[t,ks,:,:].mask      = msk.maskPac
                    x3Binip[t,ks,:,:].mask      = msk.maskPac
                    # Ind
                    depthBinii[t,ks,:,:]        = depthBini[t,ks,:,:]*1.
                    thickBinii[t,ks,:,:]        = thickBini[t,ks,:,:]*1.
                    x1Binii[t,ks,:,:]           = x1Bini[t,ks,:,:]*1.
                    x2Binii[t,ks,:,:]           = x2Bini[t,ks,:,:]*1.
                    x3Binii[t,ks,:,:]           = x3Bini[t,ks,:,:]*1.
                    depthBinii[t,ks,:,:].mask   = msk.maskInd
                    thickBinii[t,ks,:,:].mask   = msk.maskInd
                    x1Binii[t,ks,:,:].mask      = msk.maskInd
                    x2Binii[t,ks,:,:].mask      = msk.maskInd
                    x3Binii[t,ks,:,:].mask      = msk.maskInd
            # Free memory
            del(dy, ty, x1y, x2y, x3y); gc.collect()

            # Global
            depthBini   = maskVal(depthBini, valmask)
            thickBini   = maskVal(thickBini, valmask)
            x1Bini      = maskVal(x1Bini, valmask)
            x2Bini      = maskVal(x2Bini, valmask)
            x3Bini      = maskVal(x3Bini, valmask)
            # Atl
            depthBinia  = maskVal(depthBinia, valmask)
            thickBinia  = maskVal(thickBinia, valmask)
            x1Binia     = maskVal(x1Binia, valmask)
            x2Binia     = maskVal(x2Binia, valmask)
            x3Binia     = maskVal(x3Binia, valmask)
            # Pac
            depthBinip  = maskVal(depthBinip, valmask)
            thickBinip  = maskVal(thickBinip, valmask)
            x1Binip     = maskVal(x1Binip, valmask)
            x2Binip     = maskVal(x2Binip, valmask)
            x3Binip     = maskVal(x3Binip, valmask)
            # Ind
            depthBinii  = maskVal(depthBinii, valmask)
            thickBinii  = maskVal(thickBinii, valmask)
            x1Binii     = maskVal(x1Binii, valmask)
            x2Binii     = maskVal(x2Binii, valmask)
            x3Binii     = maskVal(x3Binii, valmask)

            depthbini  = cdm.createVariable(depthBini,  axes = [timeyr, rhoAxis, lati, loni], id = 'isondepthg')
            thickbini  = cdm.createVariable(thickBini,  axes = [timeyr, rhoAxis, lati, loni], id = 'isonthickg')
            x1bini     = cdm.createVariable(x1Bini   ,  axes = [timeyr, rhoAxis, lati, loni], id = 'thetaog')
            x2bini     = cdm.createVariable(x2Bini   ,  axes = [timeyr, rhoAxis, lati, loni], id = 'sog')
            x3bini     = cdm.createVariable(x3Bini   ,  axes = [timeyr, rhoAxis, lati, loni], id = 'hvoog')
            if tc == 0:
                depthbini.long_name  = 'Depth of isopycnal'
                depthbini.units      = 'm'
                thickbini.long_name  = 'Thickness of isopycnal'
                thickbini.units      = 'm'
                x1bini.long_name     = thetaoLongName
                x1bini.units         = 'C'
                x2bini.long_name     = soLongName
                x2bini.units         = soUnits
                x3bini.long_name     = 'Meridional Oceanic Transport'
                x3bini.units         = 'm/s * m'

            tozi = timc.clock()

            # Compute zonal mean
            # Global
            depthBinz   = cdu.averager(depthBini,   axis = 3)
            thickBinz   = cdu.averager(thickBini,   axis = 3)
            x1Binz      = cdu.averager(x1Bini,      axis = 3)
            x2Binz      = cdu.averager(x2Bini,      axis = 3)
            x3Binz      = cdu.averager(x3Bini,      axis = 3)
            # Atl
            depthBinza  = cdu.averager(depthBinia,  axis = 3)
            thickBinza  = cdu.averager(thickBinia,  axis = 3)
            x1Binza     = cdu.averager(x1Binia,     axis = 3)
            x2Binza     = cdu.averager(x2Binia,     axis = 3)
            x3Binza     = cdu.averager(x3Binia,     axis = 3)
            # Pac
            depthBinzp  = cdu.averager(depthBinip,  axis = 3)
            thickBinzp  = cdu.averager(thickBinip,  axis = 3)
            x1Binzp     = cdu.averager(x1Binip,     axis = 3)
            x2Binzp     = cdu.averager(x2Binip,     axis = 3)
            x3Binzp     = cdu.averager(x3Binip,     axis = 3)
            # Ind
            depthBinzi  = cdu.averager(depthBinii,  axis = 3)
            thickBinzi  = cdu.averager(thickBinii,  axis = 3)
            x1Binzi     = cdu.averager(x1Binii,     axis = 3)
            x2Binzi     = cdu.averager(x2Binii,     axis = 3)
            x3Binzi     = cdu.averager(x3Binii,     axis = 3)
            # Compute volume of isopycnals
            volBinz     = thickBinz  * area.areazt
            volBinza    = thickBinza * area.areazta
            volBinzp    = thickBinzp * area.areaztp
            volBinzi    = thickBinzi * area.areazti

            toziz = timc.clock()

            # Compute annual persistence of isopycnal bins (from their thickness): 'persist' array
            #  = percentage of time bin is occupied during each year (annual bowl if % < 100)
            for t in range(nyrtc):
                tpe0 = timc.clock()
                idxvm = npy.ma.ones([12, N_s+1, latN, lonN], dtype='float32')*valmask
                inim = t*12
                finm = t*12 + 12
                idxvm = 1-mv.masked_values(thickBin[inim:finm,:,:,:], valmask).mask
                persist[t,:,:,:] = cdu.averager(idxvm, axis = 0) * 100.
                #persist[t,:,:,:] = npy.ma.sum(idxvm, axis = 0)/12. * 100. # numpy version same CPU
                # Shallowest persistent ocean index: p_top (2D)
                maskp = persist[t,:,:,:]*1. ; maskp[...] = valmask
                maskp = mv.masked_values(persist[t,:,:,:] >= 99., 1.).mask
                maskp = npy.ma.reshape(maskp, (N_s+1, latN*lonN))
                p_top = maskp.argmax(axis=0)
                # Define properties on bowl (= shallowest persistent ocean)
                ptopdepth = npy.ma.ones([latN*lonN], dtype='float32')*valmask
                ptopsigma,ptoptemp,ptopsalt = [npy.ma.ones(npy.shape(ptopdepth)) for _ in range(3)]
                tpe1 = timc.clock()
                # Creat array of 1 on bowl and 0 elsewhere
                maskp = (maskp-npy.roll(maskp,1,axis=0))*maskp
                depthBintmp = npy.ma.reshape(depthBin[t,...],(N_s+1, latN*lonN))
                x1Bintmp    = npy.ma.reshape(x1Bin[t,...],(N_s+1, latN*lonN))
                x2Bintmp    = npy.ma.reshape(x2Bin[t,...],(N_s+1, latN*lonN))
                x3Bintmp    = npy.ma.reshape(x3Bin[t,...],(N_s+1, latN*lonN))
                ptopdepth   = cdu.averager(depthBintmp*maskp,axis=0,action='sum')
                ptoptemp    = cdu.averager(x1Bintmp*maskp,axis=0,action='sum')
                ptopsalt    = cdu.averager(x2Bintmp*maskp,axis=0,action='sum')
                ptophvm    = cdu.averager(x3Bintmp*maskp,axis=0,action='sum')

                del (depthBintmp,x1Bintmp,x2Bintmp,x3Bintmp); gc.collect()
                tpe2 = timc.clock()


                ptopsigma = ptopdepth*0. + rhoAxis[p_top] # to keep mask of ptopdepth
                ptopdepth = npy.ma.reshape(ptopdepth, (latN, lonN))
                ptopsigma = npy.ma.reshape(ptopsigma, (latN, lonN))
                ptoptemp  = npy.ma.reshape(ptoptemp , (latN, lonN))
                ptopsalt  = npy.ma.reshape(ptopsalt , (latN, lonN))
                ptophvm  = npy.ma.reshape(ptophvm , (latN, lonN))

                # Create variables to attribute right axis for zonal mean
                ptopdepth = cdm.createVariable(ptopdepth, axes = [ingrid], id = 'ptopdepth')
                ptopsigma = cdm.createVariable(ptopsigma, axes = [ingrid], id = 'ptopsigma')
                ptoptemp  = cdm.createVariable(ptoptemp , axes = [ingrid], id = 'ptopthetao')
                ptopsalt  = cdm.createVariable(ptopsalt , axes = [ingrid], id = 'ptopso')
                ptophvm  = cdm.createVariable(ptophvm , axes = [ingrid], id = 'ptophvm')

                # Mask persist where value is zero
                persist._FillValue = valmask
                persist = mv.masked_where(persist <= 1.e-6, persist)
                persbin = cdm.createVariable(persist, axes = rhoAxesList, id = 'isonpers')

                # Interpolate to target grid and create basin variables
                #
                tpe3 = timc.clock()
                for ks in range(N_s+1):
                    persisti [t,ks,:,:]         = regridObj(persbin[t,ks,:,:])
                    persisti [t,ks,:,:].mask    = msk.maski
                    persistia[t,ks,:,:]         = persisti[t,ks,:,:]*1.
                    persistia[t,ks,:,:].mask    = msk.maskAtl
                    persistip[t,ks,:,:]         = persisti[t,ks,:,:]*1.
                    persistip[t,ks,:,:].mask    = msk.maskPac
                    persistii[t,ks,:,:]         = persisti[t,ks,:,:]*1.
                    persistii[t,ks,:,:].mask    = msk.maskInd

                persisti    = maskVal(persisti,  valmask)
                persistia   = maskVal(persistia, valmask)
                persistip   = maskVal(persistip, valmask)
                persistii   = maskVal(persistii, valmask)
                tpe4 = timc.clock()
                # Compute zonal mean (2D)
                persistiz   = cdu.averager(persisti,  axis = 3)
                persistiza  = cdu.averager(persistia, axis = 3)
                persistizp  = cdu.averager(persistip, axis = 3)
                persistizi  = cdu.averager(persistii, axis = 3)
                # Persistence * thickness (used to compute % of column that is persistent - see below)
                persistv[t,:,:,:]           = persisti[t,:,:,:] * thickBini[t,:,:,:]
                persistv                    = maskVal(persistv, valmask)
                tpe5 = timc.clock()

                # Depth, temperature and salinity on bowl (i.e. at shallowest persistent ocean) (2D)
                ptopdepthi[t,:,:]           = regridObj(ptopdepth)
                ptopsigmai[t,:,:]           = regridObj(ptopsigma)
                ptoptempi[t,:,:]            = regridObj(ptoptemp)
                ptopsalti[t,:,:]            = regridObj(ptopsalt)
                ptophvmi[t,:,:]             = regridObj(ptopsalt)
                ptopdepthi[t,:,:].mask      = msk.maski
                ptopsigmai[t,:,:].mask      = msk.maski
                ptoptempi[t,:,:].mask       = msk.maski
                ptopsalti[t,:,:].mask       = msk.maski
                ptophvmi[t,:,:].mask        = msk.maski
                ptopdepthia[t,:,:]          = ptopdepthi[t,:,:]*1.
                ptopdepthip[t,:,:]          = ptopdepthi[t,:,:]*1.
                ptopdepthii[t,:,:]          = ptopdepthi[t,:,:]*1.
                ptopdepthia[t,:,:].mask     = msk.maskAtl
                ptopdepthip[t,:,:].mask     = msk.maskPac
                ptopdepthii[t,:,:].mask     = msk.maskInd
                ptopsigmaia[t,:,:]          = ptopsigmai[t,:,:]*1.
                ptopsigmaip[t,:,:]          = ptopsigmai[t,:,:]*1.
                ptopsigmaii[t,:,:]          = ptopsigmai[t,:,:]*1.
                ptopsigmaia[t,:,:].mask     = msk.maskAtl
                ptopsigmaip[t,:,:].mask     = msk.maskPac
                ptopsigmaii[t,:,:].mask     = msk.maskInd
                ptoptempia[t,:,:]           = ptoptempi[t,:,:]*1.
                ptoptempip[t,:,:]           = ptoptempi[t,:,:]*1.
                ptoptempii[t,:,:]           = ptoptempi[t,:,:]*1.
                ptoptempia[t,:,:].mask      = msk.maskAtl
                ptoptempip[t,:,:].mask      = msk.maskPac
                ptoptempii[t,:,:].mask      = msk.maskInd
                ptopsaltia[t,:,:]           = ptopsalti[t,:,:]*1.
                ptopsaltip[t,:,:]           = ptopsalti[t,:,:]*1.
                ptopsaltii[t,:,:]           = ptopsalti[t,:,:]*1.
                ptopsaltia[t,:,:].mask      = msk.maskAtl
                ptopsaltip[t,:,:].mask      = msk.maskPac
                ptopsaltii[t,:,:].mask      = msk.maskInd
                ptophvmia[t,:,:]           = ptophvmi[t,:,:]*1.
                ptophvmip[t,:,:]           = ptophvmi[t,:,:]*1.
                ptophvmii[t,:,:]           = ptophvmi[t,:,:]*1.
                ptophvmia[t,:,:].mask      = msk.maskAtl
                ptophvmip[t,:,:].mask      = msk.maskPac
                ptophvmii[t,:,:].mask      = msk.maskInd

                ptopdepthi  = maskVal(ptopdepthi,  valmask)
                ptopdepthia = maskVal(ptopdepthia, valmask)
                ptopdepthip = maskVal(ptopdepthip, valmask)
                ptopdepthii = maskVal(ptopdepthii, valmask)
                ptopsigmai  = maskVal(ptopsigmai,  valmask)
                ptopsigmaia = maskVal(ptopsigmaia, valmask)
                ptopsigmaip = maskVal(ptopsigmaip, valmask)
                ptopsigmaii = maskVal(ptopsigmaii, valmask)
                ptoptempi   = maskVal(ptoptempi,   valmask)
                ptoptempia  = maskVal(ptoptempia,  valmask)
                ptoptempip  = maskVal(ptoptempip,  valmask)
                ptoptempii  = maskVal(ptoptempii,  valmask)
                ptopsalti   = maskVal(ptopsalti,   valmask)
                ptopsaltia  = maskVal(ptopsaltia,  valmask)
                ptopsaltip  = maskVal(ptopsaltip,  valmask)
                ptopsaltii  = maskVal(ptopsaltii,  valmask)
                ptophvmi   = maskVal(ptophvmi,   valmask)
                ptophvmia  = maskVal(ptophvmia,  valmask)
                ptophvmip  = maskVal(ptophvmip,  valmask)
                ptophvmii  = maskVal(ptophvmii,  valmask)

                # Free memory
                del(persbin,ptopdepth,ptopsigma,ptoptemp,ptopsalt,ptophvm) ; gc.collect()
                # Create cdms2 transient variables
                ptopdepthi  = cdm.createVariable(ptopdepthi,  axes = [timeyr, lati, loni], id = 'ptopdepthi')
                ptopsigmai  = cdm.createVariable(ptopsigmai,  axes = [timeyr, lati, loni], id = 'ptopsigmai')
                ptoptempi   = cdm.createVariable(ptoptempi,   axes = [timeyr, lati, loni], id = 'ptopthetaoi')
                ptopsalti   = cdm.createVariable(ptopsalti,   axes = [timeyr, lati, loni], id = 'ptopsoi')
                ptophvmi   = cdm.createVariable(ptophvmi,   axes = [timeyr, lati, loni], id = 'ptophvoi')
                ptopdepthia = cdm.createVariable(ptopdepthia, axes = [timeyr, lati, loni], id = 'ptopdepthia')
                ptopsigmaia = cdm.createVariable(ptopsigmaia, axes = [timeyr, lati, loni], id = 'ptopsigmaia')
                ptoptempia  = cdm.createVariable(ptoptempia,  axes = [timeyr, lati, loni], id = 'ptopthetaoia')
                ptopsaltia  = cdm.createVariable(ptopsaltia,  axes = [timeyr, lati, loni], id = 'ptopsoia')
                ptophvmia  = cdm.createVariable(ptophvmia,  axes = [timeyr, lati, loni], id = 'ptophvoia')
                ptopdepthip = cdm.createVariable(ptopdepthip, axes = [timeyr, lati, loni], id = 'ptopdepthip')
                ptopsigmaip = cdm.createVariable(ptopsigmaip, axes = [timeyr, lati, loni], id = 'ptopsigmaip')
                ptoptempip  = cdm.createVariable(ptoptempip,  axes = [timeyr, lati, loni], id = 'ptopthetaoip')
                ptopsaltip  = cdm.createVariable(ptopsaltip,  axes = [timeyr, lati, loni], id = 'ptopsoip')
                ptophvmip  = cdm.createVariable(ptophvmip,  axes = [timeyr, lati, loni], id = 'ptophvoip')
                ptopdepthii = cdm.createVariable(ptopdepthii, axes = [timeyr, lati, loni], id = 'ptopdepthii')
                ptopsigmaii = cdm.createVariable(ptopsigmaii, axes = [timeyr, lati, loni], id = 'ptopsigmaii')
                ptoptempii  = cdm.createVariable(ptoptempii,  axes = [timeyr, lati, loni], id = 'ptopthetaoii')
                ptopsaltii  = cdm.createVariable(ptopsaltii,  axes = [timeyr, lati, loni], id = 'ptopsoii')
                ptophvmii  = cdm.createVariable(ptophvmii,  axes = [timeyr, lati, loni], id = 'ptophvoii')
                
                # Compute zonal mean of bowl variables (1D)
                ptopdiz     = cdu.averager(ptopdepthi,  axis = 2)
                ptopdiza    = cdu.averager(ptopdepthia, axis = 2)
                ptopdizp    = cdu.averager(ptopdepthip, axis = 2)
                ptopdizi    = cdu.averager(ptopdepthii, axis = 2)
                ptopriz     = cdu.averager(ptopsigmai,  axis = 2)
                ptopriza    = cdu.averager(ptopsigmaia, axis = 2)
                ptoprizp    = cdu.averager(ptopsigmaip, axis = 2)
                ptoprizi    = cdu.averager(ptopsigmaii, axis = 2)
                ptoptiz     = cdu.averager(ptoptempi,   axis = 2)
                ptoptiza    = cdu.averager(ptoptempia,  axis = 2)
                ptoptizp    = cdu.averager(ptoptempip,  axis = 2)
                ptoptizi    = cdu.averager(ptoptempii,  axis = 2)
                ptopsiz     = cdu.averager(ptopsalti,   axis = 2)
                ptopsiza    = cdu.averager(ptopsaltia,  axis = 2)
                ptopsizp    = cdu.averager(ptopsaltip,  axis = 2)
                ptopsizi    = cdu.averager(ptopsaltii,  axis = 2)
                ptopviz     = cdu.averager(ptophvmi,   axis = 2)
                ptopviza    = cdu.averager(ptophvmia,  axis = 2)
                ptopvizp    = cdu.averager(ptophvmip,  axis = 2)
                ptopvizi    = cdu.averager(ptophvmii,  axis = 2)

                tpe6 = timc.clock()
                # Compute volume/temp/salinity of persistent ocean (global, per basin) (1D)
                persvp = persisti[t,:,:,:]*1. ; persvp.mask[...] = persisti.mask[t,:,:,:]
                persvp = npy.floor(persvp/98.)
                persvp = cdm.createVariable(persvp, axes = [rhoAxis, lati, loni], id = 'toto')
                persvp._FillValue = valmask ; persvp = maskVal(persvp, valmask)
                # volume (integral of depth * area)
                thickrij       = thickBini.data[t,...]*(1-thickBini.mask[t,...])
                temprij        = x1Bini.data[t,...]*(1-thickBini.mask[t,...])
                salrij         = x2Bini.data[t,...]*(1-thickBini.mask[t,...])
                hvmrij         = x3Bini.data[t,...]*(1-thickBini.mask[t,...])
                voltotij       = npy.ma.sum(thickrij, axis=0)
                voltot         = npy.ma.sum(voltotij*area.areai)
                volpersxy      = npy.ma.sum(persvp.data*thickrij, axis=0)
                volpersist [t] = npy.ma.sum(npy.ma.reshape(volpersxy*area.areai, (Nji*Nii)))
                volpersista[t] = npy.ma.sum(npy.ma.reshape(volpersxy*area.areaia,(Nji*Nii)))
                volpersistp[t] = npy.ma.sum(npy.ma.reshape(volpersxy*area.areaip,(Nji*Nii)))
                volpersisti[t] = npy.ma.sum(npy.ma.reshape(volpersxy*area.areaii,(Nji*Nii)))
                # Temp and salinity (average)
                    # add espilon to avoid diving by zero on land points
                tempersxy      = npy.ma.sum(persvp.data*temprij*thickrij, axis=0)/(volpersxy+0.0001)
                tempersist [t] = npy.ma.sum(npy.ma.reshape(tempersxy*area.areai, (Nji*Nii)))/area.areait
                tempersista[t] = npy.ma.sum(npy.ma.reshape(tempersxy*area.areaia,(Nji*Nii)))/area.areaita
                tempersistp[t] = npy.ma.sum(npy.ma.reshape(tempersxy*area.areaip,(Nji*Nii)))/area.areaitp
                tempersisti[t] = npy.ma.sum(npy.ma.reshape(tempersxy*area.areaii,(Nji*Nii)))/area.areaiti

                salpersxy      = npy.ma.sum(persvp.data*salrij*thickrij, axis=0)/(volpersxy+0.0001)
                salpersist [t] = npy.ma.sum(npy.ma.reshape(salpersxy*area.areai, (Nji*Nii)))/area.areait
                salpersista[t] = npy.ma.sum(npy.ma.reshape(salpersxy*area.areaia,(Nji*Nii)))/area.areaita
                salpersistp[t] = npy.ma.sum(npy.ma.reshape(salpersxy*area.areaip,(Nji*Nii)))/area.areaitp
                salpersisti[t] = npy.ma.sum(npy.ma.reshape(salpersxy*area.areaii,(Nji*Nii)))/area.areaiti

                hvmpersxy      = npy.ma.sum(persvp.data*hvmrij*thickrij, axis=0)/(volpersxy+0.0001)
                hvmpersist [t] = npy.ma.sum(npy.ma.reshape(hvmpersxy*area.areai, (Nji*Nii)))/area.areait
                hvmpersista[t] = npy.ma.sum(npy.ma.reshape(hvmpersxy*area.areaia,(Nji*Nii)))/area.areaita
                hvmpersistp[t] = npy.ma.sum(npy.ma.reshape(hvmpersxy*area.areaip,(Nji*Nii)))/area.areaitp
                hvmpersisti[t] = npy.ma.sum(npy.ma.reshape(hvmpersxy*area.areaii,(Nji*Nii)))/area.areaiti
                
                if debug:
                    print ' Integral persistent values:',voltot,volpersist[t],volpersista[t]
                    print '   %', volpersist[t]/voltot*100., volpersista[t]/volpersist[t]*100.
                    print '  area global/atl/pac/ind ', area.areait, area.areaita, area.areaitp, area.areaiti
                    print '   T , S glob ',tempersist[t] , salpersist[t]
                    print '   T , S atl  ',tempersista[t], salpersista[t]
                    print '   T , S pac  ',tempersistp[t], salpersistp[t]
                    print '   T , S ind  ',tempersisti[t], salpersisti[t]
                del(volpersxy,tempersxy,salpersxy)
                tpe7 = timc.clock()
                # CPU analysis
                if cpuan:
                    print ' Persistence CPU analysis t/nyrtc = ',t,nyrtc
                    print '    cpu1 = ', tpe1 - tpe0
                    print '    cpu2 = ', tpe2 - tpe1
                    print '    cpu3 = ', tpe3 - tpe2
                    print '    cpu4 = ', tpe4 - tpe3
                    print '    cpu5 = ', tpe5 - tpe4
                    print '    cpu6 = ', tpe6 - tpe5
                    print '    cpu7 = ', tpe7 - tpe6

            #
            # end of loop on t <==

            # Compute % of persistent ocean on the vertical
            persistm                = (cdu.averager(persistv, axis = 1)/cdu.averager(thickBini, axis = 1))
            persistm._FillValue     = valmask
            persistm                = mv.masked_where(persistm > valmask/10, persistm)
            persistm.mask           = msk.maski

            # Write % of persistent ocean, depth/temp/salinity of bowl 3D (time, lat, lon)
            persim = cdm.createVariable(persistm  , axes = [timeyr,lati,loni], id = 'persistmxy')
            ptopd  = cdm.createVariable(ptopdepthi, axes = [timeyr,lati,loni], id = 'ptopdepthxy')
            ptopt  = cdm.createVariable(ptoptempi , axes = [timeyr,lati,loni], id = 'ptopthetaoxy')
            ptops  = cdm.createVariable(ptopsalti , axes = [timeyr,lati,loni], id = 'ptopsoxy')
            ptophv  = cdm.createVariable(ptophvmi , axes = [timeyr,lati,loni], id = 'ptophvoxy')
            ptopsig  = cdm.createVariable(ptopsigmai , axes = [timeyr,lati,loni], id = 'ptopsigmaxy')

            # Write volume/temp/salinity of persistent ocean 1D (time)
            # Collapse onto basin axis
            if 'timeBasinAxesList' not in locals():
                timeBasinList = basinTimeList
                timeBasinList[0] = timeyr ; # Replace monthly with annual
                timeBasinAxesList = basinAxesList
                timeBasinAxesList[0] = timeyr ; # Replace monthly with annual
                timeBasinAxesList[2] = lati ; # Replace lat with regrid target
                timeBasinRhoAxesList = basinRhoAxesList
                timeBasinRhoAxesList[0] = timeyr ; # Replace monthly with annual
                timeBasinRhoAxesList[3] = lati ; # Replace lat with regrid target
            newshape    = list(persistiz.shape) ; newshape.insert(1,1)
            persistiz   = npy.ma.reshape(persistiz,newshape)
            persistiza  = npy.ma.reshape(persistiza,newshape)
            persistizp  = npy.ma.reshape(persistizp,newshape)
            persistizi  = npy.ma.reshape(persistizi,newshape)
            dbpz        = npy.ma.concatenate((persistiz,persistiza,persistizp,persistizi),axis=1)
            del(persistiz,persistiza,persistizp,persistizi) ; gc.collect()
            dbpz        = cdm.createVariable(dbpz,axes=timeBasinRhoAxesList,id='isonpers')

            newshape    = list(ptopdiz.shape) ; newshape.insert(1,1)
            ptopdiz     = npy.ma.reshape(ptopdiz,newshape)
            ptopdiza    = npy.ma.reshape(ptopdiza,newshape)
            ptopdizp    = npy.ma.reshape(ptopdizp,newshape)
            ptopdizi    = npy.ma.reshape(ptopdizi,newshape)
            dbpdz       = npy.ma.concatenate((ptopdiz,ptopdiza,ptopdizp,ptopdizi),axis=1)
            del(ptopdiz,ptopdiza,ptopdizp,ptopdizi) ; gc.collect()
            dbpdz       = cdm.createVariable(dbpdz,axes=timeBasinAxesList,id='ptopdepth')

            newshape    = list(ptopriz.shape) ; newshape.insert(1,1)
            ptopriz     = npy.ma.reshape(ptopriz,newshape)
            ptopriza    = npy.ma.reshape(ptopriza,newshape)
            ptoprizp    = npy.ma.reshape(ptoprizp,newshape)
            ptoprizi    = npy.ma.reshape(ptoprizi,newshape)
            dbprz       = npy.ma.concatenate((ptopriz,ptopriza,ptoprizp,ptoprizi),axis=1)
            del(ptopriz,ptopriza,ptoprizp,ptoprizi) ; gc.collect()
            dbprz       = cdm.createVariable(dbprz,axes=timeBasinAxesList,id='ptopsigma')

            newshape    = list(ptoptiz.shape) ; newshape.insert(1,1)
            ptoptiz     = npy.ma.reshape(ptoptiz,newshape)
            ptoptiza    = npy.ma.reshape(ptoptiza,newshape)
            ptoptizp    = npy.ma.reshape(ptoptizp,newshape)
            ptoptizi    = npy.ma.reshape(ptoptizi,newshape)

            dbptz       = npy.ma.concatenate((ptoptiz,ptoptiza,ptoptizp,ptoptizi),axis=1)
            del(ptoptiz,ptoptiza,ptoptizp,ptoptizi) ; gc.collect()
            dbptz       = cdm.createVariable(dbptz,axes=timeBasinAxesList,id='ptopthetao')

            newshape    = list(ptopsiz.shape) ; newshape.insert(1,1)
            ptopsiz     = npy.ma.reshape(ptopsiz,newshape)
            ptopsiza    = npy.ma.reshape(ptopsiza,newshape)
            ptopsizp    = npy.ma.reshape(ptopsizp,newshape)
            ptopsizi    = npy.ma.reshape(ptopsizi,newshape)
            dbpsz       = npy.ma.concatenate((ptopsiz,ptopsiza,ptopsizp,ptopsizi),axis=1)
            del(ptopsiz,ptopsiza,ptopsizp,ptopsizi) ; gc.collect()
            dbpsz       = cdm.createVariable(dbpsz,axes=timeBasinAxesList,id='ptopso')

            newshape    = list(volpersist.shape) ; newshape.insert(1,1)
            volperw     = npy.ma.reshape(volpersist*1.e-12 ,newshape)
            volperwa    = npy.ma.reshape(volpersista*1.e-12,newshape)
            volperwp    = npy.ma.reshape(volpersistp*1.e-12,newshape)
            volperwi    = npy.ma.reshape(volpersisti*1.e-12,newshape)
            volper      = npy.ma.concatenate((volperw,volperwa,volperwp,volperwi),axis=1)
            del(volperw,volperwa,volperwp,volperwi) ; gc.collect()
            volper      = cdm.createVariable(volper,axes=timeBasinList,id='volpers')
            temperw  = npy.ma.reshape(tempersist ,newshape)
            temperwa = npy.ma.reshape(tempersista,newshape)
            temperwp = npy.ma.reshape(tempersistp,newshape)
            temperwi = npy.ma.reshape(tempersisti,newshape)
            temper      = npy.ma.concatenate((temperw,temperwa,temperwp,temperwi),axis=1)
            del(temperw,temperwa,temperwp,temperwi) ; gc.collect()
            temper      = cdm.createVariable(temper,axes=timeBasinList,id='tempers')
            salperw  = npy.ma.reshape(salpersist ,newshape)
            salperwa = npy.ma.reshape(salpersista,newshape)
            salperwp = npy.ma.reshape(salpersistp,newshape)
            salperwi = npy.ma.reshape(salpersisti,newshape)
            salper      = npy.ma.concatenate((salperw,salperwa,salperwp,salperwi),axis=1)
            del(salperw,salperwa,salperwp,salperwi) ; gc.collect()
            salper      = cdm.createVariable(salper,axes=timeBasinList,id='salpers')
            hvmperw  = npy.ma.reshape(hvmpersist ,newshape)
            hvmperwa = npy.ma.reshape(hvmpersista,newshape)
            hvmperwp = npy.ma.reshape(hvmpersistp,newshape)
            hvmperwi = npy.ma.reshape(hvmpersisti,newshape)
            hvmper      = npy.ma.concatenate((hvmperw,hvmperwa,hvmperwp,hvmperwi),axis=1)
            del(hvmperw,hvmperwa,hvmperwp,hvmperwi) ; gc.collect()
            hvmper      = cdm.createVariable(hvmper,axes=timeBasinList,id='hvmpers')

            if tc == 0:
                # Global attributes
                dbpz.long_name      = 'Zonal persistence of isopycnal bins'
                dbpz.units          = '% of time'
                #
                persim.long_name    = 'Fraction of persistence on isopycnal bins'
                persim.units        = '% of column'
                ptopd.long_name     = 'Depth of shallowest persistent ocean on ison'
                ptopd.units         = 'm'
                ptopt.long_name     = 'Temp. of shallowest persistent ocean on ison'
                ptopt.units         = 'degrees_C'
                ptops.long_name     = 'Salinity of shallowest persistent ocean on ison'
                ptops.units         = soUnits
                ptopsig.long_name     = 'Density of shallowest persistent ocean on ison'
                ptopsig.units         = 'sigma_n'
                #
                dbpdz.long_name     = 'Zonal depth of shallowest persistent ocean on ison'
                dbpdz.units         = 'm'
                dbprz.long_name     = 'Zonal rhon of shallowest persistent ocean on ison'
                dbprz.units         = 'sigma_n'
                dbptz.long_name     = 'Zonal Temp. of shallowest persistent ocean on ison'
                dbptz.units         = 'degrees_C'
                dbpsz.long_name     = 'Zonal Salinity of shallowest persistent ocean on ison'
                dbpsz.units         = soUnits
                #
                volper.long_name    = 'Volume of persistent ocean'
                volper.units        = '1.e12 m^3'
                temper.long_name    = 'Temperature of persistent ocean'
                temper.units        = 'degrees_C'
                salper.long_name    = 'Salinity of persistent ocean'
                salper.units        = soUnits
                hvmper.long_name    = 'meridional transport of persistent ocean'
                hvmper.units        = 'm/s * m'
            # Write & append
            outFile_f.write(depthbini.astype('float32'), extend = 1, index = (trmin-tmin)/12) ; # Write out 4D variable first depth,rhon,lat,lon are written together
            outFile_f.write(thickbini.astype('float32'), extend = 1, index = (trmin-tmin)/12)
            outFile_f.write(x1bini.astype('float32')   , extend = 1, index = (trmin-tmin)/12)
            outFile_f.write(x2bini.astype('float32')   , extend = 1, index = (trmin-tmin)/12)
            outFile_f.write(x3bini.astype('float32')   , extend = 1, index = (trmin-tmin)/12)
            outFile_f.write(persim.astype('float32') , extend = 1, index = (trmin-tmin)/12) ; # Write out 3D variable first depth,lat,lon are written together
            outFile_f.write(ptopd.astype('float32')  , extend = 1, index = (trmin-tmin)/12)
            outFile_f.write(ptopt.astype('float32')  , extend = 1, index = (trmin-tmin)/12)
            outFile_f.write(ptops.astype('float32')  , extend = 1, index = (trmin-tmin)/12)
            outFile_f.write(ptopsig.astype('float32')  , extend = 1, index = (trmin-tmin)/12)
            outFile_f.write(dbpz.astype('float32')   , extend = 1, index = (trmin-tmin)/12)
            del(persim,ptopd,ptopt,ptops,dbpz) ; gc.collect()
            outFile_f.write(dbpdz.astype('float32')  , extend = 1, index = (trmin-tmin)/12)
            outFile_f.write(dbprz.astype('float32')  , extend = 1, index = (trmin-tmin)/12)
            outFile_f.write(dbptz.astype('float32')  , extend = 1, index = (trmin-tmin)/12)
            outFile_f.write(dbpsz.astype('float32')  , extend = 1, index = (trmin-tmin)/12)
            del(dbpdz,dbprz,dbptz,dbpsz) ; gc.collect()

            outFile_f.write(volper.astype('float32') , extend = 1, index = (trmin-tmin)/12)
            outFile_f.write(temper.astype('float32') , extend = 1, index = (trmin-tmin)/12)
            outFile_f.write(salper.astype('float32') , extend = 1, index = (trmin-tmin)/12)
            outFile_f.write(hvmper.astype('float32') , extend = 1, index = (trmin-tmin)/12)
            #
            tozp = timc.clock()
            #
            # Init zonal mean output variables
            # Collapse onto basin axis
            newshape    = list(depthBinz.shape) ; newshape.insert(1,1)
            depthBinz   = npy.ma.reshape(depthBinz,newshape)
            depthBinza  = npy.ma.reshape(depthBinza,newshape)
            depthBinzp  = npy.ma.reshape(depthBinzp,newshape)
            depthBinzi  = npy.ma.reshape(depthBinzi,newshape)
            dbz         = npy.ma.concatenate((depthBinz,depthBinza,depthBinzp,depthBinzi),axis=1)
            del(depthBinz,depthBinza,depthBinzp,depthBinzi) ; gc.collect()
            dbz         = cdm.createVariable(dbz,axes=timeBasinRhoAxesList,id='isondepth')
            thickBinz   = npy.ma.reshape(thickBinz,newshape)
            thickBinza  = npy.ma.reshape(thickBinza,newshape)
            thickBinzp  = npy.ma.reshape(thickBinzp,newshape)
            thickBinzi  = npy.ma.reshape(thickBinzi,newshape)
            tbz         = npy.ma.concatenate((thickBinz,thickBinza,thickBinzp,thickBinzi),axis=1)
            del(thickBinz,thickBinza,thickBinzp,thickBinzi) ; gc.collect()
            tbz         = cdm.createVariable(tbz,axes=timeBasinRhoAxesList,id='isonthick')
            volBinz     = npy.ma.reshape(volBinz,newshape)
            volBinza    = npy.ma.reshape(volBinza,newshape)
            volBinzp    = npy.ma.reshape(volBinzp,newshape)
            volBinzi    = npy.ma.reshape(volBinzi,newshape)
            vbz         = npy.ma.concatenate((volBinz,volBinza,volBinzp,volBinzi),axis=1)*1.e-12
            del(volBinz,volBinza,volBinzp,volBinzi) ; gc.collect()
            vbz         = cdm.createVariable(vbz,axes=timeBasinRhoAxesList,id='isonvol')
            x1Binz      = npy.ma.reshape(x1Binz,newshape)
            x1Binza     = npy.ma.reshape(x1Binza,newshape)
            x1Binzp     = npy.ma.reshape(x1Binzp,newshape)
            x1Binzi     = npy.ma.reshape(x1Binzi,newshape)
            x1bz        = npy.ma.concatenate((x1Binz,x1Binza,x1Binzp,x1Binzi),axis=1)
            del(x1Binz,x1Binza,x1Binzp,x1Binzi) ; gc.collect()
            x1bz        = cdm.createVariable(x1bz,axes=timeBasinRhoAxesList,id='isonthetao')
            x2Binz      = npy.ma.reshape(x2Binz,newshape)
            x2Binza     = npy.ma.reshape(x2Binza,newshape)
            x2Binzp     = npy.ma.reshape(x2Binzp,newshape)
            x2Binzi     = npy.ma.reshape(x2Binzi,newshape)
            x2bz        = npy.ma.concatenate((x2Binz,x2Binza,x2Binzp,x2Binzi),axis=1)
            del(x2Binz,x2Binza,x2Binzp,x2Binzi) ; gc.collect()
            x2bz        = cdm.createVariable(x2bz,axes=timeBasinRhoAxesList,id='isonso')

            x3Binz      = npy.ma.reshape(x3Binz,newshape)
            x3Binza     = npy.ma.reshape(x3Binza,newshape)
            x3Binzp     = npy.ma.reshape(x3Binzp,newshape)
            x3Binzi     = npy.ma.reshape(x3Binzi,newshape)
            x3bz        = npy.ma.concatenate((x3Binz,x3Binza,x3Binzp,x3Binzi),axis=1)
            del(x3Binz,x3Binza,x3Binzp,x3Binzi) ; gc.collect()
            x3bz        = cdm.createVariable(x3bz,axes=timeBasinRhoAxesList,id='isonso')
            
            if tc == 0:
                # Global attributes
                dbz.long_name   = 'Zonal depth of isopycnal'
                dbz.units       = 'm'
                tbz.long_name   = 'Zonal thickness of isopycnal'
                tbz.units       = 'm'
                vbz.long_name   = 'Volume of isopycnal'
                vbz.units       = '1.e12 m^3'
                x1bz.long_name  = thetaoLongName
                x1bz.units      = 'degrees_C'
                x2bz.long_name  = soLongName
                x2bz.units      = soUnits
                x3bz.long_name  = 'Merdional ocean transport'
                x3bz.units      = 'm/s * m'
                
                # Cleanup
            # Write & append
            outFile_f.write(dbz.astype('float32'),   extend = 1, index = (trmin-tmin)/12)
            outFile_f.write(tbz.astype('float32'),   extend = 1, index = (trmin-tmin)/12)
            outFile_f.write(vbz.astype('float32'),   extend = 1, index = (trmin-tmin)/12)
            outFile_f.write(x1bz.astype('float32'),  extend = 1, index = (trmin-tmin)/12)
            outFile_f.write(x2bz.astype('float32'),  extend = 1, index = (trmin-tmin)/12)
            outFile_f.write(x3bz.astype('float32'),  extend = 1, index = (trmin-tmin)/12)
            del(dbz,tbz,vbz,x1bz,x2bz,x3bz) ; gc.collect()

        # Write/append to file
        if mthout:
            outFileMon_f.write(depthBin.astype('float32'), extend = 1, index = trmin-tmin)
            outFileMon_f.write(thickBin.astype('float32'), extend = 1, index = trmin-tmin)
            outFileMon_f.write(x1Bin.astype('float32'),    extend = 1, index = trmin-tmin)
            outFileMon_f.write(x2Bin.astype('float32'),    extend = 1, index = trmin-tmin)
            outFileMon_f.write(x3Bin.astype('float32'),    extend = 1, index = trmin-tmin)
            del(depthBin,thickBin,x1Bin,x2Bin,x3Bin) ; gc.collect()
            outFileMon_f.sync()
        outFile_f.sync()
        tozf = timc.clock()
        print '   CPU of chunk inits         =', tucz0-tuc
        print '   CPU of density bining      =', ticz0-tucz0
        print '   CPU of masking and var def =', ticz-ticz0
        if tcdel >= 12:
            print '   CPU of annual mean compute =', toz-ticz
            print '   CPU of interpolation       =', tozi-toz
            print '   CPU of zonal mean          =', toziz-tozi
            print '   CPU of persistence compute =', tozp-toziz
        print '   CPU of chunk               =', tozf-tuc
        print '   Max memory use',resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1.e6,'GB'

    ######################################### END CHUNK #########################################
    # end loop on tc <===
    print '   CPU of inits       =', tin1-ti0
    print '     CPU inits detail =', tur-ti0, tmsk-tur, tarea-tmsk, tinit-tarea, tintrp-tinit, tin1-tintrp
    print ' [ Time stamp',(timc.strftime("%d/%m/%Y %H:%M:%S")),']'
    print ' Max memory use',resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1.e6,'GB'
    ratio =  12.*float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)/float(grdsize*tmax)
    print ' Ratio to grid*nyears',ratio,'kB/unit(size*nyears)'
    print ' CPU use, elapsed', timc.clock() - ti0, timeit.default_timer() - te0
    ratio = 1.e6*(timc.clock() - ti0)/float(grdsize*tmax)
    print ' Ratio to grid*nyears',ratio,'1.e-6 sec/unit(size*nyears)'

    ft.close()
    fs.close()
    # Global attributes
    globalAttWrite(outFile_f,options=None) ; # Use function to write standard global atts
    # Write binDensity version
    eosNeutralPath = str(eosNeutral.__code__).split(' ')[6]
    eosNeutralPath = replace(replace(eosNeutralPath,'"',''),',','') ; # Clean scraped path
    if(eosNeutralPath.find('/') < 0): # no path because in current
       eosNeutralPath = os.getcwd()+'/'
    eosNeutralPath='/home/nillod/Dev/Density_bin/Density_bining/binDensity.py' ## TODO: REMOVE
    outFile_f.binDensity_version = ' '.join(getGitInfo(eosNeutralPath)[0:3])
    outFile_f.close()
    if mthout:
        # Global attributes
        globalAttWrite(outFileMon_f,options=None) ; # Use function to write standard global atts
        # Write binDensity version
        outFileMon_f.binDensity_version = ' '.join(getGitInfo(eosNeutralPath)[0:3])
        outFileMon_f.close()
        print ' Wrote file: ',outFileMon
    if tcdel >= 12:
        print ' Wrote file: ',outFile

    # Cleanup variables
        del(timeBasinAxesList,timeBasinRhoAxesList) ; gc.collect()
    #for a in locals().iterkeys():
        #print a

if __name__ == "__main__":
    modelSo = '/data/nillod/Density_bining/so_Omon_IPSL-CM5A-LR_historical_r1i1p1_185001-185412.nc'
    modelThetao = '/data/nillod/Density_bining/thetao_Omon_IPSL-CM5A-LR_historical_r1i1p1_185001-185412.nc'
    modelVo = '/data/nillod/Density_bining/vo_Omon_IPSL-CM5A-LR_historical_r1i1p1_185001-185412.nc'
    modelAreacello = '/prodigfs/project/CMIP5/main/IPSL/IPSL-CM5A-LR/piControl/fx/ocean/fx/r0i0p0/latest/areacello/areacello_fx_IPSL-CM5A-LR_piControl_r0i0p0.nc'
    outfileDensity = '/data/nillod/Density_bining/out/cmip5.IPSL-VLR0.historical.rip.mon.ocean.Omon.density.MP2.nc'
    
    # Start Binning
    densityBin(modelThetao,modelSo,modelVo,modelAreacello,outfileDensity,timeint='1,12')

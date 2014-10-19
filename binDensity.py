#!/bin/env python
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
                    - TODO:
test

@author: durack1
"""

import ESMP,gc,os,resource,timeit ; #argparse,sys
import cdms2 as cdm
from cdms2 import CdmsRegrid
import cdutil as cdu
from durolib import fixVarUnits
import MV2 as mv
import numpy as npy
from string import replace
import time as timc

# Turn off numpy warnings
npy.seterr(all='ignore') ; # Cautious use of this turning all error reporting off - shouldn't be an issue as using masked arrays

# Function definitions

def maskVal(field,valmask):
    '''
    The maskVal() function applies a mask to an array provided
    
    Author:    Eric Guilyardi : Eric.Guilyardi@locean-ipsl.upmc.fr
    Co-author: Paul J. Durack : pauldurack@llnl.gov : @durack1.
    
    Created on Sun Sep 14 21:13:30 2014

    Inputs:
    ------
    - field     - 1D/2D/3D array
    - vamask    - 1D scalar of mask value
    
    Output:
    - field     - 1D/2D/3D masked array
    
    Usage:
    ------
    >>> from binDensity import maskVal
    >>> maskedVariable = maskVal(unMaskedVariable,valmask)

    Notes:
    -----
    - PJD 15 Sep 2014 - 
    '''
    field [npy.isnan(field.data)] = valmask
    field._FillValue = valmask
    field = mv.masked_where(field > valmask/10, field)
    return field


# Compute area of grid cells on earth
def computeArea(lon,lat):
    '''
    The computeArea() function calculates grid cell area assuming values are
    cell mid-points and formula area = R^2(lon2-lon1)*(sin(lat2) - sin(lat1))
    
    Author:    Eric Guilyardi : Eric.Guilyardi@locean-ipsl.upmc.fr
    Co-author: Paul J. Durack : pauldurack@llnl.gov : @durack1.
    
    Created on Sun Sep 14 21:13:30 2014

    Inputs:
    ------
    - lon   - 1D longitude  - >0, <360
    - lat   - 1D latitude   - >-90, <90
    
    Output:
    - area(lon,lat)     - 2D area array     - m^-^2
    
    Usage:
    ------
    >>> from binDensity import computeArea
    >>> computeArea(lon,lat)

    Notes:
    -----
    - PJD 15 Sep 2014 - 
    '''
    radius = 6371000. ; # Earth radius (metres)
    radconv = npy.pi/180.
    N_i = int(lon.shape[0])
    N_j = int(lat.shape[0])
    area = npy.ma.ones([N_j, N_i], dtype='float32')*0.
    lonr = lon[:] * radconv
    latr = lat[:] * radconv
    #loop
    for i in range(1,N_i-1):
        lonm1 = (lonr[i-1] + lonr[i]  )*0.5
        lonp1 = (lonr[i]   + lonr[i+1])*0.5
        for j in range(1,N_j-1):
            latm1 = (latr[j-1] + latr[j]  )*0.5
            latp1 = (latr[j]   + latr[j+1])*0.5
            area[j,i] = npy.float(radius**2 * (lonp1 - lonm1) * (npy.sin(latp1) - npy.sin(latm1)))
        # North and south bounds
        latm1 = ((-90.*radconv) + latr[0] )*0.5
        latp1 = (latr[0]        + latr[1] )*0.5
        area[0,i] = npy.float(radius**2 * (lonp1 - lonm1) * (npy.sin(latp1) - npy.sin(latm1)))
        latm1 = (latr[N_j-2] + latr[N_j-1])*0.5
        latp1 = (latr[N_j-1] + (90.*radconv)  )*0.5
        area[N_j-1,i] = npy.float(radius**2 * (lonp1 - lonm1) * (npy.sin(latp1) - npy.sin(latm1)))
    # East and west bounds
    area[:,0]     = area[:,1]
    area[:,N_i-1] = area[:,N_i-2] 

    return area


def eosNeutral(pottemp,salt):
    '''
    The eosNeutral() function takes potential temperature and salinity arguments
    and calculates approximate neutral density (gamma_a) which is returned as a
    variable. The function uses the McDougall & Jackett (2005) equation of state
    
    McDougall, T. J. and D. R. Jackett (2005) The material derivative of neutral
    density. Journal of Marine Research, 63 (1), pp 159-185. doi: 10.1357/0022240053693734    
    
    Author:    Eric Guilyardi : Eric.Guilyardi@locean-ipsl.upmc.fr
    Co-author: Paul J. Durack : pauldurack@llnl.gov : @durack1.
    
    Created on Sun Sep 14 21:13:30 2014

    Inputs:
    ------
    - pottemp(time,lev,lat,lon)     - 4D potential temperature  - deg_C
    - salt(time,lev,lat,lon)        - 4D salinity               - PSS-78
    
    Output:
    - rho(time,lev,lat,lon)         - 4D neutral density array  - kg m^-^3
    
    Usage:
    ------
    >>> from binDensity import eosNeutral
    >>> eosNeutral(pottemp,salt)
    >>> eosNeutral(20.,35.) ; # Check value 1024.5941675119673

    Notes:
    -----
    - PJD 14 Sep 2014 - 
    '''
    zt = pottemp
    zs = salt
    # neutral density
    zsr     = npy.ma.sqrt(zs)
    zr1     = ( ( -4.3159255086706703e-4*zt+8.1157118782170051e-2 )*zt+2.2280832068441331e-1 )*zt+1002.3063688892480
    zr2     = ( -1.7052298331414675e-7*zs-3.1710675488863952e-3*zt-1.0304537539692924e-4 )*zs
    zr3     = ( ( (-2.3850178558212048e-9*zt -1.6212552470310961e-7 )*zt+7.8717799560577725e-5 )*zt+4.3907692647825900e-5 )*zt + 1.0
    zr4     = ( ( -2.2744455733317707e-9*zt*zt+6.0399864718597388e-6)*zt-5.1268124398160734e-4 )*zs
    zr5     = ( -1.3409379420216683e-9*zt*zt-3.6138532339703262e-5)*zs*zsr
    zrho    = ( zr1 + zr2 ) / ( zr3 + zr4 + zr5 )
    return zrho 


def rhonGrid(rho_min,rho_int,rho_max,del_s1,del_s2):
    '''
    The rhonGrid() function computes grid for density variables 
    
    Author:    Eric Guilyardi : Eric.Guilyardi@locean-ipsl.upmc.fr
    Co-author: Paul J. Durack : pauldurack@llnl.gov : @durack1.
    
    Created on Sun Sep 14 21:13:30 2014

    Inputs:
    ------
    - rho_min    - scalar - minimum density                       (e.g. 18)
    - rho_int    - scalar - intermediate density                  (e.g. 26)
    - rho_max    - scalar - maximum density                       (2.g. 28)
    - del_s1     - scalar - delta_rho between rho_min and rho_int (e.g. 0.2)
    - del_s2     - scalar - delta_rho between rho_mintand rho_max (e.g. 0.1)
    
    Output:
    - s_s        - 1D array - Density grid
    - s_sax      - 1D array - Density grid for plotting axis (adds the last interval e.g. 28-28.1)
    - del_s      - 1D array - delta_rho
    - N_s        - integer  - dimension of density grid
    
    Usage:
    ------
    >>> from binDensity import rhonGrid
    >>> rhonGrid(rho_min,rho_int,rho_max,del_s1,del_s2)

    Notes:
    -----
    - PJD 14 Sep 2014 - 
    - EG  23 Sep 2014 - documentation 
    '''
    s_s1 = npy.arange(rho_min, rho_int, del_s1, dtype = npy.float32)
    s_s2 = npy.arange(rho_int, rho_max, del_s2, dtype = npy.float32)
    s_s  = npy.concatenate([s_s1, s_s2])
    N_s1 = len(s_s1)
    N_s2 = len(s_s2)
    N_s  = len(s_s)
    del_s = npy.concatenate([npy.tile(del_s1, N_s1), npy.tile(del_s2, N_s2)])
    s_sax = npy.append(s_s, s_s[N_s-1]+del_s2) # make axis
    return s_s, s_sax, del_s, N_s


def densityBin(fileT,fileS,fileFx,outFile,debug=True,timeint='all',mthout=False):
    '''
    The densityBin() function takes file and variable arguments and creates
    density persistence fields which are written to a specified outfile
    Author:    Eric Guilyardi : Eric.Guilyardi@locean-ipsl.upmc.fr
    Co-author: Paul J. Durack : pauldurack@llnl.gov : @durack1.
    
    Created on Sun Sep 14 21:13:30 2014

    Inputs:
    ------
    - fileT(time,lev,lat,lon)   - 4D potential temperature array
    - fileS(time,lev,lat,lon)   - 4D salinity array
    - fileFx(lat,lon)           - 2D array containing the cell area values
    - outFile(str)              - output file with full path specified.
    - debug <optional>          - boolean value
    - timeint <optional>        - specify temporal step for binning <init_idx>,<ncount>
    - mthout <optional>         - write annual data (False) or all monthly data (True)

    Usage:
    ------
    >>> from binDensity import densityBin
    >>> densityBin(fileT,fileS,fileFx,'./out.nc')

    Notes:
    -----
    - PJD 14 Sep 2014   - Initial function rewrite
    - PJD 18 Sep 2014   - Added gc.collect calls
    - PJD 18 Sep 2014   - Added so/thetao units check (fixVarUnits)
    - PJD 18 Sep 2014   - Renamed temp to thetao for clarity
    - EG  23 Sep 2014   - Clean up and documentation
    - PJD 23 Sep 2014   - Added so/thetao variable handles so variable attributes can be copied without a heavy memory overhead
    - EG  29 Sep 2014   - remove NaN in array inits
    '''

    # Keep track of time (CPU and elapsed)
    ti0 = timc.clock()
    te0 = timeit.default_timer()

    # netCDF compression
    comp = 1 ; # 0 for no compression
    cdm.setNetcdfShuffleFlag(comp)
    cdm.setNetcdfDeflateFlag(comp)
    cdm.setNetcdfDeflateLevelFlag(comp)
    cdm.setAutoBounds('on')
    
    # == Inits
    npy.set_printoptions(precision = 2)
    
    # Determine file name from inputs
    modeln = fileT.split('/')[-1].split('.')[1]

    # Declare and open files for writing too
    #/work/durack1/Shared/data_density/140915/cmip5.ACCESS1-0.historical.r1i1p1.mo.ocn.Omon.density.ver-1.nc
    outFile = replace(outFile,'.mo.','.an.')
    if os.path.isfile(outFile):
        os.remove(outFile)
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
        print 'Debug - File names:'
        print 'thetao:',fileT
        print 'so    :',fileS
        debugp = True
    else:
        debugp = False
    
    # Open files to read
    ft      = cdm.open(fileT)
    fs      = cdm.open(fileS)
    timeax  = ft.getAxis('time')
    
    # Need to add test to ensure thetao and so are equivalent sized (times equal)
    
    # Dates to read
    if timeint == 'all':
        tmin = 0
        tmax = timeax.shape[0]
    else:
        tmin = int(timeint.split(',')[0]) - 1
        tmax = tmin + int(timeint.split(',')[1])
    
    if debugp:
        print ' Debug mode'
    
    # Define temperature and salinity arrays
    thetao      = ft('thetao', time = slice(0,1))
    thetao_h    = ft['thetao'] ; # Create variable handle
    so          = fs('so', time = slice(0,1))
    so_h        = fs['so'] ; # Create variable handle
    
    # Read file attributes to carry on to output files
    list_file   = ft.attributes.keys()
    file_dic    = {}
    for i in range(0,len(list_file)):
        file_dic[i]=list_file[i],ft.attributes[list_file[i] ]
    
    # Read masking value
    valmask = so._FillValue
    
    # Read time and grid
    lon     = thetao_h.getLongitude()
    lat     = thetao_h.getLatitude()
    depth   = thetao_h.getLevel()
    bounds  = ft('lev_bnds')
    ingrid  = thetao_h.getGrid()
    
    # Read cell area
    ff      = cdm.open(fileFx)
    area    = ff('areacello')
    ff.close()
    
    # Define dimensions
    N_i = int(lon.shape[1])
    N_j = int(lon.shape[0])
    N_z = int(depth.shape[0])
    
    # Define sigma grid with zoom on higher densities
    rho_min = 19
    rho_int = 26
    rho_max = 28.5
    del_s1  = 0.2
    del_s2  = 0.1
    s_s, s_sax, del_s, N_s = rhonGrid(rho_min, rho_int, rho_max, del_s1, del_s2)
    s_s = npy.tile(s_s, N_i*N_j).reshape(N_i*N_j,N_s).transpose() # make 3D for matrix computation
    
    # ---------------------
    #  Init density bining
    # ---------------------
    # test point  
    itest = 80 
    jtest = 60
    ijtest = jtest*N_i + itest
    
    # Define time read interval (as function of 3D array size)
    # TODO: review to optimize
    grdsize = N_i * N_j * N_z

    # define number of months in each chunk
    if grdsize <= 1.e6:
        tcdel = min(120, tmax)
    elif grdsize <= 1.e7:
        tcdel = min(24, tmax)
    #tcdel = min(24, tmax) # faster than higher tcdel ?
    nyrtc = tcdel/12
    tcmax = (tmax-tmin)/tcdel ; # number of time chunks
    print ' ==> model:', modeln,' (grid size:', grdsize,')'
    print ' ==> time interval: ', tmin, tmax - 1
    print ' ==> size of time chunk, number of time chunks (memory optimization) :', tcdel, tcmax

    # inits
    # depth profiles:
    z_zt = depth[:]
    z_zw = bounds.data[:,0]
    max_depth_ocean = 6000. # maximum depth of ocean
    
    # File output inits
    s_axis              = cdm.createAxis(s_sax, id = 'rhon')
    s_axis.long_name    = 'Neutral density'
    s_axis.units        = 'kg m-3 (anomaly, minus 1000)'
    s_axis.designateLevel()

    # output arrays for each chunk
    tmp         = npy.ma.ones([tcdel, N_s+1, N_j*N_i], dtype='float32')*valmask
    depth_bin   = tmp.copy()
    depth_bin   = maskVal(depth_bin, valmask)
    thick_bin   = tmp.copy()
    thick_bin   = maskVal(thick_bin, valmask)
    x1_bin      = tmp.copy()
    x1_bin      = maskVal(x1_bin, valmask)
    x2_bin      = tmp.copy() ; del(tmp) ; gc.collect()
    x2_bin      = maskVal(x2_bin, valmask)
    
    # Target horizonal grid for interp 
    gridFile    = '140807_WOD13_masks.nc'
    gridFile_f  = cdm.open(gridFile)
    maskg       = gridFile_f('basinmask3')
    outgrid     = maskg.getGrid()
    maski       = maskg.mask ; # Global mask
    # Regional masks
    maskAtl = maski*1 ; maskAtl[...] = True
    idxa = npy.argwhere(maskg == 1).transpose()
    maskAtl[idxa[0],idxa[1]] = False
    maskPac = maski*1 ; maskPac[...] = True
    idxp = npy.argwhere(maskg == 2).transpose()
    maskPac[idxp[0],idxp[1]] = False
    maskInd = maski*1 ; maskInd[...] = True
    idxi = npy.argwhere(maskg == 3).transpose()
    maskInd[idxi[0],idxi[1]] = False
    
    loni    = maskg.getLongitude()
    lati    = maskg.getLatitude()
    Nii     = int(loni.shape[0])
    Nji     = int(lati.shape[0])
    # Compute area of target grid and zonal and global sums
    areai = computeArea(loni[:], lati[:])
    gridFile_f.close()
    areai.mask = maski
    areaia = areai*1.
    areaia.mask = maskAtl
    areaip = areai*1.
    areaip.mask = maskPac
    areaii = areai*1.
    areaii.mask = maskInd
    areazt  = cdu.averager(areai , axis=1, action='sum')
    areazta = cdu.averager(areaia, axis=1, action='sum')
    areaztp = cdu.averager(areaip, axis=1, action='sum')
    areazti = cdu.averager(areaii, axis=1, action='sum')
    areait  = cdu.averager(npy.reshape(areai ,(Nji*Nii)), action='sum')
    areaita = cdu.averager(npy.reshape(areaia,(Nji*Nii)), action='sum')
    areaitp = cdu.averager(npy.reshape(areaip,(Nji*Nii)), action='sum')
    areaiti = cdu.averager(npy.reshape(areaii,(Nji*Nii)), action='sum')
    
    # Interpolation init (regrid)
    ESMP.ESMP_Initialize()
    regridObj = CdmsRegrid(ingrid,outgrid,depth_bin.dtype,missing=valmask,regridMethod='linear',regridTool='esmf')
    
    # Global arrays on target grid init
    depthBini = npy.ma.ones([nyrtc, N_s+1, Nji, Nii], dtype='float32')*valmask 
    thickBini,x1Bini,x2Bini = [npy.ma.ones(npy.shape(depthBini)) for _ in range(3)]
    # Basin
    depthBinia,thickBinia,x1Binia,x2Binia,depthBinip,thickBinip,\
    x1Binip,x2Binip,depthBinii,thickBinii,x1Binii,x2Binii = [npy.ma.ones(npy.shape(depthBini)) for _ in range(12)]
    # Persistence arrays on original grid
    persist     = npy.ma.ones([nyrtc, N_s+1, N_j, N_i], dtype='float32')*valmask
    persisti,persistia,persistip,persistii,persistv = [npy.ma.ones(npy.shape(depthBini)) for _ in range(5)]
    # Persistence arrays on target grid
    persistm   = npy.ma.ones([nyrtc, Nji, Nii], dtype='float32')*valmask
    ptopdepthi,ptopsigmai,ptoptempi,ptopsalti = [npy.ma.ones(npy.shape(persistm)) for _ in range(4)]
    # Basin
    ptopdepthia,ptopsigmaia,ptoptempia,ptopsaltia,ptopdepthip,ptopsigmaip,\
    ptoptempip,ptopsaltip,ptopdepthii,ptopsigmaii,ptoptempii,ptopsaltii = [npy.ma.ones(npy.shape(persistm)) for _ in range(12)]
    # Volume/thetao/so of persistent ocean
    
    volpersist = npy.ma.ones([nyrtc], dtype='float32')*valmask
    volpersista,volpersistp,volpersisti,tempersist,tempersista,tempersistp,tempersisti,salpersist,salpersista,salpersistp,salpersisti = [npy.ma.ones(npy.shape(volpersist)) for _ in range(11)]
    #
    # -----------------------------------------
    #  Density bining loop (on time chunks tc)
    # -----------------------------------------
    for tc in range(tcmax):
        tuc     = timc.clock()
        # read tcdel month by tcdel month to optimise memory
        trmin   = tmin + tc*tcdel ; # define as function of tc and tcdel
        trmax   = tmin + (tc+1)*tcdel ; # define as function of tc and tcdel
        print ' --> time chunk (bounds) = ',tc, '/',tcmax,' (',trmin,trmax-1,')', modeln
        thetao  = ft('thetao', time = slice(trmin,trmax))
        so      = fs('so'    , time = slice(trmin,trmax))
        time    = thetao.getTime()
        
        # Test variable units
        [so,soFixed] = fixVarUnits(so,'so',True)#,'logfile.txt')
        if soFixed:
            print '     so:     units corrected'
        [thetao,thetaoFixed] = fixVarUnits(thetao,'thetao',True)#,'logfile.txt')
        if thetaoFixed:
            print '     thetao: units corrected'        
        
        # Compute neutral density 
        rhon = eosNeutral(thetao,so)-1000.
        
        # reorganise i,j dims in single dimension data (speeds up loops)
        thetao  = npy.reshape(thetao, (tcdel, N_z, N_i*N_j))
        so      = npy.reshape(so  , (tcdel, N_z, N_i*N_j))
        rhon    = npy.reshape(rhon, (tcdel, N_z, N_i*N_j))
        
        # init output arrays for bined fields
        depth_bin[...] = valmask
        thick_bin[...] = valmask
        x1_bin[...]    = valmask
        x2_bin[...]    = valmask
        
        # Loop on time within chunk tc
        for t in range(trmax-trmin): 
            # x1 contents on vertical (not yet implemented - may be done to ensure conservation)
            x1_content  = thetao.data[t] 
            x2_content  = so.data[t] 
            vmask_3D    = mv.masked_values(thetao.data[t], valmask).mask 
            # find non-masked points
            nomask      = npy.equal(vmask_3D[0],0)
            # init arrays for this time chunk
            z_s         = npy.ones((N_s+1, N_i*N_j))*valmask
            c1_s        = npy.ones((N_s+1, N_i*N_j))*valmask
            c2_s        = npy.ones((N_s+1, N_i*N_j))*valmask
            t_s         = npy.ones((N_s+1, N_i*N_j))*valmask
            szmin       = npy.ones((N_i*N_j))*valmask
            szmax       = npy.ones((N_i*N_j))*valmask
            i_min       = npy.ones((N_i*N_j))*0
            i_max       = npy.ones((N_i*N_j))*0
            delta_rho   = npy.ones((N_i*N_j))*valmask
            # find bottom level at each lat/lon point
            i_bottom            = vmask_3D.argmax(axis=0)-1
            z_s [N_s, nomask]   = z_zw[i_bottom[nomask]+1] ; # Cell depth limit
            c1_s[N_s, nomask]   = x1_content[N_z-1,nomask] ; # Cell bottom temperature/salinity
            c2_s[N_s, nomask]   = x2_content[N_z-1,nomask] ; # Cell bottom tempi_profilerature/salinity
            # init arrays as a function of depth = f(z)
            s_z = rhon.data[t]
            c1_z = x1_content
            c2_z = x2_content
            # Extract a strictly increasing sub-profile
            i_min[nomask]           = s_z.argmin(axis=0)[nomask]
            i_max[nomask]           = s_z.argmax(axis=0)[nomask]-1
            i_min[i_min > i_max]    = i_max[i_min > i_max]
            # Test on bottom minus surface stratification to check that it is larger than delta_rho
            delta_rho[nomask]           = s_z[i_bottom[nomask],nomask] - s_z[0,nomask]
            i_min[delta_rho < del_s1]   = 0
            i_max[delta_rho < del_s1]   = i_bottom[delta_rho < del_s1]
            
            # General case
            # find min/max of density for each z profile
            for i in range(N_i*N_j):
                if nomask[i]:
                    szmin[i] = s_z[i_min[i],i]
                    szmax[i] = s_z[i_max[i],i]
                else:
                    szmin[i] = 0.
                    szmax[i] = rho_max+10.            
            # Find indices between density min and density max
            #
            # Construct arrays of szm/c1m/c2m = s_z[i_min[i]:i_max[i],i] and valmask otherwise
            # same for zzm from z_zt 
            szm = s_z*1. ;  szm[...] = valmask
            zzm = s_z*1. ;  zzm[...] = valmask
            c1m = c1_z*1. ; c1m[...] = valmask
            c2m = c2_z*1. ; c2m[...] = valmask
            
            for k in range(N_z):
                k_ind = i_min*1.; k_ind[:] = valmask
                k_ind = npy.argwhere( (k >= i_min) & (k <= i_max))
                szm[k,k_ind] = s_z [k,k_ind]
                c1m[k,k_ind] = c1_z[k,k_ind]
                c2m[k,k_ind] = c2_z[k,k_ind]
                zzm[k,k_ind] = z_zt[k]
                
            # interpolate depth(z) (=z_zt) to depth(s) at s_s densities (=z_s) using density(z) (=s_z)
            # TODO: no loop 
            for i in range(N_i*N_j):
                if nomask[i]:
                    z_s [0:N_s,i] = npy.interp(s_s[:,i], szm[:,i], zzm[:,i], right = valmask) ; # consider spline           
                    c1_s[0:N_s,i] = npy.interp(z_s[0:N_s,i], zzm[:,i], c1m[:,i], right = valmask) 
                    c2_s[0:N_s,i] = npy.interp(z_s[0:N_s,i], zzm[:,i], c2m[:,i], right = valmask) 
            # if level in s_s has lower density than surface, isopycnal is put at surface (z_s = 0)
            inds = npy.argwhere(s_s < szmin).transpose()
            z_s [inds[0],inds[1]] = 0.
            c1_s[inds[0],inds[1]] = valmask
            c2_s[inds[0],inds[1]] = valmask
            # if level of s_s has higher density than bottom density, 
            # isopycnal is set to bottom (z_s = z_zw[i_bottom])
            inds = npy.argwhere(s_s > szmax).transpose()
            z_s [inds[0],inds[1]] = z_s[N_s,inds[1]]
            c1_s[inds[0],inds[1]] = valmask
            c2_s[inds[0],inds[1]] = valmask
            idzmc1 = npy.argwhere(c1_s == valmask).transpose()
            z_s [idzmc1[0],idzmc1[1]] = valmask
            # Thickness of isopycnal
            t_s [0,:] = z_s [0,:] 
            t_s [1:N_s,:] = z_s[1:N_s,:]-z_s[0:N_s-1,:]
            inds = npy.argwhere( (t_s <= 0.) ^ (t_s >= max_depth_ocean)).transpose()
            t_s [inds[0],inds[1]] = valmask
            t_s [idzmc1[0],idzmc1[1]] = valmask  
            if debugp and t < 0:
                i = ijtest
                print ' s_s[i]', s_s[:,i]
                print ' szm[i]', szm[:,i]
                print ' zzm[i]', zzm[:,i]
                print ' z_s[i]', z_s[:,i]
                print ' t_s[i]', t_s[:,i]
             # assign to final arrays
            depth_bin[t,:,:] = z_s
            thick_bin[t,:,:] = t_s
            x1_bin[t,:,:]    = c1_s
            x2_bin[t,:,:]    = c2_s
        #
        # end of loop on t <===      
        #      
        # Free memory 
        del(rhon, x1_content, x2_content, vmask_3D, szm, zzm, c1m, c2m, z_s, c1_s, c2_s, t_s, inds, c1_z, c2_z) ; gc.collect()
        # Reshape i*j back to i,j
        depth_bino = npy.reshape(depth_bin, (tcdel, N_s+1, N_j, N_i))
        thick_bino = npy.reshape(thick_bin, (tcdel, N_s+1, N_j, N_i))
        x1_bino    = npy.reshape(x1_bin,    (tcdel, N_s+1, N_j, N_i))
        x2_bino    = npy.reshape(x2_bin,    (tcdel, N_s+1, N_j, N_i))
        
        # Wash mask (from temp) over variables
        maskb           = mv.masked_values(x1_bino, valmask).mask
        depth_bino.mask = maskb
        x1_bino.mask    = maskb
        x2_bino.mask    = maskb
        depth_bino      = maskVal(depth_bino, valmask)
        x1_bino         = maskVal(x1_bino,    valmask)
        x2_bino         = maskVal(x2_bino,    valmask)

        maskt           = mv.masked_values(thick_bino, valmask).mask
        thick_bino.mask = maskt
        thick_bino      = maskVal(thick_bino, valmask)
        #
        del (maskb, maskt)
        if debugp and (tc < 0):
            # test write
            i = itest
            j = jtest
            #print 'ind = ',ind
            print 'test point',i,j, area[j,i]
            print 'lon,lat',lon[j,i],lat[j,i]
            print 'depth_bin', depth_bino[0,:,j,i]
            print 'thick_bin', thick_bino[0,:,j,i]
            print 'x1_bin', x1_bino[0,:,j,i]
            print 'x2_bin', x2_bino[0,:,j,i]
        #
        # Output files as netCDF
        # Def variables 
        depthBin = cdm.createVariable(depth_bino, axes = [time, s_axis, ingrid], id = 'isondepth')
        thickBin = cdm.createVariable(thick_bino, axes = [time, s_axis, ingrid], id = 'isonthick')
        x1Bin    = cdm.createVariable(x1_bino   , axes = [time, s_axis, ingrid], id = 'thetao')
        x2Bin    = cdm.createVariable(x2_bino   , axes = [time, s_axis, ingrid], id = 'so')
        if mthout:
            if tc == 0:
                depthBin.long_name  = 'Depth of isopycnal'
                depthBin.units      = 'm'
                thickBin.long_name  = 'Thickness of isopycnal'
                thickBin.units      = 'm'
                x1Bin.long_name     = thetao_h.long_name
                x1Bin.units         = 'C'
                x2Bin.long_name     = so_h.long_name
                x2Bin.units         = so_h.units
                outFileMon_f.write(area) ; # Added area so isonvol can be computed
                # write global attributes (inherited from thetao file)
                for i in range(0,len(file_dic)):
                    dm = file_dic[i]
                    setattr(outFileMon_f,dm[0],dm[1])
                    post_txt = 'Density bining via densit_bin.py using delta_sigma = '+str(del_s1)+' and '+str(del_s2)
                    setattr(outFileMon_f , 'Post_processing_history', post_txt)
                    setattr(outFile_f, 'Post_processing_history', post_txt)
        
        # -------------------------------------------------------------
        #  Compute annual mean, persistence, make zonal mean and write
        # -------------------------------------------------------------
        ticz = timc.clock()
        if tcdel >= 12:
            # Annual mean
            dy  = cdu.averager(npy.reshape(depthBin, (nyrtc, 12, N_s+1, N_j, N_i)), axis=1)
            ty  = cdu.averager(npy.reshape(thickBin, (nyrtc, 12, N_s+1, N_j, N_i)), axis=1)
            x1y = cdu.averager(npy.reshape(x1Bin,    (nyrtc, 12, N_s+1, N_j, N_i)), axis=1)
            x2y = cdu.averager(npy.reshape(x2Bin,    (nyrtc, 12, N_s+1, N_j, N_i)), axis=1)
            # create annual time axis
            timeyr          = cdm.createAxis(dy.getAxis(0))
            timeyr.id       = 'time'
            timeyr.units    = time.units
            timeyr.designateTime()
    
            dy   = cdm.createVariable(dy  , axes = [timeyr, s_axis, ingrid], id = 'isondy')
            ty   = cdm.createVariable(ty  , axes = [timeyr, s_axis, ingrid], id = 'isonty')
            x1y  = cdm.createVariable(x1y , axes = [timeyr, s_axis, ingrid], id = 'isonx1y')
            x2y  = cdm.createVariable(x2y , axes = [timeyr, s_axis, ingrid], id = 'isonx2y')
            
            toz = timc.clock()
                
            # Interpolate onto common grid
            for t in range(nyrtc):
                for ks in range(N_s+1):
                    # Global
                    depthBini[t,ks,:,:]         = regridObj(dy [t,ks,:,:])
                    thickBini[t,ks,:,:]         = regridObj(ty [t,ks,:,:])
                    x1Bini[t,ks,:,:]            = regridObj(x1y[t,ks,:,:])
                    x2Bini[t,ks,:,:]            = regridObj(x2y[t,ks,:,:])
                    depthBini[t,ks,:,:].mask    = maski
                    thickBini[t,ks,:,:].mask    = maski
                    x1Bini[t,ks,:,:].mask       = maski
                    x2Bini[t,ks,:,:].mask       = maski
                    # Atl
                    depthBinia[t,ks,:,:]        = depthBini[t,ks,:,:]*1.
                    thickBinia[t,ks,:,:]        = thickBini[t,ks,:,:]*1.
                    x1Binia[t,ks,:,:]           = x1Bini[t,ks,:,:]*1.
                    x2Binia[t,ks,:,:]           = x2Bini[t,ks,:,:]*1.
                    depthBinia[t,ks,:,:].mask   = maskAtl
                    thickBinia[t,ks,:,:].mask   = maskAtl
                    x1Binia[t,ks,:,:].mask      = maskAtl
                    x2Binia[t,ks,:,:].mask      = maskAtl
                    # Pac
                    depthBinip[t,ks,:,:]        = depthBini[t,ks,:,:]*1.
                    thickBinip[t,ks,:,:]        = thickBini[t,ks,:,:]*1.
                    x1Binip[t,ks,:,:]           = x1Bini[t,ks,:,:]*1.
                    x2Binip[t,ks,:,:]           = x2Bini[t,ks,:,:]*1.
                    depthBinip[t,ks,:,:].mask   = maskPac
                    thickBinip[t,ks,:,:].mask   = maskPac
                    x1Binip[t,ks,:,:].mask      = maskPac
                    x2Binip[t,ks,:,:].mask      = maskPac
                    # Ind
                    depthBinii[t,ks,:,:]        = depthBini[t,ks,:,:]*1.
                    thickBinii[t,ks,:,:]        = thickBini[t,ks,:,:]*1.
                    x1Binii[t,ks,:,:]           = x1Bini[t,ks,:,:]*1.
                    x2Binii[t,ks,:,:]           = x2Bini[t,ks,:,:]*1.
                    depthBinii[t,ks,:,:].mask   = maskInd
                    thickBinii[t,ks,:,:].mask   = maskInd
                    x1Binii[t,ks,:,:].mask      = maskInd
                    x2Binii[t,ks,:,:].mask      = maskInd
            # Free memory
            del(dy, ty, x1y, x2y); gc.collect()
            # Global
            depthBini   = maskVal(depthBini, valmask)
            thickBini   = maskVal(thickBini, valmask)
            x1Bini      = maskVal(x1Bini, valmask)
            x2Bini      = maskVal(x2Bini, valmask)
            # Atl
            depthBinia  = maskVal(depthBinia, valmask)
            thickBinia  = maskVal(thickBinia, valmask)
            x1Binia     = maskVal(x1Binia, valmask)
            x2Binia     = maskVal(x2Binia, valmask)
            # Pac
            depthBinip  = maskVal(depthBinip, valmask)
            thickBinip  = maskVal(thickBinip, valmask)
            x1Binip     = maskVal(x1Binip, valmask)
            x2Binip     = maskVal(x2Binip, valmask)
            # Ind
            depthBinii  = maskVal(depthBinii, valmask)
            thickBinii  = maskVal(thickBinii, valmask)
            x1Binii     = maskVal(x1Binii, valmask)
            x2Binii     = maskVal(x2Binii, valmask)

            depthbini  = cdm.createVariable(depthBini,  axes = [timeyr, s_axis, lati, loni], id = 'isondepthg')
            thickbini  = cdm.createVariable(thickBini,  axes = [timeyr, s_axis, lati, loni], id = 'isonthickg')
            x1bini     = cdm.createVariable(x1Bini   ,  axes = [timeyr, s_axis, lati, loni], id = 'thetaog')
            x2bini     = cdm.createVariable(x2Bini   ,  axes = [timeyr, s_axis, lati, loni], id = 'sog')
            if tc == 0:
                depthbini.long_name  = 'Depth of isopycnal'
                depthbini.units      = 'm'
                thickbini.long_name  = 'Thickness of isopycnal'
                thickbini.units      = 'm'
                x1bini.long_name     = thetao_h.long_name
                x1bini.units         = 'C'
                x2bini.long_name     = so_h.long_name
                x2bini.units         = so_h.units
     
            tozi = timc.clock()

            # Compute zonal mean
            # Global
            depthBinz   = cdu.averager(depthBini,   axis = 3)
            thickBinz   = cdu.averager(thickBini,   axis = 3)
            x1Binz      = cdu.averager(x1Bini,      axis = 3)
            x2Binz      = cdu.averager(x2Bini,      axis = 3)
            # Atl
            depthBinza  = cdu.averager(depthBinia,  axis = 3)
            thickBinza  = cdu.averager(thickBinia,  axis = 3)
            x1Binza     = cdu.averager(x1Binia,     axis = 3)
            x2Binza     = cdu.averager(x2Binia,     axis = 3)
            # Pac
            depthBinzp  = cdu.averager(depthBinip,  axis = 3)
            thickBinzp  = cdu.averager(thickBinip,  axis = 3)
            x1Binzp     = cdu.averager(x1Binip,     axis = 3)
            x2Binzp     = cdu.averager(x2Binip,     axis = 3)
            # Ind
            depthBinzi  = cdu.averager(depthBinii,  axis = 3)
            thickBinzi  = cdu.averager(thickBinii,  axis = 3)
            x1Binzi     = cdu.averager(x1Binii,     axis = 3)
            x2Binzi     = cdu.averager(x2Binii,     axis = 3)
            # Compute volume of isopycnals
            volBinz     = thickBinz  * areazt
            volBinza    = thickBinza * areazta
            volBinzp    = thickBinzp * areaztp
            volBinzi    = thickBinzi * areazti

            toziz = timc.clock()

            # Compute annual persistence of isopycnal bins (from their thickness): 'persist' array
            #  = percentage of time bin is occupied during each year (annual bowl if % < 100)
            for t in range(nyrtc):
                idxvm = npy.ma.ones([12, N_s+1, N_j, N_i], dtype='float32')*valmask 
                inim = t*12
                finm = t*12 + 12
                idxvm = 1-mv.masked_values(thick_bino[inim:finm,:,:,:], valmask).mask 
                persist[t,:,:,:] = cdu.averager(idxvm, axis = 0) * 100.
                # Shallowest persistent ocean index: p_top (2D)
                maskp = persist[t,:,:,:]*1. ; maskp[...] = valmask
                maskp = mv.masked_values(persist[t,:,:,:] >= 99., 1.).mask
                maskp = npy.reshape(maskp, (N_s+1, N_j*N_i))
                p_top = maskp.argmax(axis=0) 
                del(maskp) ; gc.collect()
                # Define properties on bowl (= shallowest persistent ocean)
                ptopdepth = npy.ma.ones([N_j*N_i], dtype='float32')*valmask 
                ptopsigma,ptoptemp,ptopsalt = [npy.ma.ones(npy.shape(ptopdepth)) for _ in range(3)]
                # TODO: can we remove the loop ?
                for i in range(N_j*N_i): 
                    ptopdepth[i]    = depth_bin[t,p_top[i],i]
                    ptoptemp[i]     = x1_bin[t,p_top[i],i]
                    ptopsalt[i]     = x2_bin[t,p_top[i],i]
                ptopsigma = ptopdepth*0. + s_sax [p_top] # to keep mask of ptopdepth
                ptopdepth = npy.reshape(ptopdepth, (N_j, N_i))
                ptopsigma = npy.reshape(ptopsigma, (N_j, N_i))
                ptoptemp  = npy.reshape(ptoptemp , (N_j, N_i))
                ptopsalt  = npy.reshape(ptopsalt , (N_j, N_i))
                # Create variables to attribute right axis for zonal mean
                ptopdepth = cdm.createVariable(ptopdepth, axes = [ingrid], id = 'ptopdepth')
                ptopsigma = cdm.createVariable(ptopsigma, axes = [ingrid], id = 'ptopsigma')
                ptoptemp  = cdm.createVariable(ptoptemp , axes = [ingrid], id = 'ptoptemp')
                ptopsalt  = cdm.createVariable(ptopsalt , axes = [ingrid], id = 'ptopsalt')
                # Mask persist where value is zero
                persist._FillValue = valmask
                persist = mv.masked_where(persist <= 1.e-6, persist)
                persbin = cdm.createVariable(persist, axes = [timeyr, s_axis, ingrid], id = 'isonpers') 
                # Interpolate to target grid and create basin variables 
                # TODO: can we remove the loop ?
                for ks in range(N_s+1):
                    persisti [t,ks,:,:]         = regridObj(persbin[t,ks,:,:])
                    persisti [t,ks,:,:].mask    = maski
                    persistia[t,ks,:,:]         = persisti[t,ks,:,:]*1.
                    persistia[t,ks,:,:].mask    = maskAtl
                    persistip[t,ks,:,:]         = persisti[t,ks,:,:]*1.
                    persistip[t,ks,:,:].mask    = maskPac
                    persistii[t,ks,:,:]         = persisti[t,ks,:,:]*1.
                    persistii[t,ks,:,:].mask    = maskInd

                persisti    = maskVal(persisti,  valmask)
                persistia   = maskVal(persistia, valmask)
                persistip   = maskVal(persistip, valmask)
                persistii   = maskVal(persistii, valmask)
                # Compute zonal mean (2D)
                persistiz   = cdu.averager(persisti,  axis = 3)
                persistiza  = cdu.averager(persistia, axis = 3)
                persistizp  = cdu.averager(persistip, axis = 3)
                persistizi  = cdu.averager(persistii, axis = 3)
                # Persistence * thickness (used to compute % of column that is persistent - see below)
                persistv[t,:,:,:]           = persisti[t,:,:,:] * thickBini[t,:,:,:]
                persistv                    = maskVal(persistv, valmask)
                # Depth, temperature and salinity on bowl (i.e. at shallowest persistent ocean) (2D) 
                ptopdepthi[t,:,:]           = regridObj(ptopdepth)
                ptopsigmai[t,:,:]           = regridObj(ptopsigma)
                ptoptempi[t,:,:]            = regridObj(ptoptemp)
                ptopsalti[t,:,:]            = regridObj(ptopsalt)
                ptopdepthi[t,:,:].mask      = maski
                ptopsigmai[t,:,:].mask      = maski
                ptoptempi[t,:,:].mask       = maski
                ptopsalti[t,:,:].mask       = maski
                ptopdepthia[t,:,:]          = ptopdepthi[t,:,:]*1.
                ptopdepthip[t,:,:]          = ptopdepthi[t,:,:]*1.
                ptopdepthii[t,:,:]          = ptopdepthi[t,:,:]*1.
                ptopdepthia[t,:,:].mask     = maskAtl
                ptopdepthip[t,:,:].mask     = maskPac
                ptopdepthii[t,:,:].mask     = maskInd
                ptopsigmaia[t,:,:]          = ptopsigmai[t,:,:]*1.
                ptopsigmaip[t,:,:]          = ptopsigmai[t,:,:]*1.
                ptopsigmaii[t,:,:]          = ptopsigmai[t,:,:]*1.
                ptopsigmaia[t,:,:].mask     = maskAtl
                ptopsigmaip[t,:,:].mask     = maskPac
                ptopsigmaii[t,:,:].mask     = maskInd
                ptoptempia[t,:,:]           = ptoptempi[t,:,:]*1.
                ptoptempip[t,:,:]           = ptoptempi[t,:,:]*1.
                ptoptempii[t,:,:]           = ptoptempi[t,:,:]*1.
                ptoptempia[t,:,:].mask      = maskAtl
                ptoptempip[t,:,:].mask      = maskPac
                ptoptempii[t,:,:].mask      = maskInd
                ptopsaltia[t,:,:]           = ptopsalti[t,:,:]*1.
                ptopsaltip[t,:,:]           = ptopsalti[t,:,:]*1.
                ptopsaltii[t,:,:]           = ptopsalti[t,:,:]*1.
                ptopsaltia[t,:,:].mask      = maskAtl
                ptopsaltip[t,:,:].mask      = maskPac
                ptopsaltii[t,:,:].mask      = maskInd
    
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
                # Free memory
                del(persbin, ptopdepth, ptopsigma, ptoptemp, ptopsalt) ; gc.collect()
                # Create cdms2 transient variables
                ptopdepthi  = cdm.createVariable(ptopdepthi,  axes = [timeyr, lati, loni], id = 'ptopdepthi')
                ptopsigmai  = cdm.createVariable(ptopsigmai,  axes = [timeyr, lati, loni], id = 'ptopsigmai')
                ptoptempi   = cdm.createVariable(ptoptempi,   axes = [timeyr, lati, loni], id = 'ptoptempi')
                ptopsalti   = cdm.createVariable(ptopsalti,   axes = [timeyr, lati, loni], id = 'ptopsalti')
                ptopdepthia = cdm.createVariable(ptopdepthia, axes = [timeyr, lati, loni], id = 'ptopdepthia')
                ptopsigmaia = cdm.createVariable(ptopsigmaia, axes = [timeyr, lati, loni], id = 'ptopsigmaia')
                ptoptempia  = cdm.createVariable(ptoptempia,  axes = [timeyr, lati, loni], id = 'ptoptempia')
                ptopsaltia  = cdm.createVariable(ptopsaltia,  axes = [timeyr, lati, loni], id = 'ptopsaltia')
                ptopdepthip = cdm.createVariable(ptopdepthip, axes = [timeyr, lati, loni], id = 'ptopdepthip')
                ptopsigmaip = cdm.createVariable(ptopsigmaip, axes = [timeyr, lati, loni], id = 'ptopsigmaip')
                ptoptempip  = cdm.createVariable(ptoptempip,  axes = [timeyr, lati, loni], id = 'ptoptempip')
                ptopsaltip  = cdm.createVariable(ptopsaltip,  axes = [timeyr, lati, loni], id = 'ptopsaltip')
                ptopdepthii = cdm.createVariable(ptopdepthii, axes = [timeyr, lati, loni], id = 'ptopdepthii')
                ptopsigmaii = cdm.createVariable(ptopsigmaii, axes = [timeyr, lati, loni], id = 'ptopsigmaii')
                ptoptempii  = cdm.createVariable(ptoptempii,  axes = [timeyr, lati, loni], id = 'ptoptempii')
                ptopsaltii  = cdm.createVariable(ptopsaltii,  axes = [timeyr, lati, loni], id = 'ptopsaltii')
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
                # Compute volume/temp/salinity of persistent ocean (global, per basin) (1D)
                persvp = persisti[t,:,:,:]*1. ; persvp[...] = valmask
                persvp = mv.masked_values(persisti[t,:,:,:] >= 99., 1.).mask
                persvp = cdm.createVariable(persvp, axes = [s_axis, lati, loni], id = 'toto')
                persvp._FillValue = valmask ; persvp = maskVal(persvp, valmask)
                  # volume (integral of depth * area)
                voltot         = cdu.averager(thickBini[t,...], axis=0, action = 'sum')
                voltot         = cdu.averager(npy.reshape(voltot*areai ,(Nji*Nii)), action = 'sum')
                volpersxy      = cdu.averager(persvp*thickBini[t,...], axis=0, action = 'sum')
                volpersxy.mask = maski    ; volpersxy  = maskVal(volpersxy , valmask)
                #volpersxya = volpersxy*1.
                #volpersxya.mask = maskAtl ; volpersxya = maskVal(volpersxya, valmask)
                #volpersxyp = volpersxy*1.
                #volpersxyp.mask = maskPac ; volpersxyp = maskVal(volpersxyp, valmask)
                #volpersxyi = volpersxy*1.
                #volpersxyi.mask = maskInd ; volpersxyi = maskVal(volpersxyi, valmask)
                volpersist [t] = cdu.averager(npy.reshape(volpersxy*areai ,(Nji*Nii)), action = 'sum')
                volpersista[t] = cdu.averager(npy.reshape(volpersxy*areaia,(Nji*Nii)), action = 'sum')
                volpersistp[t] = cdu.averager(npy.reshape(volpersxy*areaip,(Nji*Nii)), action = 'sum')
                volpersisti[t] = cdu.averager(npy.reshape(volpersxy*areaii,(Nji*Nii)), action = 'sum')
                  # Temp and salinity (average)
                tempersxy       = cdu.averager(persvp*x1Bini[t,...], axis=0)
                tempersxy.mask  = maski  ; tempersxy  = maskVal(tempersxy , valmask)
                #tempersxya      = tempersxy*1.
                #tempersxya.mask = maskAtl; tempersxya = maskVal(tempersxya, valmask)
                #tempersxyp      = tempersxy*1.
                #tempersxyp.mask = maskPac; tempersxyp = maskVal(tempersxyp, valmask)
                #tempersxyi      = tempersxy*1.
                #tempersxyi.mask = maskInd; tempersxyi = maskVal(tempersxyi, valmask)
                tempersist [t] = cdu.averager(npy.reshape(tempersxy*areai ,(Nji*Nii)), action='sum')/areait
                tempersista[t] = cdu.averager(npy.reshape(tempersxy*areaia,(Nji*Nii)), action='sum')/areaita
                tempersistp[t] = cdu.averager(npy.reshape(tempersxy*areaip,(Nji*Nii)), action='sum')/areaitp
                tempersisti[t] = cdu.averager(npy.reshape(tempersxy*areaii,(Nji*Nii)), action='sum')/areaiti
                salpersxy       = cdu.averager(persvp*x2Bini[t,...], axis=0)
                salpersxy.mask  = maski  ; salpersxy  = maskVal(salpersxy , valmask)
                #salpersxya      = salpersxy*1.
                #salpersxya.mask = maskAtl; salpersxya = maskVal(salpersxya, valmask)
                #salpersxyp      = salpersxy*1.
                #salpersxyp.mask = maskPac; salpersxyp = maskVal(salpersxyp, valmask)
                #salpersxyi      = salpersxy*1.
                #salpersxyi.mask = maskInd; salpersxyi = maskVal(salpersxyi, valmask)
                salpersist [t] = cdu.averager(npy.reshape(salpersxy*areai ,(Nji*Nii)), action='sum')/areait
                salpersista[t] = cdu.averager(npy.reshape(salpersxy*areaia,(Nji*Nii)), action='sum')/areaita
                salpersistp[t] = cdu.averager(npy.reshape(salpersxy*areaip,(Nji*Nii)), action='sum')/areaitp
                salpersisti[t] = cdu.averager(npy.reshape(salpersxy*areaii,(Nji*Nii)), action='sum')/areaiti
                print volpersist.data[t], volpersist.data[t]/voltot.data*100.,tempersist.data[t], salpersist.data[t]
                del(volpersxy,tempersxy,salpersxy)
            #
            # end of loop on t <==
            #
            # Compute % of persistent ocean on the vertical
            persistm                = (cdu.averager(persistv, axis = 1)/cdu.averager(thickBini, axis = 1))
            persistm._FillValue     = valmask
            persistm                = mv.masked_where(persistm > valmask/10, persistm)
            persistm.mask           = maski
            
            # Write zonal mean persistence 3D (time, rhon, lat)
            dbpz   = cdm.createVariable(persistiz , axes = [timeyr, s_axis, lati], id = 'isonpers')
            dbpza  = cdm.createVariable(persistiza, axes = [timeyr, s_axis, lati], id = 'isonpersa')
            dbpzp  = cdm.createVariable(persistizp, axes = [timeyr, s_axis, lati], id = 'isonpersp')
            dbpzi  = cdm.createVariable(persistizi, axes = [timeyr, s_axis, lati], id = 'isonpersi')
            # Write zonal depth of bowl 2D (time, lat)
            dbpdz  = cdm.createVariable(ptopdiz   , axes = [timeyr, lati], id = 'ptopdepth')
            dbpdza = cdm.createVariable(ptopdiza  , axes = [timeyr, lati], id = 'ptopdeptha')
            dbpdzp = cdm.createVariable(ptopdizp  , axes = [timeyr, lati], id = 'ptopdepthp')
            dbpdzi = cdm.createVariable(ptopdizi  , axes = [timeyr, lati], id = 'ptopdepthi')
            # Write zonal density of bowl 2D (time, lat)
            dbprz  = cdm.createVariable(ptopriz   , axes = [timeyr, lati], id = 'ptopsigma')
            dbprza = cdm.createVariable(ptopriza  , axes = [timeyr, lati], id = 'ptopsigmaa')
            dbprzp = cdm.createVariable(ptoprizp  , axes = [timeyr, lati], id = 'ptopsigmap')
            dbprzi = cdm.createVariable(ptoprizi  , axes = [timeyr, lati], id = 'ptopsigmai')
            # Write zonal temperature of bowl 2D (time, lat)
            dbptz  = cdm.createVariable(ptoptiz   , axes = [timeyr, lati], id = 'ptoptemp')
            dbptza = cdm.createVariable(ptoptiza  , axes = [timeyr, lati], id = 'ptoptempa')
            dbptzp = cdm.createVariable(ptoptizp  , axes = [timeyr, lati], id = 'ptoptempp')
            dbptzi = cdm.createVariable(ptoptizi  , axes = [timeyr, lati], id = 'ptoptempi')
            # Write zonal salinity bowl 2D (time, lat)
            dbpsz  = cdm.createVariable(ptopsiz   , axes = [timeyr, lati], id = 'ptopsalt')
            dbpsza = cdm.createVariable(ptopsiza  , axes = [timeyr, lati], id = 'ptopsalta')
            dbpszp = cdm.createVariable(ptopsizp  , axes = [timeyr, lati], id = 'ptopsaltp')
            dbpszi = cdm.createVariable(ptopsizi  , axes = [timeyr, lati], id = 'ptopsalti')
            # Write % of persistent ocean, depth/temp/salinity of bowl 3D (time, lat, lon)
            persim = cdm.createVariable(persistm  , axes = [timeyr, lati, loni], id = 'persim')
            ptopd  = cdm.createVariable(ptopdepthi, axes = [timeyr, lati, loni], id = 'ptopdepth2')
            ptopt  = cdm.createVariable(ptoptempi , axes = [timeyr, lati, loni], id = 'ptoptemp2')
            ptops  = cdm.createVariable(ptopsalti , axes = [timeyr, lati, loni], id = 'ptopsalt2')
            # Write volume/temp/salinity of persistent ocean 1D (time)
            volper  = cdm.createVariable(volpersist *1.e-12, axes = [timeyr], id = 'volpers')
            volpera = cdm.createVariable(volpersista*1.e-12, axes = [timeyr], id = 'volpersa')
            volperp = cdm.createVariable(volpersistp*1.e-12, axes = [timeyr], id = 'volpersp')
            volperi = cdm.createVariable(volpersisti*1.e-12, axes = [timeyr], id = 'volpersi')
            temper  = cdm.createVariable(tempersist        , axes = [timeyr], id = 'tempers')
            tempera = cdm.createVariable(tempersista       , axes = [timeyr], id = 'tempersa')
            temperp = cdm.createVariable(tempersistp       , axes = [timeyr], id = 'tempersp')
            temperi = cdm.createVariable(tempersisti       , axes = [timeyr], id = 'tempersi')
            salper  = cdm.createVariable(salpersist        , axes = [timeyr], id = 'salpers')
            salpera = cdm.createVariable(salpersista       , axes = [timeyr], id = 'salpersa')
            salperp = cdm.createVariable(salpersistp       , axes = [timeyr], id = 'salpersp')
            salperi = cdm.createVariable(salpersisti       , axes = [timeyr], id = 'salpersi')
            if tc == 0:
                # Global attributes
                dbpz.long_name      = 'zonal persistence of isopycnal bins'
                dbpz.units          = '% of time'
                dbpza.long_name     = 'Atl. zonal persistence of isopycnal bins'
                dbpza.units         = '% of time'
                dbpzp.long_name     = 'Pac. zonal persistence of isopycnal bins'
                dbpzp.units         = '% of time'
                dbpzi.long_name     = 'Ind. zonal persistence of isopycnal bins'
                dbpzi.units         = '% of time'
                #
                persim.long_name    = 'Fraction of persistence on isopycnal bins'
                persim.units        = '% of column'
                ptopd.long_name     = 'Depth of shallowest persistent ocean on ison'
                ptopd.units         = 'm'
                ptopt.long_name     = 'Temp. of shallowest persistent ocean on ison'
                ptopt.units         = 'degrees_C'   
                ptops.long_name     = 'Salinity of shallowest persistent ocean on ison'
                ptops.units         = so_h.units
                #
                dbpdz.long_name     = 'Zonal depth of shallowest persistent ocean on ison'
                dbpdz.units         = 'm'
                dbprz.long_name     = 'Zonal rhon of shallowest persistent ocean on ison'
                dbprz.units         = 'sigma_n'
                dbptz.long_name     = 'Zonal Temp. of shallowest persistent ocean on ison'
                dbptz.units         = 'degrees_C'   
                dbpsz.long_name     = 'Zonal Salinity of shallowest persistent ocean on ison'
                dbpsz.units         = so_h.units  
                dbpdza.long_name    = 'Atl. zonal depth of shallowest persistent ocean on ison'
                dbpdza.units        = 'm'
                dbprza.long_name    = 'Atl. Zonal rhon of shallowest persistent ocean on ison'
                dbprza.units        = 'sigma_n'
                dbptza.long_name    = 'Atl. zonal Temp. of shallowest persistent ocean on ison'
                dbptza.units        = 'degrees_C'   
                dbpsza.long_name    = 'Atl. Zonal Salinity of shallowest persistent ocean on ison'
                dbpsza.units        = so_h.units  
                dbpdzp.long_name    = 'Pac. zonal depth of shallowest persistent ocean on ison'
                dbpdzp.units        = 'm'
                dbprzp.long_name    = 'Pac. Zonal rhon of shallowest persistent ocean on ison'
                dbprzp.units        = 'sigma_n'
                dbptzp.long_name    = 'Pac. zonal Temp. of shallowest persistent ocean on ison'
                dbptzp.units        = 'degrees_C'   
                dbpszp.long_name    = 'Pac. zonal Salinity of shallowest persistent ocean on ison'
                dbpszp.units        = so_h.units  
                dbpdzi.long_name    = 'Ind. zonal depth of shallowest persistent ocean on ison'
                dbpdzi.units        = 'm'
                dbprzi.long_name    = 'Ind. Zonal rhon of shallowest persistent ocean on ison'
                dbprzi.units        = 'sigma_n'
                dbptzi.long_name    = 'Ind. zonal Temp. of shallowest persistent ocean on ison'
                dbptzi.units        = 'degrees_C'   
                dbpszi.long_name    = 'Ind. zonal Salinity of shallowest persistent ocean on ison'
                dbpszi.units        = so_h.units  
    
                volper.long_name    = 'Volume of persistent ocean'
                volper.units        = '1.e12 m^3'
                volpera.long_name   = 'Atl. Volume of persistent ocean'
                volpera.units       = '1.e12 m^3'
                volperp.long_name   = 'Pac. Volume of persistent ocean'
                volperp.units       = '1.e12 m^3'
                volperi.long_name   = 'Ind. Volume of persistent ocean'
                volperi.units       = '1.e12 m^3'
                temper.long_name    = 'Temperature of persistent ocean'
                temper.units        = 'degrees_C'
                tempera.long_name   = 'Atl. Temperature of persistent ocean'
                tempera.units       = 'degrees_C'
                temperp.long_name   = 'Pac. Temperature of persistent ocean'
                temperp.units       = 'degrees_C'
                temperi.long_name   = 'Ind. Temperature of persistent ocean'
                temperi.units       = 'degrees_C'
                salper.long_name    = 'Salinity of persistent ocean'
                salper.units        = so_h.units 
                salpera.long_name   = 'Atl. Salinity of persistent ocean'
                salpera.units       = so_h.units 
                salperp.long_name   = 'Pac. Salinity of persistent ocean'
                salperp.units       = so_h.units 
                salperi.long_name   = 'Ind. Salinity of persistent ocean'
                salperi.units       = so_h.units 

            # Write & append
            outFile_f.write(depthbini, extend = 1, index = (trmin-tmin)/12) ; # Write out 4D variable first depth,rhon,lat,lon are written together
            outFile_f.write(thickbini, extend = 1, index = (trmin-tmin)/12) 
            outFile_f.write(x1bini   , extend = 1, index = (trmin-tmin)/12) 
            outFile_f.write(x2bini   , extend = 1, index = (trmin-tmin)/12) 
            outFile_f.write(persim , extend = 1, index = (trmin-tmin)/12) ; # Write out 3D variable depth,lat,lon are written together
            outFile_f.write(ptopd  , extend = 1, index = (trmin-tmin)/12)
            outFile_f.write(ptopt  , extend = 1, index = (trmin-tmin)/12)
            outFile_f.write(ptops  , extend = 1, index = (trmin-tmin)/12)
            outFile_f.write(dbpz   , extend = 1, index = (trmin-tmin)/12)
            outFile_f.write(dbpza  , extend = 1, index = (trmin-tmin)/12)
            outFile_f.write(dbpzp  , extend = 1, index = (trmin-tmin)/12)
            outFile_f.write(dbpzi  , extend = 1, index = (trmin-tmin)/12)
            outFile_f.write(dbpdz  , extend = 1, index = (trmin-tmin)/12)
            outFile_f.write(dbprz  , extend = 1, index = (trmin-tmin)/12)
            outFile_f.write(dbptz  , extend = 1, index = (trmin-tmin)/12)
            outFile_f.write(dbpsz  , extend = 1, index = (trmin-tmin)/12)
            outFile_f.write(dbpdza , extend = 1, index = (trmin-tmin)/12)
            outFile_f.write(dbprza , extend = 1, index = (trmin-tmin)/12)
            outFile_f.write(dbptza , extend = 1, index = (trmin-tmin)/12)
            outFile_f.write(dbpsza , extend = 1, index = (trmin-tmin)/12)
            outFile_f.write(dbpdzp , extend = 1, index = (trmin-tmin)/12)
            outFile_f.write(dbprzp , extend = 1, index = (trmin-tmin)/12)
            outFile_f.write(dbptzp , extend = 1, index = (trmin-tmin)/12)
            outFile_f.write(dbpszp , extend = 1, index = (trmin-tmin)/12)
            outFile_f.write(dbpdzi , extend = 1, index = (trmin-tmin)/12)
            outFile_f.write(dbprzi , extend = 1, index = (trmin-tmin)/12)
            outFile_f.write(dbptzi , extend = 1, index = (trmin-tmin)/12)
            outFile_f.write(dbpszi , extend = 1, index = (trmin-tmin)/12)
            outFile_f.write(volper , extend = 1, index = (trmin-tmin)/12)
            outFile_f.write(volpera, extend = 1, index = (trmin-tmin)/12)
            outFile_f.write(volperp, extend = 1, index = (trmin-tmin)/12)
            outFile_f.write(volperi, extend = 1, index = (trmin-tmin)/12)
            outFile_f.write(temper , extend = 1, index = (trmin-tmin)/12)
            outFile_f.write(tempera, extend = 1, index = (trmin-tmin)/12)
            outFile_f.write(temperp, extend = 1, index = (trmin-tmin)/12)
            outFile_f.write(temperi, extend = 1, index = (trmin-tmin)/12)
            outFile_f.write(salper , extend = 1, index = (trmin-tmin)/12)
            outFile_f.write(salpera, extend = 1, index = (trmin-tmin)/12)
            outFile_f.write(salperp, extend = 1, index = (trmin-tmin)/12)
            outFile_f.write(salperi, extend = 1, index = (trmin-tmin)/12)
            #
            tozp = timc.clock()
            #
            # Init zonal mean output variables
            # Global
            dbz  = cdm.createVariable(depthBinz,        axes = [timeyr, s_axis, lati], id = 'isondepth')
            tbz  = cdm.createVariable(thickBinz,        axes = [timeyr, s_axis, lati], id = 'isonthick')
            vbz  = cdm.createVariable(volBinz*1.e-12,   axes = [timeyr, s_axis, lati], id = 'isonvol')
            x1bz = cdm.createVariable(x1Binz,           axes = [timeyr, s_axis, lati], id = 'thetao')
            x2bz = cdm.createVariable(x2Binz,           axes = [timeyr, s_axis, lati], id = 'so')
            # Atl
            dbza  = cdm.createVariable(depthBinza,      axes = [timeyr, s_axis, lati], id = 'isondeptha')
            tbza  = cdm.createVariable(thickBinza,      axes = [timeyr, s_axis, lati], id = 'isonthicka')
            vbza  = cdm.createVariable(volBinza*1.e-12, axes = [timeyr, s_axis, lati], id = 'isonvola')
            x1bza = cdm.createVariable(x1Binza,         axes = [timeyr, s_axis, lati], id = 'thetaoa')
            x2bza = cdm.createVariable(x2Binza,         axes = [timeyr, s_axis, lati], id = 'soa')
            # Pac
            dbzp  = cdm.createVariable(depthBinzp,      axes = [timeyr, s_axis, lati], id = 'isondepthp')
            tbzp  = cdm.createVariable(thickBinzp,      axes = [timeyr, s_axis, lati], id = 'isonthickp')
            vbzp  = cdm.createVariable(volBinzp*1.e-12, axes = [timeyr, s_axis, lati], id = 'isonvolp')
            x1bzp = cdm.createVariable(x1Binzp,         axes = [timeyr, s_axis, lati], id = 'thetaop')
            x2bzp = cdm.createVariable(x2Binzp,         axes = [timeyr, s_axis, lati], id = 'sop')
            # Ind
            dbzi  = cdm.createVariable(depthBinzi,      axes = [timeyr, s_axis, lati], id = 'isondepthi')
            tbzi  = cdm.createVariable(thickBinzi,      axes = [timeyr, s_axis, lati], id = 'isonthicki')
            vbzi  = cdm.createVariable(volBinzi*1.e-12, axes = [timeyr, s_axis, lati], id = 'isonvoli')
            x1bzi = cdm.createVariable(x1Binzi,         axes = [timeyr, s_axis, lati], id = 'thetaoi')
            x2bzi = cdm.createVariable(x2Binzi,         axes = [timeyr, s_axis, lati], id = 'soi')
            if tc == 0:
                # Global attributes
                dbz.long_name   = 'Global zonal depth of isopycnal'
                dbz.units       = 'm'
                tbz.long_name   = 'Global zonal thickness of isopycnal'
                tbz.units       = 'm'
                vbz.long_name   = 'Volume of isopycnal'
                vbz.units       = '1.e12 m^3'
                x1bz.long_name  = thetao_h.long_name
                x1bz.units      = 'degrees_C'
                x2bz.long_name  = so_h.long_name
                x2bz.units      = so_h.units
                # Atl
                dbza.long_name  = 'Atl. zonal depth of isopycnal'
                dbza.units      = dbz.units
                tbza.long_name  = 'Atl. zonal thickness of isopycnal'
                tbza.units      = 'm'
                vbza.long_name  = 'Atl. volume of isopycnal'
                vbza.units      = '1.e12 m^3'
                x1bza.long_name = thetao_h.long_name
                x1bza.units     = 'degrees_C'
                x2bza.long_name = so_h.long_name
                x2bza.units     = so_h.units
                # Pac
                dbzp.long_name  = 'Pac. zonal depth of isopycnal'
                dbzp.units      = dbz.units
                tbzp.long_name  = 'Pac. zonal thickness of isopycnal'
                tbzp.units      = 'm'
                vbzp.long_name  = 'Pac. volume of isopycnal'
                vbzp.units      = '1.e12 m^3'
                x1bzp.long_name = thetao_h.long_name
                x1bzp.units     = 'degrees_C'
                x2bzp.long_name = so_h.long_name
                x2bzp.units     = so_h.units
                # Ind
                dbzi.long_name  = 'Ind. zonal depth of isopycnal'
                dbzi.units      = dbz.units
                tbzi.long_name  = 'Ind. zonal thickness of isopycnal'
                tbzi.units      = 'm'
                vbzi.long_name  = 'Ind. volume of isopycnal'
                vbzi.units      = '1.e12 m^3'
                x1bzi.long_name = thetao_h.long_name
                x1bzi.units     = 'degrees_C'
                x2bzi.long_name = so_h.long_name
                x2bzi.units     = so_h.units
                # Cleanup
            # Write & append
            outFile_f.write(dbz,   extend = 1, index = (trmin-tmin)/12)
            outFile_f.write(tbz,   extend = 1, index = (trmin-tmin)/12)
            outFile_f.write(vbz,   extend = 1, index = (trmin-tmin)/12)
            outFile_f.write(x1bz,  extend = 1, index = (trmin-tmin)/12)
            outFile_f.write(x2bz,  extend = 1, index = (trmin-tmin)/12)
            # Atl
            outFile_f.write(dbza,  extend = 1, index = (trmin-tmin)/12)
            outFile_f.write(tbza,  extend = 1, index = (trmin-tmin)/12)
            outFile_f.write(vbza,  extend = 1, index = (trmin-tmin)/12)
            outFile_f.write(x1bza, extend = 1, index = (trmin-tmin)/12)
            outFile_f.write(x2bza, extend = 1, index = (trmin-tmin)/12)
            # Pac
            outFile_f.write(dbzp,  extend = 1, index = (trmin-tmin)/12)
            outFile_f.write(tbzp,  extend = 1, index = (trmin-tmin)/12)
            outFile_f.write(vbzp,  extend = 1, index = (trmin-tmin)/12)
            outFile_f.write(x1bzp, extend = 1, index = (trmin-tmin)/12)
            outFile_f.write(x2bzp, extend = 1, index = (trmin-tmin)/12)
            # Ind
            outFile_f.write(dbzi,  extend = 1, index = (trmin-tmin)/12)
            outFile_f.write(tbzi,  extend = 1, index = (trmin-tmin)/12)
            outFile_f.write(vbzi,  extend = 1, index = (trmin-tmin)/12)
            outFile_f.write(x1bzi, extend = 1, index = (trmin-tmin)/12)
            outFile_f.write(x2bzi, extend = 1, index = (trmin-tmin)/12)
    
        # Write/append to file
        if mthout:
            outFileMon_f.write(depthBin, extend = 1, index = trmin-tmin)
            outFileMon_f.write(thickBin, extend = 1, index = trmin-tmin)
            outFileMon_f.write(x1Bin,    extend = 1, index = trmin-tmin)
            outFileMon_f.write(x2Bin,    extend = 1, index = trmin-tmin)
        
        tozf = timc.clock()
        print '   CPU of density bining      =', ticz-tuc
        if tcdel >= 12:
            print '   CPU of annual mean compute =', toz-ticz
            print '   CPU of interpolation       =', tozi-toz
            print '   CPU of zonal mean          =', toziz-tozi
            print '   CPU of persistence compute =', tozp-toziz
        print '   CPU of chunk               =', tozf-tuc
    
    # end loop on tc <===
    print ' [ Time stamp',(timc.strftime("%d/%m/%Y %H:%M:%S")),']'
    print ' Max memory use',resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1.e6,'GB'
    ratio =  12.*float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)/float(grdsize*tmax)
    print ' Ratio to grid*nyears',ratio,'kB/unit(size*nyears)'
    print ' CPU use, elapsed', timc.clock() - ti0, timeit.default_timer() - te0
    ratio = 1.e6*(timc.clock() - ti0)/float(grdsize*tmax)
    print ' Ratio to grid*nyears',ratio,'1.e-6 sec/unit(size*nyears)'
    
    ft.close()
    fs.close()
    if mthout:
        outFileMon_f.close()
        print ' Wrote file: ',outFileMon
    outFile_f.close()
    if tcdel >= 12:
        print ' Wrote file: ',outFile

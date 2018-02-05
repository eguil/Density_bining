'''
 libDensity.py contains all the functions needed by binDensity.py
'''


import os,sys,gc,glob
import cdms2 as cdm
import cdutil as cdu
import MV2 as mv
import numpy as npy
from string import replace
from genutil import statistics
from libToE import findToE
import time as timc


'''
    Define class with all info about input specifications (grid, varname...)
    
    Author:    Nicolas Lebas : nicolas.lebas@locean-ipsl.upmc.fr
    
    Created on Oct 02 11:13:30 2017
    
    Input:
    ------
    - varGrid    - variable with grid infos
    - modeln     - model name
'''
class NCInSpec:
    def __init__(self, varGrid, modeln):
        ## Model
        self.model = modeln
        
        ## Grid
        # Horizontal grid
        self.ingrid  = varGrid.getGrid()
        # Get grid objects
        self.axesList = varGrid.getAxisList()
        # Define dimensions
        self.lon    = varGrid.shape[3]
        self.lat    = varGrid.shape[2]
        self.depth  = varGrid.shape[1]
        # Depth profiles
        z_zt = []
        z_zw = []
        # Read masking value
        try:
            self.valmask = varGrid.missing_value
            if self.valmask is None:
                print 'EC-EARTH missing_value fix'
                self.valmask = 1.e20
        except Exception,err:
            print 'Exception: ',err
            if 'EC-EARTH' == self.model:
                print 'EC-EARTH missing_value fix'
                self.valmask = 1.e20

        ## Variables
        self.thetaoLongName = ''
        self.thetaoUnits = ''
        self.soLongName = ''
        self.soUnits = ''
        self.voLongName = ''
        self.voUnits = ''


'''
Define class with all info about mask for density treatments and outputs
'''
class Mask:
    def __init__(self, gridFile, tgVname):
        gridFile_f  = cdm.open(gridFile)
        maskg       = gridFile_f(tgVname)
        self.outgrid     = maskg.getGrid()
        self.maski       = maskg.mask ; # Global mask
        ### Regional masks ###
        # Atl
        idxa = npy.argwhere(maskg == 1).transpose()
        self.maskAtl = self.maski*1 ; self.maskAtl[...] = True
        self.maskAtl[idxa[0],idxa[1]] = False
        # Pac
        self.maskPac = self.maski*1 ; self.maskPac[...] = True
        idxp = npy.argwhere(maskg == 2).transpose()
        self.maskPac[idxp[0],idxp[1]] = False
        # Ind
        self.maskInd = self.maski*1 ; self.maskInd[...] = True
        idxi = npy.argwhere(maskg == 3).transpose()
        self.maskInd[idxi[0],idxi[1]] = False
        
        self.lonI = maskg.getLongitude()
        self.latI = maskg.getLatitude()
        gridFile_f.close()


'''
Define class with all info about area for density treatments and outputs
'''
class Area:
    def __init__(self, fileArea, loni, lati, masks):
        ff = cdm.open(fileArea)
        self.area    = ff('areacello')
        ff.close()
        
        Nii     = len(loni)
        Nji     = len(lati)
        # Compute area of target grid and zonal and global sums
        self.areai = computeArea(loni[:], lati[:])
        self.areai.mask = masks['glob']
        self.areaia = self.areai*1. ; self.areaia.mask = masks['Atl']
        self.areaip = self.areai*1. ; self.areaip.mask = masks['Pac']
        self.areaii = self.areai*1. ; self.areaii.mask = masks['Ind']
        self.areazt  = npy.ma.sum(self.areai , axis=1)
        self.areazta = npy.ma.sum(self.areaia, axis=1)
        self.areaztp = npy.ma.sum(self.areaip, axis=1)
        self.areazti = npy.ma.sum(self.areaii, axis=1)
        self.areait  = npy.ma.sum(npy.reshape(self.areai ,(Nji*Nii)))
        self.areaita = npy.ma.sum(npy.reshape(self.areaia,(Nji*Nii)))
        self.areaitp = npy.ma.sum(npy.reshape(self.areaip,(Nji*Nii)))
        self.areaiti = npy.ma.sum(npy.reshape(self.areaii,(Nji*Nii)))


def createAxisRhoBassin(s_sax):
    rhoAxis                 = cdm.createAxis(s_sax,bounds=None,id='lev')
    rhoAxis.positive        = 'down'
    rhoAxis.long_name       = 'ocean neutral density coordinate'
    rhoAxis.standard_name   = 'lev'
    rhoAxis.units           = 'kg m-3'
    rhoAxis.units_long      = 'kg m-3 (anomaly, minus 1000)'
    rhoAxis.axis            = 'Z'
    rhoAxis.designateLevel()
    # Define basin output axis
    basinAxis               = cdm.createAxis([0,1,2,3],bounds=None,id='basin')
    basinAxis.long_name     = 'ocean basin index'
    basinAxis.standard_name = 'basin'
    basinAxis.units         = 'basin index'
    basinAxis.units_long    = '0: global_ocean 1: atlantic_ocean; 2: pacific_ocean; 3: indian_ocean'
    basinAxis.axis          = 'B'
    # Create rho axis list
    #    rhoAxesList             = [axesList[0],rhoAxis,axesList[2],axesList[3]] ; # time, rho, lat, lon
    # Create basin-zonal axes lists
    #    basinTimeList           = [axesList[0],basinAxis] ; # time, basin
    #    basinAxesList           = [axesList[0],basinAxis,axesList[2]] ; # time, basin, lat
    #    basinRhoAxesList        = [axesList[0],basinAxis,rhoAxis,axesList[2]] ; # time, basin, rho, lat
    return (rhoAxis, basinAxis)


def maskVal(field,valmask):
    '''
    The maskVal() function applies a mask to an array provided
    
    Author:    Eric Guilyardi : Eric.Guilyardi@locean-ipsl.upmc.fr
    Co-author: Paul J. Durack : pauldurack@llnl.gov : @durack1.
    
    Created on Sun Sep 14 21:13:30 2014

    Inputs:
    ------
    - field     - 1D/2D/3D array
    - valmask    - 1D scalar of mask value
    
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
    lonN = int(lon.shape[0])
    latN = int(lat.shape[0])
    area = npy.ma.ones([latN, lonN], dtype='float32')*0.
    lonr = lon[:] * radconv
    latr = lat[:] * radconv
    #loop
    for i in range(1,lonN-1):
        lonm1 = (lonr[i-1] + lonr[i]  )*0.5
        lonp1 = (lonr[i]   + lonr[i+1])*0.5
        for j in range(1,latN-1):
            latm1 = (latr[j-1] + latr[j]  )*0.5
            latp1 = (latr[j]   + latr[j+1])*0.5
            area[j,i] = npy.float(radius**2 * (lonp1 - lonm1) * (npy.sin(latp1) - npy.sin(latm1)))
        # North and south bounds
        latm1 = ((-90.*radconv) + latr[0] )*0.5
        latp1 = (latr[0]        + latr[1] )*0.5
        area[0,i] = npy.float(radius**2 * (lonp1 - lonm1) * (npy.sin(latp1) - npy.sin(latm1)))
        latm1 = (latr[latN-2] + latr[latN-1])*0.5
        latp1 = (latr[latN-1] + (90.*radconv)  )*0.5
        area[latN-1,i] = npy.float(radius**2 * (lonp1 - lonm1) * (npy.sin(latp1) - npy.sin(latm1)))
    # East and west bounds
    area[:,0]     = area[:,1]
    area[:,lonN-1] = area[:,lonN-2]

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
    - PJD 14 Sep 2014 - Rewrote as function
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

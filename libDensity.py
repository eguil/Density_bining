import os,sys,gc,glob
import cdms2 as cdm
import cdutil as cdu
import MV2 as mv
import numpy as npy
from string import replace
from genutil import statistics
from libToE import findToE
import time as timc


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


def mmeAveMsk2D(listFiles, years, inDir, outDir, outFile, timeInt, mme, timeBowl, ToeType, debug=True):
    '''
    The mmeAveMsk2D() function averages rhon/lat density bined files with differing masks
    It ouputs
     - the MME
     - a percentage of non-masked bins
     - the sign agreement of period2-period1 differences
     - ToE per run and for MME

    Author:    Eric Guilyardi : Eric.Guilyardi@locean-ipsl.upmc.fr

    Created on Tue Nov 25 13:56:20 CET 2014

    Inputs:
    -------
    - listFiles(str)         - the list of files to be averaged
    - years(t1,t2)           - years for slice read
    - inDir[](str)           - input directory where files are stored (add histnat as inDir[1] for ToE)
    - outDir(str)            - output directory
    - outFile(str)           - output file
    - timeInt(2xindices)     - indices of init period to compare with (e.g. [1,20])
    - mme(bool)              - multi-model mean (will read in single model ensemble stats)
    - timeBowl               - either time 'mean' or time 'max' bowl used to mask out bowl
    - ToeType(str)           - ToE type ('F': none, 'histnat')
                               -> requires running first mm+mme without ToE to compute Stddev
    - debug <optional>       - boolean value

    Notes:
    -----
    - EG 25 Nov 2014   - Initial function write
    - EG 27 Nov 2014   - Rewrite with loop on variables
    - EG 06 Dec 2014   - Added agreement on difference with init period - save as <var>Agree
    - EG 07 Dec 2014   - Read bowl to remove points above bowl - save as <var>Bowl
    - EG 19 Apr 2016   - ToE computation (just for 2D files)
    - EG 07 Oct 2016   - add 3D file support
    - EG 21 Nov 2016   - move 3D support to new function
    - EG 10 jan 2017   - added timeBowl option

    - TODO :
                 - remove loops
                 - add computation of ToE per model (toe 1 and toe 2) see ticket #50
                 - add isonhtc (see ticket #48)
    '''

    # CDMS initialisation - netCDF compression
    comp = 1 # 0 for no compression
    cdm.setNetcdfShuffleFlag(comp)
    cdm.setNetcdfDeflateFlag(comp)
    cdm.setNetcdfDeflateLevelFlag(comp)
    cdm.setAutoBounds('on')
    # Numpy initialisation
    npy.set_printoptions(precision=2)

    if debug:
        debug = True
    else:
        debug = False
    # File dim and grid inits
    t1 = years[0]
    t2 = years[1]
    if t2 <= 0:
        useLastYears = True
        t2 = -t2
    else:
        useLastYears = False
    t10 = t1
    t20 = t2
    # Bound of period average to remove
    peri1 = timeInt[0]
    peri2 = timeInt[1]
    fi      = cdm.open(inDir[0]+'/'+listFiles[0])
    isond0  = fi['isondepth'] ; # Create variable handle
    # Get grid objects
    axesList = isond0.getAxisList()
    sigmaGrd = isond0.getLevel()
    latN = isond0.shape[3]
    levN = isond0.shape[2]
    basN = isond0.shape[1]
    varsig='ptopsigma'

    # Declare and open files for writing
    if os.path.isfile(outDir+'/'+outFile):
        os.remove(outDir+'/'+outFile)
    outFile_f = cdm.open(outDir+'/'+outFile,'w')

    # Testing mme with less models
    #listFiles=listFiles[0:4]

    #timN = isond0.shape[0]
    timN = t2-t1
    runN = len(listFiles)

    print ' Number of members:',len(listFiles)

    valmask = isond0.missing_value
    varList = ['isondepth','isonpers','isonso','isonthetao','isonthick','isonvol']
    varFill = [0.,0.,valmask,valmask,0.,0.]
    # init arrays (2D rho/lat)
    percent  = npy.ma.ones([runN,timN,basN,levN,latN], dtype='float32')*0.
    #minbowl  = npy.ma.ones([basN,latN], dtype='float32')*1000.
    varbowl  = npy.ma.ones([runN,timN,basN,latN], dtype='float32')*1.
    #varList = ['isondepth']
    #print ' !!! ### Testing one variable ###'
    #varList = ['isonthetao']

    # init time axis
    time       = cdm.createAxis(npy.float32(range(timN)))
    time.id    = 'time'
    time.units = 'years since 1861'
    time.designateTime()
    # init ensemble axis
    ensembleAxis       = cdm.createAxis(npy.float32(range(runN)))
    ensembleAxis.id    = 'members'
    ensembleAxis.units = 'N'

    # loop on variables
    for iv,var in enumerate(varList):

        # Array inits (2D rho/lat 3D rho/lat/lon)
            #shapeR = [basN,levN,latN]
        isonvar  = npy.ma.ones([runN,timN,basN,levN,latN], dtype='float32')*valmask
        vardiff,varbowl2D = [npy.ma.ones(npy.ma.shape(isonvar)) for _ in range(2)]
        varstd,varToE1,varToE2 =  [npy.ma.ones([runN,basN,levN,latN], dtype='float32')*valmask for _ in range(3)]
        varones  = npy.ma.ones([runN,timN,basN,levN,latN], dtype='float32')*1.

        print ' Variable ',iv, var
        # loop over files to fill up array
        for i,file in enumerate(listFiles):
            ft      = cdm.open(inDir[0]+'/'+file)
            model = file.split('.')[1]
            timeax  = ft.getAxis('time')
            file1d = replace(inDir[0]+'/'+file,'2D','1D')
            if os.path.isfile(file1d):
                f1d = cdm.open(file1d)
            else:
                print 'ERROR:',file1d,'missing (if mme, run 1D first)'
                sys.exit(1)
            tmax = timeax.shape[0]
            if i == 0:
                tmax0 = tmax
            #adapt [t1,t2] time bounds to piControl last NN years
            if useLastYears:
                t1 = tmax-t20
                t2 = tmax
            else:
                if tmax != tmax0:
                    print 'wrong time axis: exiting...'
                    return

            # read array
            # loop over time/density for memory management
            for it in range(timN):
                t1r = t1 + it
                t2r = t1r + 1
                isonRead = ft(var,time = slice(t1r,t2r))
                if varFill[iv] != valmask:
                    isonvar[i,it,...] = isonRead.filled(varFill[iv])
                else:
                    isonvar[i,it,...] = isonRead
            # compute percentage of non-masked points accros MME
            if iv == 0:
                maskvar = mv.masked_values(isonRead.data,valmask).mask
                percent[i,...] = npy.float32(npy.equal(maskvar,0))
            if mme:
                # if mme then just accumulate Bowl, Agree fields
                varst = var+'Agree'
                vardiff[i,...] = ft(varst,time = slice(t1,t2))
                varb = var+'Bowl'
                varbowl2D[i,...] = ft(varb,time = slice(t1,t2))
            else:
                # Compute difference with average of first initN years
                varinit = cdu.averager(isonvar[i,peri1:peri2,...],axis=0)
                for t in range(timN):
                    vardiff[i,t,...] = isonvar[i,t,...] - varinit
                vardiff[i,...].mask = isonvar[i,...].mask
                # Read bowl and truncate 2D field above bowl
                if iv == 0:
                    bowlRead = f1d(varsig,time = slice(t1,t2))
                    varbowl[i,...] = bowlRead
                # Compute Stddev
                varstd[i,...] = npy.ma.std(isonvar[i,...], axis=0)
                # Compute ToE
                if ToeType == 'histnat':
                    # Read mean and Std dev from histnat
                    if i == 0:
                        filehn  = glob.glob(inDir[1]+'/cmip5.'+model+'.*zon2D*')[0]
                        #filehn = replace(outFile,'historical','historicalNat')
                        fthn = cdm.open(filehn)
                        varmeanhn = fthn(var)
                        varst = var+'Std'
                        varmaxstd = fthn(varst)
                    toemult = 1.
                    signal = npy.reshape(isonvar[i,...]-varmeanhn,(timN,basN*levN*latN))
                    noise = npy.reshape(varmaxstd,(basN*levN*latN))
                    varToE1[i,...] = npy.reshape(findToE(signal, noise, toemult),(basN,levN,latN))
                    toemult = 2.
                    varToE2[i,...] = npy.reshape(findToE(signal, noise, toemult),(basN,levN,latN))
            ft.close()
            f1d.close()
        # <-- end of loop on files

        # Compute percentage of bin presence
        # Only keep points where percent > 50%
        if iv == 0:
            percenta = (cdu.averager(percent,axis=0))*100.
            percenta = mv.masked_less(percenta, 50)
            percentw = cdm.createVariable(percenta, axes = [time,axesList[1],axesList[2],axesList[3]], id = 'isonpercent')
            percentw._FillValue = valmask
            percentw.long_name = 'percentage of MME bin'
            percentw.units     = '%'
            outFile_f.write(percentw.astype('float32'))

        # Sign of difference
        if mme:
            vardiffsgSum = cdu.averager(vardiff, axis=0)
            vardiffsgSum = cdm.createVariable(vardiffsgSum , axes =[time,axesList[1],axesList[2],axesList[3]] , id = 'foo')
            vardiffsgSum = maskVal(vardiffsgSum, valmask)
            vardiffsgSum.mask = percentw.mask
        else:
            vardiffsg = npy.copysign(varones,vardiff)
            # average signs
            vardiffsgSum = cdu.averager(vardiffsg, axis=0)
            vardiffsgSum = mv.masked_greater(vardiffsgSum, 10000.)
            vardiffsgSum.mask = percentw.mask
            vardiffsgSum._FillValue = valmask

        # average variable accross members
        isonVarAve = cdu.averager(isonvar, axis=0)
        isonVarAve = cdm.createVariable(isonVarAve , axes =[time,axesList[1],axesList[2],axesList[3]] , id = 'foo')
        # mask
        if varFill[iv] == valmask:
            isonVarAve = maskVal(isonVarAve, valmask)

        isonVarAve.mask = percentw.mask

        # Only keep points with rhon >  bowl-delta_rho
        delta_rho = 0.
        if mme: # start from average of <var>Agree
            isonVarBowl = cdu.averager(varbowl2D, axis=0)
            isonVarBowl = cdm.createVariable(isonVarBowl , axes =[time,axesList[1],axesList[2],axesList[3]] , id = 'foo')
            isonVarBowl = maskVal(isonVarBowl, valmask)
            isonVarBowl.mask = percentw.mask
            # Compute intermodel stddev
            isonVarStd = statistics.std(varbowl2D, axis=0)
            isonVarStd = cdm.createVariable(isonVarStd , axes =[time,axesList[1],axesList[2],axesList[3]] , id = 'foo')
            isonVarStd = maskVal(isonVarStd, valmask)
            isonVarStd.mask = percentw.mask
            if iv == 0:
                # Read mulitmodel sigma on bowl and average in time
                file1d  =  replace(outDir+'/'+outFile,'2D','1D')
                if os.path.isfile(file1d):
                    f1d = cdm.open(file1d)
                else:
                    print 'ERROR:',file1d,'missing (if mme, run 1D first)'
                    sys.exit(1)
                bowlRead = f1d(varsig,time = slice(t1,t2))
                f1d.close()
                siglimit = cdu.averager(bowlRead, axis=0)  - delta_rho
            # TODO: remove loop by building global array with 1/0
            for il in range(latN):
                for ib in range(basN):
                    #if ib == 2:
                    #    print il, siglimit[ib,il]
                    if siglimit[ib,il] < valmask/1000.:
                         # if mme bowl density defined, mask above bowl
                        index = (npy.argwhere(sigmaGrd[:] >= siglimit[ib,il]))
                        isonVarBowl [:,ib,0:index[0],il].mask = True
                        isonVarStd  [:,ib,0:index[0],il].mask = True
                        vardiffsgSum[:,ib,0:index[0],il].mask = True
                    else:
                        # mask all points
                        isonVarBowl [:,ib,:,il].mask = True
                        isonVarStd  [:,ib,:,il].mask = True
                        vardiffsgSum[:,ib,:,il].mask = True
        else:
            isonVarBowl = isonVarAve*1. # start from variable
            isonVarStd  = isonVarAve*1. # start from variable
            if iv == 0:
                siglimit = cdu.averager(varbowl, axis=0) # average accross members
                # Average bowl in time
                if timeBowl == 'mean':
                    siglimit = cdu.averager(siglimit, axis=0) - delta_rho
                # or take largest sigma over time
                else:
                    siglimit = npy.ma.max(siglimit, axis=0) - delta_rho
            # TODO: remove loop by building global array with 1/0
            for il in range(latN):
                for ib in range(basN):
                    if siglimit[ib,il] < valmask/1000.:
                        # if bowl density defined, mask above bowl
                        index = (npy.argwhere(sigmaGrd[:] >= siglimit[ib,il]))
                        isonVarBowl[:,ib,0:index[0],il].mask = True
                        vardiffsgSum[:,ib,0:index[0],il].mask = True
                    else:
                        # mask all points
                        vardiffsgSum[:,ib,:,il].mask = True

            isonVarBowl = maskVal(isonVarBowl, valmask)
            # Find max of Std dev of all members
            isonVarStd = npy.ma.max(varstd, axis=0)
            # mask
            if varFill[iv] == valmask:
                isonVarStd = maskVal(isonVarStd, valmask)

        # Write
        isonave = cdm.createVariable(isonVarAve, axes = [time,axesList[1],axesList[2],axesList[3]], id = isonRead.id)
        isonave.long_name = isonRead.long_name
        isonave.units     = isonRead.units
        isonavediff = cdm.createVariable(vardiffsgSum, axes = [time,axesList[1],axesList[2],axesList[3]], id = isonRead.id+'Agree')
        isonavediff.long_name = isonRead.long_name
        isonavediff.units     = isonRead.units
        isonavebowl = cdm.createVariable(isonVarBowl, axes = [time,axesList[1],axesList[2],axesList[3]], id = isonRead.id+'Bowl')
        isonavebowl.long_name = isonRead.long_name
        isonavebowl.units     = isonRead.units
        if not mme:
            isonmaxstd = cdm.createVariable(isonVarStd, axes = [axesList[1],axesList[2],axesList[3]], id = isonRead.id+'Std')
            isonmaxstd.long_name = isonRead.long_name
            isonmaxstd.units     = isonRead.units

        outFile_f.write(    isonave.astype('float32'))
        outFile_f.write(isonavediff.astype('float32'))
        outFile_f.write(isonavebowl.astype('float32'))
        if not mme:
            outFile_f.write( isonmaxstd.astype('float32'))

        if ToeType == 'histnat':
            isontoe1 = cdm.createVariable(varToE1, axes = [ensembleAxis,axesList[1],axesList[2],axesList[3]], id = isonRead.id+'ToE1')
            isontoe1.long_name = 'ToE 1 for '+isonRead.long_name
            isontoe1.units     = 'Year'
            isontoe2 = cdm.createVariable(varToE2, axes = [ensembleAxis,axesList[1],axesList[2],axesList[3]], id = isonRead.id+'ToE2')
            isontoe2.long_name = 'ToE 2 for '+isonRead.long_name
            isontoe2.units     = 'Year'
            outFile_f.write(isontoe1.astype('float32'))
            outFile_f.write(isontoe2.astype('float32'))

        if mme:
            isonvarstd = cdm.createVariable(isonVarStd , axes =[time,axesList[1],axesList[2],axesList[3]] , id = isonRead.id+'ModStd')
            isonvarstd.long_name = isonRead.long_name+' intermodel std'
            isonvarstd.units     = isonRead.units
            outFile_f.write(isonvarstd.astype('float32'))

    # <--- end of loop on variables 

    outFile_f.close()
    fi.close()

def mmeAveMsk3D(listFiles, years, inDir, outDir, outFile, timeInt, mme, ToeType, debug=True):
    '''
    The mmeAveMsk3D() function averages rhon/lat density bined files with differing masks
    It ouputs
     - the MME
     - a percentage of non-masked bins
     - the sign agreement of period2-period1 differences
     - ToE per run and for MME

    Author:    Eric Guilyardi : Eric.Guilyardi@locean-ipsl.upmc.fr

    Created on Tue Nov 21 2016

    Inputs:
    -------
    - listFiles(str)         - the list of files to be averaged
    - years(t1,t2)           - years for slice read
    - inDir[](str)           - input directory where files are stored (add histnat as inDir[1] for ToE)
    - outDir(str)            - output directory
    - outFile(str)           - output file
    - timeInt(2xindices)     - indices of init period to compare with (e.g. [1,20])
    - mme(bool)              - multi-model mean (will read in single model ensemble stats)
    - ToeType(str)           - ToE type ('F': none, 'histnat')
                               -> requires running first mm+mme without ToE to compute Stddev
    - debug <optional>       - boolean value

    Notes:
    -----
    - EG 21 Nov 2016   - Initial function write

    - TODO :
                 - add computation of ToE per model (toe 1 and toe 2) see ticket #50
                 - add isonhtc (see ticket #48)
    '''

    # CDMS initialisation - netCDF compression
    comp = 1 # 0 for no compression
    cdm.setNetcdfShuffleFlag(comp)
    cdm.setNetcdfDeflateFlag(comp)
    cdm.setNetcdfDeflateLevelFlag(comp)
    cdm.setAutoBounds('on')
    # Numpy initialisation
    npy.set_printoptions(precision=2)

    if debug:
        debug = True
    else:
        debug = False
    # File dim and grid inits
    t1 = years[0]
    t2 = years[1]
    # Bound of period average to remove
    peri1 = timeInt[0]
    peri2 = timeInt[1]
    fi    = cdm.open(inDir[0]+'/'+listFiles[0])
    # Switch if only variables below the bowl are present/treated
    nobowl = True
    if nobowl:
        isond0 = fi['isondepthgBowl'] ; # Create variable handle
    else:
        isond0 = fi['isondepthg'] ; # Create variable handle
    # Get grid objects
    axesList = isond0.getAxisList()
    sigmaGrd = isond0.getLevel()
    #time = isond0.getTime()
    lonN = isond0.shape[3]
    latN = isond0.shape[2]
    levN = isond0.shape[1]
    varsig='ptopsigmaxy'

    # Limit number of models to 3 for testing of mme
    #if mme:
    #    listFiles = listFiles[0:2]
    #    print ' !!! ### Testing 3 models ###',  listFiles

    # Declare and open files for writing
    if os.path.isfile(outDir+'/'+outFile):
        os.remove(outDir+'/'+outFile)
    outFile_f = cdm.open(outDir+'/'+outFile,'w')

    #timN = isond0.shape[0]
    timN = t2-t1
    runN = len(listFiles)

    print ' Number of members:',len(listFiles)

    valmask = isond0.missing_value

    varList = ['isondepthg','persistmxy','sog','thetaog','isonthickg']
    varFill = [valmask,valmask,valmask,valmask,valmask]
    percent  = npy.ma.ones([runN,timN,latN,lonN], dtype='float32')*0.
    varbowl  = npy.ma.ones([runN,timN,latN,lonN], dtype='float32')*1.
    #varList = ['isondepthg']
    #print ' !!! ### Testing one variable ###', varList

    # init sigma axis
    sigma = cdm.createAxis(npy.float32(range(1)))
    sigma.id = axesList[1].id
    sigma.units = axesList[1].units
    sigma.designateTime()
    # init time axis
    time       = cdm.createAxis(npy.float32(range(timN)))
    time.id    = 'time'
    time.units = 'years since 1861'
    # init ensemble axis
    ensembleAxis       = cdm.createAxis(npy.float32(range(runN)))
    ensembleAxis.id    = 'members'
    ensembleAxis.units = 'N'
    # Output axis
    sigmaList = [sigma,axesList[2],axesList[3]] ; # sigma, lat, lon
    sigmaTimeList = [sigma,time,axesList[2],axesList[3]] ; # sigma, time, lat, lon
    # init arrays
    isonvar  = npy.ma.ones([runN,timN,latN,lonN], dtype='float32')*valmask
    varbowl2D  = npy.ma.ones([runN,timN,latN,lonN], dtype='float32')*valmask
    varstd,varToE1,varToE2 =  [npy.ma.ones([runN,latN,lonN], dtype='float32')*valmask for _ in range(3)]

    # Loop on density levels (for memory management, becomes UNLIMITED axis and requires a ncpq to reorder dimensions)

    delta_ib = 1
    print ' Sigma index:'
    for ib in range(levN):
        ib1 = ib + delta_ib
        print ib,
        tim0 = timc.clock()
        # loop on variables
        for iv,var in enumerate(varList):
            if nobowl:
                varb = var+'Bowl'
            else:
                varb = var
            if ib == 0:
                print ' Variable ',iv, varb
            # loop over files to fill up array
            for i,file in enumerate(listFiles):
                tim01 = timc.clock()
                ft      = cdm.open(inDir[0]+'/'+file)
                model = file.split('.')[1]
                timeax  = ft.getAxis('time')
                if i == 0:
                    tmax0 = timeax.shape[0]
                tmax = timeax.shape[0]
                if tmax != tmax0:
                    print 'wrong time axis: exiting...'
                    return
                # read array
                isonRead = ft(varb,time = slice(t1,t2), lev = slice(ib,ib1)).squeeze()
                if varFill[iv] != valmask:
                    isonvar[i,...] = isonRead.filled(varFill[iv])
                else:
                    isonvar[i,...] = isonRead
                tim02 = timc.clock()
                # compute percentage of non-masked points accros MME
                if iv == 0:
                    maskvar = mv.masked_values(isonRead.data,valmask).mask
                    percent[i,...] = npy.float32(npy.equal(maskvar,0))
                tim03 = timc.clock()
                if mme:
                    # if mme then just accumulate Bowl, Agree and Std fields
                    #varst = var+'Agree'
                    #vardiff[i,...] = ft(varst,time = slice(t1,t2),lev = slice(ib,ib1)).squeeze()
                    isonRead = ft(varb,time = slice(t1,t2),lev = slice(ib,ib1)).squeeze()
                    varbowl2D[i,...] = isonRead
                else:
                    # Compute difference with average of first initN years
                    #varinit = cdu.averager(isonvar[i,peri1:peri2,...],axis=0)
                    #for t in range(timN):
                    #    vardiff[i,t,...] = isonvar[i,t,...] - varinit
                    #vardiff[i,...].mask = isonvar[i,...].mask
                    # Read bowl to truncate field above bowl
                    if ib == 0 and iv == 0:
                        varbowl[i,...] = ft(varsig,time = slice(t1,t2))
                        #varbowl[i,...] = bowlRead
                    # Compute Stddev
                    varstd[i,...] = npy.ma.std(isonvar[i,...], axis=0)
                    # Compute ToE
                    if ToeType == 'histnat':
                        toto=1
                        # TODO
                        # Read mean and Std dev from histnat
                        #    if i == 0:
                        #        filehn  = glob.glob(inDir[1]+'/cmip5.'+model+'.*zon2D*')[0]
                        #        #filehn = replace(outFile,'historical','historicalNat')
                        #        fthn = cdm.open(filehn)
                        #        varmeanhn = fthn(var)
                        #        varst = var+'Std'
                        #        varmaxstd = fthn(varst)
                        #    toemult = 1.
                        #    signal = npy.reshape(isonvar[i,...]-varmeanhn,(timN,basN*levN*latN))
                        #    noise = npy.reshape(varmaxstd,(basN*levN*latN))
                        #    varToE1[i,...] = npy.reshape(findToE(signal, noise, toemult),(basN,levN,latN))
                        #    toemult = 2.
                        #    varToE2[i,...] = npy.reshape(findToE(signal, noise, toemult),(basN,levN,latN))
                tim04 = timc.clock()
                ft.close()
                #print 'ib, section 1 timing',ib, tim02-tim01,tim03-tim02,tim04-tim03
            # <-- end of loop on files (i)

            tim1 = timc.clock()

            # Compute percentage of bin presence
            # Only keep points where percent > 50%
            if iv == 0:
                percenta = (cdu.averager(percent,axis=0))*100.
                percenta = mv.masked_less(percenta, 50)
                percenta = npy.reshape(percenta,[delta_ib,timN,latN,lonN])
                percentw = cdm.createVariable(percenta, axes = sigmaTimeList, id = 'isonpercent')
                percentw._FillValue = valmask
                percentw.long_name = 'percentage of MME bin'
                percentw.units     = '%'
                outFile_f.write(percentw.astype('float32'), extend = 1, index = ib)

            # Sign of difference
            #if mme:
            #    vardiffsgSum = cdu.averager(vardiff, axis=0)
            #    vardiffsgSum = cdm.createVariable(vardiffsgSum , axes = sigmaTimeList , id = 'foo')
            #    vardiffsgSum = maskVal(vardiffsgSum, valmask)
            #    vardiffsgSum.mask = percentw.mask
            #else:
            #    vardiffsg = npy.copysign(varones,vardiff)
            #    # average signs
            #    vardiffsgSum = cdu.averager(vardiffsg, axis=0)
            #    vardiffsgSum = mv.masked_greater(vardiffsgSum, 10000.)
            #    vardiffsgSum.mask = percentw.mask
            #    vardiffsgSum._FillValue = valmask

            # average variable accross members
            isonVarAve = cdu.averager(isonvar, axis=0)
            isonVarAve = npy.reshape(isonVarAve,[delta_ib,timN,latN,lonN])
            isonVarAve = cdm.createVariable(isonVarAve , axes = sigmaTimeList , id = 'foo')
            # mask
            if varFill[iv] == valmask:
                isonVarAve = maskVal(isonVarAve, valmask)

            isonVarAve.mask = percentw.mask
            tim2 = timc.clock()

            # Only keep points with rhon >  bowl-delta_rho
            delta_rho = 0.
            # mme case
            if mme: # start from average of <var>Agree
                isonVarBowl = cdu.averager(varbowl2D, axis=0)
                isonVarBowl = npy.reshape(isonVarBowl,[delta_ib,timN,latN,lonN])
                isonVarBowl = cdm.createVariable(isonVarBowl , axes = sigmaTimeList , id = 'foo')
                isonVarBowl = maskVal(isonVarBowl, valmask)
                isonVarBowl.mask = percentw.mask
                # Compute intermodel stddev
                isonVarStd = statistics.std(varbowl2D, axis=0)
                isonVarStd = npy.reshape(isonVarStd,[delta_ib,timN,latN,lonN])
                isonVarStd = cdm.createVariable(isonVarStd , axes = sigmaTimeList , id = 'foo')
                isonVarStd = maskVal(isonVarStd, valmask)
                isonVarStd.mask = percentw.mask

                # Write
                isonvarbowlw = cdm.createVariable(isonVarBowl , axes = sigmaTimeList , id = isonRead.id)
                isonvarbowlw.long_name = isonRead.long_name
                isonvarbowlw.units     = isonRead.units
                isonvarstdw = cdm.createVariable(isonVarStd , axes = sigmaTimeList , id = isonRead.id+'Std')
                isonvarstdw.long_name = isonRead.long_name+' intermodel std'
                isonvarstdw.units     = isonRead.units

                outFile_f.write(isonvarbowlw.astype('float32'), extend = 1, index = ib)
                outFile_f.write(isonvarstdw.astype('float32'), extend = 1, index = ib)

                #if ib == 0 and iv == 0:
                #    # TODO review
                #    # Read multimodel sigma on bowl and average in time
                #    file1d  =  replace(outDir+'/'+outFile,'2D','1D')
                #    if os.path.isfile(file1d):
                #        f1d = cdm.open(file1d)
                #    else:
                #        print 'ERROR:',file1d,'missing (if mme, run 2D first)'
                #        sys.exit(1)
                #    bowlRead = f1d(varsig,time = slice(t1,t2),lev = slice(ib,ib1))
                #    f1d.close()
                #    siglimit = cdu.averager(bowlRead, axis=0)  - delta_rho
                # TODO: remove loop by building global array with 1/0
                #if sw2d == 1:
                #    for il in range(latN):
                #        for ib in range(basN):
                #            #if ib == 2:
                #            #    print il, siglimit[ib,il]
                #            if siglimit[ib,il] < valmask/1000.:
                #                 # if mme bowl density defined, mask above bowl
                #                index = (npy.argwhere(sigmaGrd[:] >= siglimit[ib,il]))
                #                isonVarBowl [:,ib,0:index[0],il].mask = True
                #                isonVarStd  [:,ib,0:index[0],il].mask = True
                #                vardiffsgSum[:,ib,0:index[0],il].mask = True
                #            else:
                #                # mask all points
                #                isonVarBowl [:,ib,:,il].mask = True
                #                isonVarStd  [:,ib,:,il].mask = True
                #                vardiffsgSum[:,ib,:,il].mask = True
            # mm case
            else:
                isonVarBowl = isonVarAve*1. # start from variable
                #isonVarStd  = isonVarAve*1. # start from variable
                if ib == 0 and iv == 0:
                    # build bowl position
                    siglimit = cdu.averager(varbowl, axis=0) # average accross members
                    siglimit = npy.reshape(siglimit,[timN*latN*lonN]) - delta_rho
                if iv == 0:
                    sigarr = siglimit*1.
                    sigarr[:] = sigmaGrd[ib]
                # test
                i = 60
                j = 60
                ij = j*lonN+i
                isonVarBowl = npy.reshape(isonVarBowl,[timN*latN*lonN])
                #vardiffsgSum = npy.reshape(vardiffsgSum,[timN*latN*lonN])

                isonVarBowl.mask = npy.where(sigarr < siglimit, True, isonVarBowl.mask)
                #vardiffsgSum.mask = npy.where(sigarr < siglimit, True, vardiffsgSum.mask)

                isonVarBowl = npy.reshape(isonVarBowl,[timN,latN,lonN])
                #vardiffsgSum = npy.reshape(vardiffsgSum,[timN,latN,lonN])

                isonVarBowl = maskVal(isonVarBowl, valmask)
                #vardiffsgSum = maskVal(vardiffsgSum, valmask)
                # Find max of Std dev of all members
                isonVarStd = npy.ma.max(varstd, axis=0)
                # mask
                isonVarStd = maskVal(isonVarStd, valmask)

                # Write
                #isonave = cdm.createVariable(isonVarAve, axes = sigmaTimeList, id = isonRead.id)
                #isonave.long_name = isonRead.long_name
                #isonave.units     = isonRead.units
                #vardiffsgSum = npy.reshape(vardiffsgSum,[delta_ib,timN,latN,lonN])
                #isonavediff = cdm.createVariable(vardiffsgSum, axes = sigmaTimeList, id = isonRead.id+'Agree')
                #isonavediff.long_name = isonRead.long_name
                #isonavediff.units     = isonRead.units
                isonVarBowl = npy.reshape(isonVarBowl,[delta_ib,timN,latN,lonN])
                isonavebowl = cdm.createVariable(isonVarBowl, axes = sigmaTimeList, id = isonRead.id+'Bowl')
                isonavebowl.long_name = isonRead.long_name
                isonavebowl.units     = isonRead.units
                isonVarStd = npy.reshape(isonVarStd,[delta_ib,latN,lonN])
                isonmaxstd = cdm.createVariable(isonVarStd, axes = sigmaList, id = isonRead.id+'Std')
                isonmaxstd.long_name = isonRead.long_name
                isonmaxstd.units     = isonRead.units

                #outFile_f.write(    isonave.astype('float32'), extend = 1, index = ib)
                #outFile_f.write(isonavediff.astype('float32'), extend = 1, index = ib)
                outFile_f.write(isonavebowl.astype('float32'), extend = 1, index = ib)
                outFile_f.write(isonmaxstd.astype('float32'), extend = 1, index = ib)

            tim3 = timc.clock()

            if ToeType == 'histnat':
                isontoe1 = cdm.createVariable(varToE1, axes = [ensembleAxis,axesList[1],axesList[2],axesList[3]], id = isonRead.id+'ToE1')
                isontoe1.long_name = 'ToE 1 for '+isonRead.long_name
                isontoe1.units     = 'Year'
                isontoe2 = cdm.createVariable(varToE2, axes = [ensembleAxis,axesList[1],axesList[2],axesList[3]], id = isonRead.id+'ToE2')
                isontoe2.long_name = 'ToE 2 for '+isonRead.long_name
                isontoe2.units     = 'Year'
                outFile_f.write(isontoe1.astype('float32'), extend = 1, index = ib)
                outFile_f.write(isontoe2.astype('float32'), extend = 1, index = ib)

            tim4 = timc.clock()
        # <--- end of loop on variables

        #print 'ib, timing',ib, tim01-tim0,tim1-tim01,tim2-tim1,tim3-tim2,tim4-tim3
    # <--- end of loop on density
    print ' '

    outFile_f.close()
    fi.close()


def mmeAveMsk1D(listFiles, sw2d, years, inDir, outDir, outFile, timeInt, mme, ToeType, fullTS, debug=True):
    '''
    The mmeAveMsk1D() function averages rhon or scalar density bined files with differing masks
    It ouputs the MME and a percentage of non-masked bins
    
    Created on Tue Nov 25 13:56:20 CET 2014

    Inputs:
    -------
    - listFiles(str)         - the list of files to be averaged
    - sw2d                   - dimension of fields to consider (1 or 2)
    - years(t1,t2)           - years for slice read
    - inDir(str)             - input directory where files are stored
    - outDir(str)            - output directory
    - outFile(str)           - output file
    - timeInt(2xindices)     - indices of init period to compare with (e.g. [1,20])
    - mme(bool)              - multi-model mean (will read in single model ensemble stats)
    - FfllTS                 - 0/1: if 1, uses full time serie (ignores years(t1,t2))
    - debug <optional>       - boolean value

    Notes:
    -----
    - EG 25 Nov 2014   - Initial function write
    - EG  9 Dec 2014   - Add agreement on difference with init period - save as <var>Agree
    - EG 04 Oct 2016   - Add 3D files support

    TODO:
    ------

    '''

    # CDMS initialisation - netCDF compression
    comp = 1 ; # 0 for no compression
    cdm.setNetcdfShuffleFlag(comp)
    cdm.setNetcdfDeflateFlag(comp)
    cdm.setNetcdfDeflateLevelFlag(comp)
    cdm.setAutoBounds('on')
    # Numpy initialisation
    npy.set_printoptions(precision=2)

    if debug:
        debug = True
    else:
        debug = False
    # File dim and grid inits
    t1 = years[0]
    t2 = years[1]
    if t2 <= 0:
        useLastYears = True
        t2 = -t2
    else:
        useLastYears = False
    # Bound of period average to remove
    peri1 = timeInt[0]
    peri2 = timeInt[1]
    # Find dimension
    runN = len(listFiles)
    try:
        fi = cdm.open(inDir[0]+'/'+listFiles[0])
    except:
        print ' *** file not found ',inDir[0]+'/'+listFiles[0]
        sys.exit(' Abort')
    if sw2d == 1:
        ptopd0  = fi['ptopdepth'] ; # Create variable handle
        latN = ptopd0.shape[2]
        basN = ptopd0.shape[1]
    elif sw2d == 2:
        ptopd0  = fi['ptopdepthxy'] ; # Create variable handle
        lonN = ptopd0.shape[2]
        latN = ptopd0.shape[1]

    #timN = ptopd0.shape[0]
    timN = t2-t1
    if fullTS:
        print '  !!! Working on full Time Serie (fullTS = True)'
        timN = ptopd0.shape[0]
        t1=0
        t2=timN
    t10 = t1
    t20 = t2
    # Get grid objects
    axesList = ptopd0.getAxisList()
    # Declare and open files for writing
    if os.path.isfile(outDir+'/'+outFile):
        os.remove(outDir+'/'+outFile)
    outFile_f = cdm.open(outDir+'/'+outFile,'w')

    print ' Number of members:',len(listFiles)

    valmask = ptopd0.missing_value

    # init time axis
    time       = cdm.createAxis(npy.float32(range(timN)))
    time.id    = 'time'
    time.units = 'years since 1861'
    time.designateTime()

    # loop on variables
    # init percent array

    if sw2d == 1:
        varList = ['ptopdepth','ptopsigma','ptopso','ptopthetao','volpers','salpers','tempers']
        #varList = ['ptopdepth']
        varDim  = [1,1,1,1,0,0,0]
        percent  = npy.ma.ones([runN,timN,basN,latN], dtype='float32')*0.
    elif sw2d == 2:
        varList = ['ptopdepthxy','ptopsigmaxy','ptopsoxy','ptopthetaoxy']
        #varList = ['ptopdepthxy']
        varDim  = [2,2,2,2]
        percent  = npy.ma.ones([runN,timN,latN,lonN], dtype='float32')*0.

    varFill = [valmask,valmask,valmask,valmask,valmask,valmask,valmask,valmask,valmask]

    axis1D = [time,axesList[1],axesList[2]]
    axis0D = [time,axesList[1]]
    print ' timN = ',timN

    # loop on 1D variables
    for iv,var in enumerate(varList):
        ti0 = timc.clock()

        # Array inits
        if varDim[iv] == 2:
            isonvar = npy.ma.ones([runN,timN,latN,lonN], dtype='float32')*valmask
            vardiff = npy.ma.ones([runN,timN,latN,lonN], dtype='float32')*valmask
            varones = npy.ma.ones([runN,timN,latN,lonN], dtype='float32')*1.
            axisVar = axis1D
        elif varDim[iv] == 1:
            isonvar = npy.ma.ones([runN,timN,basN,latN], dtype='float32')*valmask
            vardiff = npy.ma.ones([runN,timN,basN,latN], dtype='float32')*valmask
            varones = npy.ma.ones([runN,timN,basN,latN], dtype='float32')*1.
            axisVar = axis1D
        else:
            isonvar = npy.ma.ones([runN,timN,basN], dtype='float32')*valmask
            vardiff = npy.ma.ones([runN,timN,basN], dtype='float32')*valmask
            varones = npy.ma.ones([runN,timN,basN], dtype='float32')*1.
            axisVar = axis0D
        print ' Variable ',iv, var, varDim[iv]
        # loop over files to fill up array
        for ic,file in enumerate(listFiles):
            #print ic,file
            ft      = cdm.open(inDir[0]+'/'+file)
            timeax  = ft.getAxis('time')
            tmax = timeax.shape[0]
            if ic == 0:
                tmax0 = tmax
            #adapt [t1,t2] time bounds to piControl last NN years
            if useLastYears:
                t1 = tmax-t20
                t2 = tmax
            else:
                if tmax != tmax0:
                    print 'wrong time axis: exiting...'
                    return
            #print 'Time dims:',ic, t1,t2,tmax
            # read array
            computeVar = True
            allVars = ft.variables.keys()
            if 'ptopsigmaxy' in allVars:
                computeVar = False
            if (var == 'ptopsigmaxy') & computeVar:
                #print '  ic = ',ic
                # reconstruct from isondepthg and ptopdepthxy

                isond = ft('isondepthg',time = slice(t1,t2))
                #print isond.data.shape, timN*latN*lonN
                itest = 94*360+150
                axesList = isond.getAxisList()
                levs = axesList[1][:]
                levN = len(levs)
                #ti02 = timc.clock()
                levs3d0  = mv.reshape(npy.tile(levs,latN*lonN),(latN*lonN,levN))
                #ti05 = timc.clock()
                isonRead = npy.ma.ones([timN,latN,lonN], dtype='float32')*valmask
                for it in range(timN): # loop on time to limit memory usage
                    levs3d = levs3d0*1.
                    depthlo = mv.reshape(vardepth[ic,it,...],latN*lonN)
                    depth3d = npy.reshape(npy.repeat(depthlo,levN),(latN*lonN,levN))
                    isond3d = mv.reshape(npy.transpose(isond.data[it,...],(1,2,0)),(latN*lonN,levN))
                    #print isond3d[itest,:]
                    isond3d[isond3d > valmask/10] = 0.
                    #print isond3d[itest,:]
                    isond3dp1 = npy.roll(isond3d,-1,axis=1)
                    isond3dp1[:,-1] = isond3d[:,-1]
                    #print isond3dp1[itest,:]
                    #levs3d[levs3d > 30. ] = 0. # to distinguish bottom masked points from surface masked points
                    #print levs3d[itest,:]
                    levs3d[(depth3d <= isond3d)] = 0.
                    #print levs3d[itest,:]
                    levs3d[(depth3d > isond3dp1)] = 0.
                    #print levs3d[itest,:]
                    #isonwrk = npy.sum(levs3d,axis=1)
                    isonwrk = npy.max(levs3d,axis=1)
                    if it < 0:
                        print ic,it
                        print depthlo[itest]
                        print isond3d[itest,:]
                        print isonwrk[itest]
                        print
                    isonRead[it,...] = mv.reshape(isonwrk,(latN,lonN))
                # <-- end of loop on time
                del (isond3d,isond3dp1); gc.collect()
                # mask with depthxy and where sigmaxy = 0
                isonRead.mask = vardepth.mask[ic,...]
                isonRead = mv.masked_where(isonRead == 0, isonRead)
                isonRead.long_name = var
                isonRead.units = 'sigma_n'
                isonRead.id = var
                del (isond,depth3d,levs3d,levs3d0,isonwrk); gc.collect()
                #ti3 = timc.clock()
                #print ti02-ti0,ti05-ti02, ti1-ti05,ti12-ti1,ti15-ti12,ti2-ti15,ti3-ti2
                #print ti3-ti0
                # write ptopsigmaxy
                if os.path.isfile(inDir[0]+'/work_ptopsigmaxy/'+file):
                    os.remove(inDir[0]+'/work_ptopsigmaxy/'+file)
                fiout = cdm.open(inDir[0]+'/work_ptopsigmaxy/'+file,'w')
                if ic == 0:
                        print ' Creating ',inDir[0]+'/work_ptopsigmaxy/'+file
                isonsigxy = cdm.createVariable(isonRead, axes = axis1D, id = 'ptopsigmaxy')
                isonsigxy.long_name = 'Density of shallowest persistent ocean on ison'
                isonsigxy.units     = 'sigma_n'
                fiout.write(isonsigxy.astype('float32'))
                fiout.close()
            else:
                # Direct read of variable
                isonRead = ft(var,time = slice(t1,t2))
            #print isonRead.shape, timN
            if varFill[iv] != valmask:
                isonvar[ic,...] = isonRead.filled(varFill[iv])
            else:
                isonvar[ic,...] = isonRead
            #print isonvar[ic,:,40,100]
            # compute percentage of non-masked points accros MME
            if iv == 0:
                maskvar = mv.masked_values(isonRead.data,valmask).mask
                percent[ic,...] = npy.float32(npy.equal(maskvar,0))
            if mme:
                # if mme then just average Bowl and Agree fields
                varst = var+'Agree'
                vardiff[ic,...] = ft(varst,time = slice(t1,t2))
            else:
                # Compute difference with average of first initN years, use mask of last month
                varinit = cdu.averager(isonvar[ic,peri1:peri2,...],axis=0)
                for tr in range(timN):
                    vardiff[ic,tr,...] = isonvar[ic,tr,...] - varinit
                vardiff[ic,...].mask = isonvar[ic,...].mask

            ft.close()
        # <-- end of loop on files
        # TODO remove masked points at longitudes 0 or 180deg for some models
        # if ptopdepthxy, keep for ptopsigmaxy computation (reconstruct from isondepthg and ptopdepthxy)
        if var =='ptopdepthxy':
            vardepth = isonvar
        # Compute percentage of bin presence
        # Only keep points where percent > 50%
        if iv == 0:
            percenta = (cdu.averager(percent,axis=0))*100.
            percenta = mv.masked_less(percenta, 50)
            percentw = cdm.createVariable(percenta, axes = axis1D, id = 'ptoppercent')
            percentw._FillValue = valmask
            percentw.long_name = 'percentage of MME bin'
            percentw.units     = '%'
            outFile_f.write(percentw.astype('float32'))
        # Sign of difference
        if mme:
            vardiffsgSum = cdu.averager(vardiff, axis=0)
            vardiffsgSum = cdm.createVariable(vardiffsgSum , axes = axisVar , id = 'foo')
            vardiffsgSum = maskVal(vardiffsgSum, valmask)
            vardiffsgSum.mask = percentw.mask
        else:
            vardiffsg = npy.copysign(varones,vardiff)
            # average signs
            vardiffsgSum = cdu.averager(vardiffsg, axis=0)
            vardiffsgSum = mv.masked_greater(vardiffsgSum, 10000.)
            vardiffsgSum.mask = percentw.mask
            vardiffsgSum._FillValue = valmask

        # average accross members
        isonVarAve = cdu.averager(isonvar, axis=0)
        isonVarAve = cdm.createVariable(isonVarAve , axes = axisVar , id = 'foo')
        # mask
        if varFill[iv] == valmask:
            isonVarAve = maskVal(isonVarAve, valmask)

        isonVarAve.mask = percentw.mask

        # Write
        isonave = cdm.createVariable(isonVarAve, axes = axisVar, id = isonRead.id)
        isonave.long_name = isonRead.long_name
        isonave.units     = isonRead.units
        isonavediff = cdm.createVariable(vardiffsgSum, axes = axisVar, id = isonRead.id+'Agree')
        isonavediff.long_name = isonRead.long_name
        isonavediff.units     = isonRead.units

        outFile_f.write(isonave.astype('float32'))
        outFile_f.write(isonavediff.astype('float32'))
        tf = timc.clock()
        print '   time var',tf-ti0
    # <--- end of loop on variables 

    outFile_f.close()
    fi.close()

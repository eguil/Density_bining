# -*- coding: utf-8 -*-

# Import
import os, glob, sys
import numpy as np
from netCDF4 import Dataset as open_ncfile
import matplotlib.pyplot as plt
# Add path to PYTHONPATH
sys.path.append("/home/ysilvy/Density_bining/Yona_analysis/programs")
sys.path.append("/home/ysilvy/Density_bining/")
from binDensity import eosNeutral, maskVal
from scipy.interpolate import griddata
import ESMF as ESMP
import cdms2 as cdm
from cdms2 import CdmsRegrid
import MV2 as mv
import time as timc

def maptogamma(fieldz,gammaz,targetrho):
    """
    From field(depth) and gamma(depth), build field(gamma)
    
    input: fieldz - (basin,depth,latitude)
           gammaz - reference gamma (basin,depth,lat)
           targetrho - target 1D density grid for mapping
           
    output: fieldgamma - (basin,gamma,latitude)
    """
    basinN = fieldz.shape[0]
    depthN = fieldz.shape[1]
    latN = fieldz.shape[2]
    
    fieldgamma = np.ma.masked_all((basinN,len(targetrho),latN))
    
    for ibasin in range(basinN):
        for ilat in range(latN):

            gamz = gammaz[ibasin,:,ilat] # Read gamma (function of depth) of water column
            field_z = fieldz[ibasin,:,ilat] # Read field values of the water column

            field_sort = np.ma.compressed(field_z) # Remove masked values
            gam_sort = gamz[np.where(field_z!=np.ma.masked)]

            if len(gam_sort) > 1:
                fieldgamma[ibasin,:,ilat] = griddata(gam_sort,field_sort,targetrho) # Grid field with target pressure grid
            else :
                fieldgamma[ibasin,:,ilat] = np.ma.masked
                
    # Mask nans
    fieldgamma[np.isnan(fieldgamma)] = np.ma.masked

    # Mask out of bounds data
    min = gamz[0]
    idx = np.argmin(np.abs(targetrho-min))
    fieldgamma[ibasin,0:idx-1,ilat] = np.ma.masked

    return fieldgamma

def readzfile(fileT,fileS,compute_gamma,targetGridFile,outdir):
    
    """ Reads a netcdf file in so and in thetao, computes gamma neutral at each point and each time step, makes annual means, interpolates horizonally to rregular grid, makes zonal means and writes to new file
    inputs:
    - fileT: netcdf file for thetao variable (function of time, depth, latitude, longitude)
    - fileS: netcdf file for so variable, same dimensions as fileT
    - compute_gamma: boolean (e.g. don't compute gamma for historicalNat)
    - targetGridFile: netcdf file with target horizontal grid for interpolation + basinmask
    - outdir: directory to write new file in
    """
    
    # == Prepare work ==
    
    t0 = timc.time()
    # Open files to read
    ft = cdm.open(fileT)
    fs = cdm.open(fileS)
    thetao_h = ft('thetao', time = slice(1,2))
    so_h = fs('so', time = slice(1,10))
    print('Opening files :\n - '+os.path.basename(fileT)+'\n - '+os.path.basename(fileS))

    # Horizontal grid
    ingrid  = thetao_h.getGrid()
    # Get grid objects
    axesList = thetao_h.getAxisList()

    depth=thetao_h.getLevel()
    # Define dimensions
    lonN    = thetao_h.shape[3]
    latN    = thetao_h.shape[2]
    depthN  = thetao_h.shape[1]
    timeax  = ft.getAxis('time')
    t1 = timc.time()
    #print(t1-t0)

    thetaoLongName = thetao_h.long_name
    soLongName = so_h.long_name
    soUnits = so_h.units

    # Target horizonal grid for interp
    gridFile_f  = cdm.open(targetGridFile)
    maskg       = gridFile_f('basinmask3')
    outgrid     = maskg.getGrid()
    maski       = maskg.mask ; # Global mask
    # Regional masks
    maskAtl = maski*1 ; maskAtl[...] = True
    idxa = np.argwhere(maskg == 1).transpose()
    maskAtl[idxa[0],idxa[1]] = False
    maskPac = maski*1 ; maskPac[...] = True
    idxp = np.argwhere(maskg == 2).transpose()
    maskPac[idxp[0],idxp[1]] = False
    maskInd = maski*1 ; maskInd[...] = True
    idxi = np.argwhere(maskg == 3).transpose()
    maskInd[idxi[0],idxi[1]] = False
    loni    = maskg.getLongitude()
    lati    = maskg.getLatitude()
    Nii     = len(loni)
    Nji     = len(lati)

    gridFile_f.close()

    t2 = timc.time()
    #print(t2-t1)
    
    # Outfile
    fname = os.path.basename(fileS)
    split = fname.split('_')
    split[0]='so_thetao_gamma'
    split[1]='Oan'
    fname_yearly='_'.join(split)
    outFile = outdir+fname_yearly
    outFile_f = cdm.open(outFile,'w')
    
    # Valmask
    valmask  = 1.e20
    
    # Define time chunks
    tmin = 0
    tmax = timeax.shape[0] # Dimension of total time axis of infiles (total number of months)
    nyrtc = 5 # nb of years per time chunk
    tcdel = 5*12 # nb of months per time chunk
    tcmax = (tmax-tmin)/tcdel ; # number of time chunks
    
    # Initialize interpolated arrays on target grid - global and basin zonal means
    so_interp = np.ma.masked_all((nyrtc, depthN, Nji, Nii),dtype='float32')
    thetao_interp = np.ma.masked_all((nyrtc, depthN, Nji, Nii),dtype='float32')
    soa_interp, sop_interp, soi_interp, thetaoa_interp, thetaop_interp, thetaoi_interp = [np.ma.masked_all((nyrtc, depthN, Nji, Nii),dtype='float32') for _ in range(6)]
    if compute_gamma:
        rhon_interp = np.ma.masked_all((nyrtc, depthN, Nji, Nii),dtype='float32')
        rhona_interp, rhonp_interp, rhoni_interp = [np.ma.masked_all((nyrtc, depthN, Nji, Nii),dtype='float32') for _ in range(3)]

    # Interpolation init (regrid)
    t=timc.time()
    ESMP.ESMP_Initialize()
    regridObj = CdmsRegrid(ingrid,outgrid,so_interp.dtype,missing=valmask,regridMethod='distwgt',regridTool='esmf', coordSys='deg', diag = {},periodicity=1)
    #print(timc.time()-t)
    
    # Define basin output axis
    basinAxis               = cdm.createAxis([0,1,2,3],bounds=None,id='basin')
    basinAxis.long_name     = 'ocean basin index'
    basinAxis.standard_name = 'basin'
    basinAxis.units         = 'basin index'
    basinAxis.units_long    = '0: global_ocean 1: atlantic_ocean; 2: pacific_ocean; 3: indian_ocean'
    basinAxis.axis          = 'B'

    # Create basin-zonal axes lists
    basinTimeList           = [axesList[0],basinAxis] ; # time, basin
    basinAxesList           = [axesList[0],basinAxis,axesList[2]] ; # time, basin, lat
    basinZAxesList        = [axesList[0],basinAxis,axesList[1],axesList[2]] ; # time, basin, depth, lat
    
    
    # ====
    # == Start loop on time chunks ==
    # ====
    
    for tc in range(tcmax): # Loop on time chunks
    #tc=0
        print('- time chunk: '+str(tc))
        # read tcdel month by tcdel month to optimise memory
        trmin   = tmin + tc*tcdel ; # define as function of tc and tcdel
        trmax   = tmin + (tc+1)*tcdel ; # define as function of tc and tcdel
        print('  months: '+str(trmin)+' '+str(trmax))
        thetao  = ft('thetao', time = slice(trmin,trmax))-273.15
        so      = fs('so'    , time = slice(trmin,trmax))

        if compute_gamma:
            # Compute neutral density
            tbfrho = timc.time()
            rhon = eosNeutral(thetao,so)-1000.
            print('  Neutral density computed in '+str(timc.time()-tbfrho))
            
        # Turn into cdms variable
        time    = thetao.getTime()
        newAxesList = [axesList[0],axesList[1],axesList[2],axesList[3]]
        newAxesList[0]  = time ; # replace time axis
        So = cdm.createVariable(so, axes=newAxesList)
        Thetao = cdm.createVariable(thetao, axes=newAxesList)
        if compute_gamma:
            Rhon = cdm.createVariable(rhon, axes=newAxesList)

        # == Compute annual mean ==
        so_temp = np.ma.reshape(So,(nyrtc,12,depthN,latN,lonN))
        thetao_temp = np.ma.reshape(Thetao,(nyrtc,12,depthN,latN,lonN))
        so_yearly = np.ma.average(so_temp,axis=1)
        thetao_yearly = np.ma.average(thetao_temp,axis=1)
        if compute_gamma:
            rhon_temp = np.ma.reshape(Rhon,(nyrtc,12,depthN,latN,lonN))
            rhon_yearly = np.ma.average(rhon_temp,axis=1)

        # Turn into cdms variable for horizontal interpolation
        so_yearly = cdm.createVariable(so_yearly)
        thetao_yearly = cdm.createVariable(thetao_yearly)
        if compute_gamma:    
            rhon_yearly = cdm.createVariable(rhon_yearly)

        # Create annual time axis
        timeyr = cdm.createAxis(np.arange(trmin/12,trmax/12),bounds=None,id='time')
        timeyr.units    = 'years'
        timeyr.designateTime()
        newAxesList[0] = timeyr # replace time axis

        so_yearly.setAxisList(newAxesList)
        thetao_yearly.setAxisList(newAxesList)
        if compute_gamma:
            rhon_yearly.setAxisList(newAxesList)
        
        # == Interpolate onto regular grid ==
        for t in range(nyrtc):
            for ks in range(depthN):
                # Global
                so_interp[t,ks,:,:]         = regridObj(so_yearly[t,ks,:,:])
                so_interp[t,ks,:,:].mask    = maski
                thetao_interp[t,ks,:,:]     = regridObj(thetao_yearly[t,ks,:,:])
                thetao_interp[t,ks,:,:].mask= maski
                if compute_gamma:
                    rhon_interp[t,ks,:,:]     = regridObj(rhon_yearly[t,ks,:,:])
                    rhon_interp[t,ks,:,:].mask= maski

                # Atl
                soa_interp[t,ks,:,:]          = so_interp[t,ks,:,:]*1.
                soa_interp[t,ks,:,:].mask     = maskAtl
                thetaoa_interp[t,ks,:,:]      = thetao_interp[t,ks,:,:]*1.
                thetaoa_interp[t,ks,:,:].mask = maskAtl
                if compute_gamma:
                    rhona_interp[t,ks,:,:]      = rhon_interp[t,ks,:,:]*1.
                    rhona_interp[t,ks,:,:].mask = maskAtl

                # Pac
                sop_interp[t,ks,:,:]          = so_interp[t,ks,:,:]*1.
                sop_interp[t,ks,:,:].mask     = maskPac
                thetaop_interp[t,ks,:,:]      = thetao_interp[t,ks,:,:]*1.
                thetaop_interp[t,ks,:,:].mask = maskPac
                if compute_gamma:
                    rhonp_interp[t,ks,:,:]      = rhon_interp[t,ks,:,:]*1.
                    rhonp_interp[t,ks,:,:].mask = maskPac

                # Ind
                soi_interp[t,ks,:,:]          = so_interp[t,ks,:,:]*1.
                soi_interp[t,ks,:,:].mask     = maskInd
                thetaoi_interp[t,ks,:,:]      = thetao_interp[t,ks,:,:]*1.
                thetaoi_interp[t,ks,:,:].mask = maskInd
                if compute_gamma:
                    rhoni_interp[t,ks,:,:]      = rhon_interp[t,ks,:,:]*1.
                    rhoni_interp[t,ks,:,:].mask = maskInd
        
        # Mask after interpolation
        so_interp = maskVal(so_interp,valmask)
        soa_interp = maskVal(soa_interp,valmask)
        sop_interp = maskVal(sop_interp,valmask)
        soi_interp = maskVal(soi_interp,valmask)
        thetao_interp = maskVal(thetao_interp,valmask)
        thetaoa_interp = maskVal(thetaoa_interp,valmask)
        thetaop_interp = maskVal(thetaop_interp,valmask)
        thetaoi_interp = maskVal(thetaoi_interp,valmask)
        if compute_gamma:
            rhon_interp = maskVal(rhon_interp,valmask)
            rhona_interp = maskVal(rhona_interp,valmask)
            rhonp_interp = maskVal(rhonp_interp,valmask)
            rhoni_interp = maskVal(rhoni_interp,valmask)
            
        # == Compute zonal means ==
        soz = np.ma.average(so_interp, axis=3) # Global
        soza = np.ma.average(soa_interp, axis=3) # Atlantic
        sozp = np.ma.average(sop_interp, axis=3) # Pacific
        sozi = np.ma.average(soi_interp, axis=3) # Indian
        thetaoz = np.ma.average(thetao_interp, axis=3)
        thetaoza = np.ma.average(thetaoa_interp, axis=3)
        thetaozp = np.ma.average(thetaop_interp, axis=3)
        thetaozi = np.ma.average(thetaoi_interp, axis=3)
        if compute_gamma:
            rhonz = np.ma.average(rhon_interp, axis=3)
            rhonza = np.ma.average(rhona_interp, axis=3)
            rhonzp = np.ma.average(rhonp_interp, axis=3)
            rhonzi = np.ma.average(rhoni_interp, axis=3)
        
        # == Write zonal means to outfile ==
        # Prepare axis
        timeBasinZAxesList = basinZAxesList
        timeBasinZAxesList[0] = timeyr ; # Replace monthly with annual
        timeBasinZAxesList[3] = lati ; # Replace lat with regrid target

        # Collapse onto basin axis
        Sz = np.ma.stack([soz,soza,sozp,sozi],axis=1)
        del(soz,soza,sozp,sozi)
        Sz = cdm.createVariable(Sz,axes=timeBasinZAxesList,id='salinity')
        Tz = np.ma.stack([thetaoz,thetaoza,thetaozp,thetaozi],axis=1)
        del(thetaoz,thetaoza,thetaozp,thetaozi)
        Tz = cdm.createVariable(Tz,axes=timeBasinZAxesList,id='temperature')
        if compute_gamma:
            Rz = np.ma.stack([rhonz,rhonza,rhonzp,rhonzi],axis=1)
            del(rhonz,rhonza,rhonzp,rhonzi)
            Rz = cdm.createVariable(Rz,axes=timeBasinZAxesList,id='density')

        if tc == 0:
            # Global attributes
            Sz.long_name  = soLongName
            Sz.units = soUnits
            Tz.long_name  = thetaoLongName
            Tz.units      = 'degrees_C'
            if compute_gamma:
                Rz.long_name  = 'Neutral density'
                Rz.units      = 'kg.m-3'

        # Write & append
        outFile_f.write(Sz.astype('float32'), extend = 1, index = (trmin-tmin)/12)
        outFile_f.write(Tz.astype('float32'), extend = 1, index = (trmin-tmin)/12)
        del(Sz,Tz) 
        if compute_gamma:
            outFile_f.write(Rz.astype('float32'), extend = 1, index = (trmin-tmin)/12)
            del(Rz)

        outFile_f.sync()
        
    # == End of loop on time chunks ==
    ft.close()
    fs.close()
    outFile_f.close()
    tf = timc.time()-t0
    print('Wrote file: '+outFile)
    print('Total time :'+str(tf)+'s ('+str(tf/60)+'mn)')

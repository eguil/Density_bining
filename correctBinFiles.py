import os,sys,gc,glob
import cdms2 as cdm
import cdutil as cdu
import MV2 as mv
import numpy as npy
from string import replace

def correctFile(model,inFile, inDir, outFile, outDir):
    '''
    Correct density binned files (undefined ptop & long 0 issue)
    CCSM4: i=139-140, j=0,145
    '''
    # CDMS initialisation - netCDF compression
    comp = 1 # 0 for no compression
    cdm.setNetcdfShuffleFlag(comp)
    cdm.setNetcdfDeflateFlag(comp)
    cdm.setNetcdfDeflateLevelFlag(comp)
    cdm.setAutoBounds('on')
    # Numpy initialisation
    npy.set_printoptions(precision=2)

    varList3D = ['isondepthg','isonthickg', 'sog','thetaog']
    varList2D = ['persistmxy','ptopdepthxy','ptopsigmaxy','ptopsoxy','ptopthetaoxy',]

    # First test, read level by level and write level by level (memory management)
    # use ncpdq -a time,lev,lat,lon to recover the dimension order

    fi = cdm.open(inDir+'/'+inFile)
    fo = cdm.open(outDir+'/'+outFile,'w')
    isondg  = fi['isondepthg'] ; # Create variable handle
        # Get grid objects
    axesList = isondg.getAxisList()
    sigmaGrd = isondg.getLevel()
    lonN = isondg.shape[3]
    latN = isondg.shape[2]
    levN = isondg.shape[1]
    timN = isondg.shape[0]
    valmask = isondg.missing_value

    if model == 'CCSM4':
        ic1 = 139
        ic2 = 140
        jcmax = 145

    for it in range(1):
        print it
        # test
        i = 90
        j = 90
        ij = j*360+i
        #print 'ij=',ij
        # 3D variables
        for iv in varList3D:
            print iv
            var = fi(iv,time = slice(it,it+1))
            outVar = var
            # Correct
            for jt in range(jcmax):
                outVar[:,:,jt,ic1] = (outVar[:,:,jt,ic1-1]+outVar[:,:,jt,ic2+1])/2
                outVar[:,:,jt,ic2] = outVar[:,:,jt,ic1]
            if iv == 'sog':
                varsog = npy.reshape(outVar,(levN,latN*lonN))
            elif iv =='thetaog':
                varthetao = npy.reshape(outVar,(levN,latN*lonN))
            elif iv =='isondepthg':
                # find indices of surface points
                vardepth = npy.min(npy.reshape(outVar,(levN,latN*lonN)),axis=0)
                print vardepth.shape
                #print idxdepth[:,ij]
                #idxdepth1 = npy.min(vardepth,axis=0).transpose()
                #idxrep = npy.reshape(npy.tile(idxdepth1,levN),(levN,latN*lonN))
                #idxsurf = npy.argwhere(idxrep == vardepth).transpose()
                #print idxrep.shape, idxrep[0,ij]
                #print idxsurf.shape
                #for te in range(idxsurf.shape[1]):
                #    if idxsurf[1,te] == ij:
                #        print te, idxsurf[0,te]

            # Write
            outVar.long_name = var.long_name
            outVar.units = var.units
            fo.write(outVar.astype('float32'), extend = 1, index = it)
            fo.sync()
        # 2D variables and correct isondepthg = 0
        for iv in varList2D:
            outVar = fi(iv,time = slice(it,it+1))
            # Correct for longitude interpolation issue
            for jt in range(jcmax):
                outVar[:,jt,ic1] = (outVar[:,jt,ic1-1]+outVar[:,jt,ic2+1])/2
                outVar[:,jt,ic2] = outVar[:,jt,ic1]
            # Correct for ptopdepthxy = 0
            if iv == 'ptopdepthxy':
                outVar.data[...] = npy.where(npy.reshape(outVar,(latN*lonN)) == 0.,vardepth,npy.reshape(outVar,(latN*lonN))).reshape(outVar.shape)[...]

            #if iv == 'ptopsoxy':


            #elif iv == 'ptopthetaoxy':
            #elif iv == 'ptopsigmaxy':

            # Write
            outVar.long_name = var.long_name
            outVar.units = var.units
            fo.write(outVar.astype('float32'), extend = 1, index = it)
            fo.sync()

    fi.close()
    fo.close()

# testing
model = 'CCSM4'
inFile = 'cmip5.CCSM4.historical24.r1i1p1.an.ocn.Omon.density.ver-v20121128.nc'
inDir = '/Users/ericg/Projets/Density_bining/Raw_testing'
outFile = 'cmip5.CCSM4.historical24.outtest.nc'
outDir = inDir

correctFile(model,inFile, inDir, outFile, outDir)
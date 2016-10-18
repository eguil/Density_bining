import gc
import cdms2 as cdm
import MV2 as mv
import numpy as npy

def correctFile(idxcorr, ncorr, inFile, inDir, outFile, outDir):
    '''
    Correct density binned files (undefined ptop & long 0 issue)
    idxcorr = [idx_i,idx_i1,jmax] indices for longitude correction - if [0,0,0] ignore
    ncorr   = number of corrections
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
    varList2D = ['ptopsoxy','ptopdepthxy','ptopsigmaxy','ptopthetaoxy','persistmxy']

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

    if ncorr == 2:
        ic1 = idxcorr[0][0]
        ic2 = idxcorr[0][1]
        jcmax = idxcorr[0][2]
        ic12 = idxcorr[1][0]
        ic22 = idxcorr[1][1]
        jcmax2 = idxcorr[1][2]
        if ic2 >= lonN-1:
            ic2 = 0
        if ic22 >= lonN-1:
            ic22 = 0
    elif ncorr == 1:
        ic1 = idxcorr[0]
        ic2 = idxcorr[1]
        jcmax = idxcorr[2]
        if ic2 >= lonN-1:
            ic2 = 0
    #print ic1,ic2,jcmax
    corr_long = True
    if ic1 == 0 and ic2 == 0 and jcmax == 0:
        corr_long = False
    testp = 10
    for it in range(timN):
        #if it/testp*testp == it:
        #    print ' year =',it
        # test
        i = 90
        j = 90
        i2d = 6
        j2d = 12
        ij = j*lonN+i
        ij2d = j2d*lonN+i2d
        #print 'ij=',ij
        # 3D variables
        for iv in varList3D:
            #print iv
            outVar = fi(iv,time = slice(it,it+1))
            # Correct for longitude interpolation issue
            if corr_long:
                for jt in range(jcmax):
                    outVar[:,:,jt,ic1] = (outVar[:,:,jt,ic1-1]+outVar[:,:,jt,ic2+1])/2
                    outVar[:,:,jt,ic2] = outVar[:,:,jt,ic1]
                if ncorr == 2:
                    for jt in range(jcmax2):
                        outVar[:,:,jt,ic12] = (outVar[:,:,jt,ic12-1]+outVar[:,:,jt,ic22+1])/2
                        outVar[:,:,jt,ic22] = outVar[:,:,jt,ic12]
            # Correct Bowl properties
            if iv =='isondepthg':
                vardepth = npy.reshape(outVar,(levN,latN*lonN))
                #print 'test'
                #print outVar[:,:,j2d,i2d]
                #print vardepth[:,ij2d]
                # find values of surface points
                vardepthBowl = npy.min(npy.reshape(outVar,(levN,latN*lonN)),axis=0)
                vardepthBowlTile = npy.repeat(vardepthBowl,levN,axis=0).reshape((latN*lonN,levN)).transpose()
                #print vardepthBowlTile.shape
                #print vardepthBowl[ij2d], vardepthBowlTile[:,ij2d]
                levs = outVar.getAxisList()[1][:]
                #print 'levs',levs
                levs3d  = mv.reshape(npy.tile(levs,latN*lonN),(latN*lonN,levN)).transpose()
                varsigmaBowl = npy.max(npy.where(vardepth == vardepthBowlTile,levs3d,0),axis=0)
                #print varsigmaBowl[ij2d],levs3d[:,ij2d]

            elif iv == 'sog':
                varsog = npy.reshape(outVar,(levN,latN*lonN))
                varsoBowl = npy.max(npy.where(vardepth == vardepthBowlTile,varsog,0),axis=0)
                #print varsoBowl[ij2d], varsog[:,ij2d]
                #print vardepth[:,ij2d],vardepthBowlTile[:,ij2d]
                del (varsog); gc.collect()
            elif iv =='thetaog':
                varthetao = npy.reshape(outVar,(levN,latN*lonN))
                varthetaoBowl = npy.max(npy.where(vardepth == vardepthBowlTile,varthetao,-1000),axis=0)
                #print varthetaoBowl[ij2d],varthetao[:,ij2d]
                del (varthetao); gc.collect()
            # Write
            fo.write(outVar.astype('float32'), extend = 1, index = it)
            fo.sync()
        del (vardepth); gc.collect()
        # 2D variables and correct isondepthg = 0
        for iv in varList2D:
            outVar = fi(iv,time = slice(it,it+1))
            # Correct for longitude interpolation issue
            if corr_long:
                for jt in range(jcmax):
                    outVar[:,jt,ic1] = (outVar[:,jt,ic1-1]+outVar[:,jt,ic2+1])/2
                    outVar[:,jt,ic2] = outVar[:,jt,ic1]
                if ncorr == 2:
                    for jt in range(jcmax2):
                        outVar[:,jt,ic12] = (outVar[:,jt,ic12-1]+outVar[:,jt,ic22+1])/2
                        outVar[:,jt,ic22] = outVar[:,jt,ic12]
            # Correct for ptopsoxy < 30
            #print 'before',outVar[:,j2d,i2d]
            if iv == 'ptopsoxy':
                testso = npy.reshape(outVar,(latN*lonN)) < 30.
                #print 'testdepth', testdepth[ij2d]
                #print npy.argwhere(testdepth)[0:10]/lonN, npy.argwhere(testdepth)[0:10]-npy.argwhere(testdepth)[0:10]/lonN*lonN
                outVar.data[...] = npy.where(testso,varsoBowl,npy.reshape(outVar,(latN*lonN))).reshape(outVar.shape)[...]
            elif iv == 'ptopdepthxy':
                outVar.data[...] = npy.where(testso,vardepthBowl,npy.reshape(outVar,(latN*lonN))).reshape(outVar.shape)[...]
            elif iv == 'ptopthetaoxy':
                outVar.data[...] = npy.where(testso,varthetaoBowl,npy.reshape(outVar,(latN*lonN))).reshape(outVar.shape)[...]
            elif iv == 'ptopsigmaxy':
                outVar.data[...] = npy.where(testso,varsigmaBowl,npy.reshape(outVar,(latN*lonN))).reshape(outVar.shape)[...]
            #print 'after',outVar[:,j2d,i2d]

            # Write
            fo.write(outVar.astype('float32'), extend = 1, index = it)
            fo.sync()

    fi.close()
    fo.close()

# testing
#model = 'CCSM4'
#idxcorr=[139,140,145]
#ncorr = 1
#inFile = 'cmip5.CCSM4.historical24.r1i1p1.an.ocn.Omon.density.ver-v20121128.nc'
#inDir = '/Users/ericg/Projets/Density_bining/Raw_testing'
#outFile = 'cmip5.CCSM4.historical24.outtest.nc'

#model = 'CanESM2'
#idxcorr=[179,180,180]
#ncorr = 1
#inFile = 'cmip5.CanESM2.historical24.r1i1p1.an.ocn.Omon.density.ver-1.nc'
#inDir = '/Users/ericg/Projets/Density_bining/Raw_testing'
#outFile = 'cmip5.CanESM2.historical24.outtest.nc'

#model = 'IPSL-CM5A-LR'
#idxcorr=[0,0,0]
#ncorr=1
#inFile = 'cmip5.IPSL-CM5A-LR.historical24.r1i1p1.an.ocn.Omon.density.ver-v20111119.nc'
#inDir = '/Users/ericg/Projets/Density_bining/Raw_testing'
#outFile = 'cmip5.IPSL-CM5A-LR.historical24.outtest.nc'


#model = 'Ishii'
#idxcorr=[[359,359,39],[180,180,180]]
#idxcorr=[359,359,39]
#ncorr = 1
#inFile = 'obs.Ishii.historical.r0i0p0.an.ocn.Omon.density.ver-1.latestX.nc'
#inDir='/Volumes/hciclad/data/Density_binning/Prod_density_obs_april16'
#outFile = 'obs.Ishii.historical.r0i0p0.an.ocn.Omon.density.ver-1.latestXCorr.nc'

#model = 'EN4'
#idxcorr=[[359,359,39],[180,180,180]]
#idxcorr=[359,359,39]
#ncorr = 2
#inFile = 'obs.EN4.historical.r0i0p0.mo.ocn.Omon.density.ver-1.latestX.nc'
#inDir='/Volumes/hciclad/data/Density_binning/Prod_density_obs_april16'
#outFile = 'obs.EN4.historical.r0i0p0.mo.ocn.Omon.density.ver-1.latestXCorr.nc'

#outDir = inDir


#correctFile(idxcorr, ncorr, inFile, inDir, outFile, outDir)
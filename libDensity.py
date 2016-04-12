import os,sys
import cdms2 as cdm
import cdutil as cdu
import MV2 as mv
import numpy as npy
from string import replace
from genutil import statistics



def defModels():
#
# List available models and properties
#
#  name, props=[Nb of hist members, nb of HistNat members, idx of common interval [1861-2005]]
#      picontrol=[length of run]
#
#  set last idx (2005) to 0 to ignore model
#
    models = [
        {'name':'ACCESS1-0'     ,'props':[2,0,11,156], 'picontrol':[0]}, # 0
        {'name':'ACCESS1-3'     ,'props':[3,0,11,156], 'picontrol':[0]}, # 1
        {'name':'bcc-csm1-1-m'  ,'props':[3,0,11,156], 'picontrol':[0]}, # 2
        {'name':'bcc-csm1-1'    ,'props':[2,1,11,156], 'picontrol':[0]}, # 3
        {'name':'BNU-ESM'       ,'props':[1,0,11,156], 'picontrol':[0]}, # 4
        {'name':'CanESM2'       ,'props':[5,5,11,156], 'picontrol':[0]}, # 5
        {'name':'CCSM4'         ,'props':[6,4,11,156], 'picontrol':[0]}, # 6
        {'name':'CESM1-BGC'     ,'props':[1,0,11,156], 'picontrol':[0]}, # 7
        {'name':'CESM1-CAM5'    ,'props':[3,2,11,156], 'picontrol':[0]}, # 8
        {'name':'CESM1-FASTCHEM','props':[3,0,11,156], 'picontrol':[0]}, # 9
        {'name':'CESM1-WACCM'   ,'props':[1,0,11,156], 'picontrol':[0]}, # 10
        {'name':'CMCC-CESM'     ,'props':[1,0,11,156], 'picontrol':[0]}, # 11
        {'name':'CMCC-CM'       ,'props':[1,0,11,156], 'picontrol':[0]}, # 12
        {'name':'CMCC-CMS'      ,'props':[1,0,11,156], 'picontrol':[0]}, # 13
        {'name':'CNRM-CM5-2'    ,'props':[1,0,11,156], 'picontrol':[0]}, # 14
        {'name':'CNRM-CM5'      ,'props':[9,6,11,156], 'picontrol':[0]}, # 15
        {'name':'CSIRO-Mk3-6-0' ,'props':[9,5,11,156], 'picontrol':[0]}, # 16
        {'name':'CSIRO-Mk3L-1-2','props':[2,0,10,155], 'picontrol':[0]}, # 17 - 145 time steps
        {'name':'EC-EARTH'      ,'props':[6,0,11,156], 'picontrol':[0]}, # 18
        {'name':'FGOALS-g2'     ,'props':[4,3,11,156], 'picontrol':[0]}, # 19
        {'name':'FGOALS-s2'     ,'props':[3,0,11,156], 'picontrol':[0]}, # 20
        {'name':'GFDL-CM2p1'    ,'props':[9,0,11,156], 'picontrol':[0]}, # 21
        {'name':'GFDL-CM3'      ,'props':[4,3, 0,146], 'picontrol':[0]}, # 22 - 145 time steps
        {'name':'GFDL-ESM2G'    ,'props':[1,0, 0,146], 'picontrol':[0]}, # 23 - 145 time steps
        {'name':'GFDL-ESM2M'    ,'props':[1,1, 0,146], 'picontrol':[0]}, # 24 - 145 time steps
        {'name':'GISS-E2-H'     ,'props':[14,11,11,156],'picontrol':[0]},# 25
        {'name':'GISS-E2-H-CC'  ,'props':[1,0,11,156], 'picontrol':[0]}, # 26
        {'name':'GISS-E2-R'     ,'props':[16,11,11,156], 'picontrol':[0]}, # 27
        {'name':'GISS-E2-R-CC'  ,'props':[1,0,11,156], 'picontrol':[0]}, # 28
        {'name':'HadCM3'        ,'props':[9,0, 1,146], 'picontrol':[0]}, # 29
        {'name':'HadGEM2-CC'    ,'props':[1,0, 1,146], 'picontrol':[0]}, # 30
        {'name':'HadGEM2-ES'    ,'props':[3,3, 1,146], 'picontrol':[0]}, # 31
        {'name':'IPSL-CM5A-LR'  ,'props':[6,3,11,156], 'picontrol':[0]}, # 32
        {'name':'IPSL-CM5A-MR'  ,'props':[3,3,11,156], 'picontrol':[0]}, # 33
        {'name':'IPSL-CM5B-LR'  ,'props':[1,0,11,156], 'picontrol':[0]}, # 34
        {'name':'MIROC-ESM-CHEM','props':[1,1,11,156], 'picontrol':[0]}, # 35
        {'name':'MIROC-ESM'     ,'props':[3,3,11,156], 'picontrol':[0]}, # 36
        {'name':'MPI-ESM-LR'    ,'props':[3,0,11,156], 'picontrol':[0]}, # 37
        {'name':'MPI-ESM-MR'    ,'props':[3,0,11,156], 'picontrol':[0]}, # 38
        {'name':'MPI-ESM-P'     ,'props':[2,0,11,156], 'picontrol':[0]}, # 39
        {'name':'NorESM1-ME'    ,'props':[1,0,11,156], 'picontrol':[0]}, # 40
        {'name':'NorESM1-M'     ,'props':[2,0,11,156], 'picontrol':[0]}, # 41
        ]

    return models

# models to check
"""
    models = [
        ]
"""

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

def findToE(signal, noise, mult):
    # define Time of Emergence (ToE) from last time index at which signal is larger than mult*noise
    # signal is [time,space]
    # noise is [space]
    # mult is float
    toe_wrk = npy.ma.ones([timN,levN*latN])*1. # init toe1_wrk array to 1
    signaltile = npy.reshape(np.tile(signal,timN),(timN,levN*latN)) # repeat timN
    toe_idx = npy.argwhere(abs(tvarha-tvarhna) >= multStd*stdvarams) # find indices of points > stdev
    toe_wrk[toe1_idx[:,0],toe1_idx[:,1]] = 0. # set points in toe1_wrk to zero
    toe = timN-npy.flipud(toe1_wrk).argmax(axis=0) # compute ToE1

    varam = np.reshape(toe1,(levN,latN))

    return toe
def mmeAveMsk2D(listFiles, years, inDir, outDir, outFile, timeInt, mme, debug=True):
    '''
    The mmeAveMsk2D() function averages rhon/lat density bined files with differing masks
    It ouputs
     - the MME
     - a percentage of non-masked bins
     - the sign agreement of period2-period1 differences

    Author:    Eric Guilyardi : Eric.Guilyardi@locean-ipsl.upmc.fr

    Created on Tue Nov 25 13:56:20 CET 2014

    Inputs:
    -------
    - listFiles(str)         - the list of files to be averaged
    - years(t1,t2)           - years for slice read
    - inDir(str)             - input directory where files are stored
    - outDir(str)            - output directory
    - outFile(str)           - output file
    - timeInt(2xindices)     - indices of init period to compare with (e.g. [1,20])
    - mme(bool)              - multi-model mean (will read in single model ensemble stats)
    - debug <optional>       - boolean value

    Notes:
    -----
    - EG 25 Nov 2014   - Initial function write
    - EG 27 Nov 2014   - Rewrite with loop on variables
    - EG 06 Dec 2014   - Added agreement on difference with init period - save as <var>Agree
    - EG 07 Dec 2014   - Read bowl to remove points above bowl - save as <var>Bowl

    - TO DO :
                 - optimization of loops
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
    # Bound of period average to remove
    peri1 = timeInt[0]
    peri2 = timeInt[1]
    fi      = cdm.open(inDir+'/'+listFiles[0])
    isond0  = fi('isondepth',time = slice(t1,t2)) ; # Create variable handle
    # Get grid objects
    axesList = isond0.getAxisList()
    sigmaGrd = isond0.getLevel()
    # Declare and open files for writing
    if os.path.isfile(outDir+'/'+outFile):
        os.remove(outDir+'/'+outFile)
    outFile_f = cdm.open(outDir+'/'+outFile,'w')

    latN = isond0.shape[3]
    levN = isond0.shape[2]
    basN = isond0.shape[1]
    timN = isond0.shape[0]
    runN = len(listFiles)

    print ' Number of members:',len(listFiles)

    valmask = isond0.missing_value

    varList = ['isondepth','isonpers','isonso','isonthetao','isonthick','isonvol']
    varFill = [0.,0.,valmask,valmask,0.,0.]
    #varList = ['isondepth']
    #varList = ['isonthetao']

    # init arrays
    percent  = npy.ma.ones([runN,timN,basN,levN,latN], dtype='float32')*0.
    minbowl  = npy.ma.ones([basN,latN], dtype='float32')*1000.
    varbowl  = npy.ma.ones([runN,timN,basN,latN], dtype='float32')*1.
    # init time axis
    time       = cdm.createAxis(npy.float32(range(timN)))
    time.id    = 'time'
    time.units = 'years since 1861'
    time.designateTime()

    # loop on variables
    for iv,var in enumerate(varList):

        # Array inits
        isonvar  = npy.ma.ones([runN,timN,basN,levN,latN], dtype='float32')*valmask
        vardiff,varbowl2D = [npy.ma.ones(npy.ma.shape(isonvar)) for _ in range(2)]
        varones  = npy.ma.ones([runN,timN,basN,levN,latN], dtype='float32')*1.
        print ' Variable ',iv, var
        # loop over files to fill up array
        for i,file in enumerate(listFiles):
            ft      = cdm.open(inDir+'/'+file)
            timeax  = ft.getAxis('time')
            file1d  =  replace(inDir+'/'+file,'2D','1D')
            if os.path.isfile(file1d):
                f1d = cdm.open(file1d)
            else:
                print 'ERROR:',file1d,'missing (if mme, run 1D first)'
                sys.exit(1)
            if i == 0:
                tmax0 = timeax.shape[0]
                tmax = timeax.shape[0]
            if tmax != tmax0:
                print 'wrong time axis: exiting...'
                return
            # read array
            isonRead = ft(var,time = slice(t1,t2))
            if varFill[iv] != valmask:
                isonvar[i,...] = isonRead.filled(varFill[iv])
            else:
                isonvar[i,...] = isonRead
            # compute percentage of non-masked points accros MME
            if iv == 0:
                maskvar = mv.masked_values(isonRead.data,valmask).mask
                percent[i,...] = npy.float32(npy.equal(maskvar,0))
            if mme:
                # if mme then just average Bowl and Agree fields
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
                    bowlRead = f1d('ptopsigma',time = slice(t1,t2))
                    varbowl[i,...] = bowlRead

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

        # average accross members
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
                bowlRead = f1d('ptopsigma',time = slice(t1,t2))
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
            if iv == 0:
                siglimit = cdu.averager(varbowl, axis=0)
                siglimit = cdu.averager(siglimit, axis=0) - delta_rho
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

        outFile_f.write(    isonave.astype('float32'))
        outFile_f.write(isonavediff.astype('float32'))
        outFile_f.write(isonavebowl.astype('float32'))
        if mme:
            isonvarstd = cdm.createVariable(isonVarStd , axes =[time,axesList[1],axesList[2],axesList[3]] , id = isonRead.id+'ModStd')
            isonvarstd.long_name = isonRead.long_name+' intermodel std'
            isonvarstd.units     = isonRead.units
            outFile_f.write(isonvarstd.astype('float32'))

    # <--- end of loop on variables 

    outFile_f.close()
    fi.close()

def mmeAveMsk1D(listFiles, years, inDir, outDir, outFile, timeInt, mme, debug=True):
    '''
    The mmeAveMsk1D() function averages rhon or scalar density bined files with differing masks
    It ouputs the MME and a percentage of non-masked bins
    
    Created on Tue Nov 25 13:56:20 CET 2014

    Inputs:
    -------
    - listFiles(str)         - the list of files to be averaged
    - years(t1,t2)           - years for slice read
    - inDir(str)             - input directory where files are stored
    - outDir(str)            - output directory
    - outFile(str)           - output file
    - timeInt(2xindices)     - indices of init period to compare with (e.g. [1,20])
    - mme(bool)              - multi-model mean (will read in single model ensemble stats)
    - debug <optional>       - boolean value

    Notes:
    -----
    - EG 25 Nov 2014   - Initial function write
    - EG  9 Dec 2014   - Add agreement on difference with init period - save as <var>Agree

    To do:
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
    # Bound of period average to remove
    peri1 = timeInt[0]
    peri2 = timeInt[1]
    fi      = cdm.open(inDir+'/'+listFiles[0])
    ptopd0  = fi('ptopdepth',time=slice(t1,t2)) ; # Create variable handle
    # Get grid objects
    axesList = ptopd0.getAxisList()
    # Declare and open files for writing
    if os.path.isfile(outDir+'/'+outFile):
        os.remove(outDir+'/'+outFile)
    outFile_f = cdm.open(outDir+'/'+outFile,'w')

    latN = ptopd0.shape[2]
    basN = ptopd0.shape[1]
    timN = ptopd0.shape[0]
    runN = len(listFiles)

    print ' Number of members:',len(listFiles)

    valmask = ptopd0.missing_value

    # loop on variables

    varList = ['ptopdepth','ptopsigma','ptopso','ptopthetao','volpers','salpers','tempers']
    #varFill = [0.,0.,valmask,valmask,0.,0.,0.,valmask,valmask]
    varFill = [valmask,valmask,valmask,valmask,valmask,valmask,valmask,valmask,valmask]
    varDim  = [1,1,1,1,0,0,0]
    #varList = ['ptopdepth']

    valmask = ptopd0.missing_value

    # init percent array
    percent  = npy.ma.ones([runN,timN,basN,latN], dtype='float32')*0.
    # init time axis
    time       = cdm.createAxis(npy.float32(range(timN)))
    time.id    = 'time'
    time.units = 'years since 1861'
    time.designateTime()

    axis1D = [time,axesList[1],axesList[2]]
    axis0D = [time,axesList[1]]

    # loop on 1D variables
    for iv,var in enumerate(varList):

        # Array inits
        if varDim[iv] == 1:
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
            ft      = cdm.open(inDir+'/'+file)
            timeax  = ft.getAxis('time')
            if ic == 0:
                tmax0 = timeax.shape[0]
                tmax = timeax.shape[0]
            if tmax != tmax0:
                print "wrong time axis: exiting..."
                return
            # read array
            isonRead = ft(var,time = slice(t1,t2))
            #print isonRead.shape, isonvar.shape, file
            if varFill[iv] != valmask:
                isonvar[ic,...] = isonRead.filled(varFill[iv])
            else:
                isonvar[ic,...] = isonRead

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
    # <--- end of loop on variables 

    outFile_f.close()
    fi.close()

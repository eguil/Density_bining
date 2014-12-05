import gc,os,resource,timeit,glob,re,math
import cdms2 as cdm
import cdutil as cdu
import MV2 as mv
import numpy as npy
from string import replace
import time as timc

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

def mmeAveMsk2D(listFiles, years, indDir, outDir, outFile, timeInt, mme, debug=True):
    '''
    The mmeAveMsk2D() function averages rhon/lat density bined files with differing masks
    It ouputs the MME, a percentage of non-masked bins, and the sign agreement of period2-period1 differences.
    
    Author:    Eric Guilyardi : Eric.Guilyardi@locean-ipsl.upmc.fr

    Created on Tue Nov 25 13:56:20 CET 2014

    Inputs:
    -------
    - listFiles(str)         - the list of files to be averaged
    - years(t1,t2)           - years for slice read
    - inDir(str)             - input directory where files are stored
    - outDir(str)            - output directory
    - outFile(str)           - output file
    - timeInt(4xindices)     - indices of period 1 and period 2 (e.g. [1,5,140,146])
    - mme(bool)              - multi-model mean (will read in single model ensemble stats)
    - debug <optional>       - boolean value

    Notes:
    -----
    - EG 25 Nov 2014   - Initial function write
    - EG 27 Nov 2014   - Rewrite with loop on variables
    - EG 06 Dec 2014   - Added agreement on difference between period 1 and period 2
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
    # Periods 1 and 2 bounds
    per1i = timeInt[0]
    per1f = timeInt[1]
    per2i = timeInt[2]
    per2f = timeInt[3]

    fi      = cdm.open(indir+'/'+listFiles[0])
    isond0  = fi('isondepth',time = slice(t1,t2)) ; # Create variable handle
    # Get grid objects
    axesList = isond0.getAxisList()
    # Declare and open files for writing
    if os.path.isfile(outDir+'/'+outFile):
        os.remove(outDir+'/'+outFile)
    outFile_f = cdm.open(outDir+'/'+outFile,'w')

    latN = isond0.shape[3]
    levN = isond0.shape[2]
    basN = isond0.shape[1]
    timN = isond0.shape[0]
    runN = len(listFiles)  

    print 'Number of runs:',runN

    valmask = isond0.missing_value

    varList = ['isondepth','isonpers','isonso','isonthetao','isonthick','isonvol']
    varFill = [0.,0.,valmask,valmask,0.,0.]
    #varList = ['isonvol']

    # init percent array
    percent  = npy.ma.ones([runN,timN,basN,levN,latN], dtype='float32')*0.
    # init time axis
    time       = cdm.createAxis(npy.float32(range(timN)))
    time.id    = 'time'
    time.units = 'years since 1851'
    time.designateTime()

    # loop on variables
    for iv,var in enumerate(varList):

        # Array inits
        isonvar  = npy.ma.ones([runN,timN,basN,levN,latN], dtype='float32')*valmask
        vardiff  = npy.ma.ones([runN,basN,levN,latN], dtype='float32')*valmask
        varones  = npy.ma.ones([runN,basN,levN,latN], dtype='float32')*1.
        print ' Variable ',iv, var
        # loop over files to fill up array
        for i,file in enumerate(listFiles):
            print i, file
            ft      = cdm.open(indir+'/'+file)
            timeax  = ft.getAxis('time')
            if i == 0:
                tmax0 = timeax.shape[0]
                tmax = timeax.shape[0]
            if tmax != tmax0:
                print 'wrong time axis: exiting...'
                return
            # read array
            isonRead = ft(var,time = slice(t1,t2))
            if mme:
                varst = var+'Agree'
                isonReadst = ft(varst,time = slice(t1,t2))
            if varFill[iv] != valmask:
                isonvar[i,...] = isonRead.filled(varFill[iv])
            else:
                isonvar[i,...] = isonRead
            ft.close()
            # compute percentage of non-masked points accros MME
            if iv == 0:
                maskvar = mv.masked_values(isonRead.data,valmask).mask
                percent[i,...] = npy.float32(npy.equal(maskvar,0))
            # Compute difference between period 1 and period 2, use mask of last month
            if mme:
                vardiff[i,...] = isonReadst
            else:
                vardiff[i,...] = cdu.averager(isonvar[i,per2i:per2f,...],axis=0) - cdu.averager(isonvar[i,per1i:per1f,...],axis=0)
                vardiff[i,...].mask = isonvar[i,-1,...].mask
                

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
        else:
            vardiffsg = npy.copysign(varones,vardiff)
            # average signs
            vardiffsgSum = cdu.averager(vardiffsg, axis=0)
            vardiffsgSum._FillValue = valmask
            vardiffsgSum.mask = percentw.mask

        # average
        isonVarAve = cdu.averager(isonvar, axis=0)
        isonVarAve = cdm.createVariable(isonVarAve , axes =[time,axesList[1],axesList[2],axesList[3]] , id = 'foo')
        # mask
        if varFill[iv] == valmask:
            isonVarAve = maskVal(isonVarAve, valmask)

        isonVarAve.mask = percentw.mask

        # Write
        isonave = cdm.createVariable(isonVarAve, axes = [time,axesList[1],axesList[2],axesList[3]], id = isonRead.id)
        isonave.long_name = isonRead.long_name
        isonave.units     = isonRead.units
        isonavediff = cdm.createVariable(vardiffsgSum, axes = [axesList[1],axesList[2],axesList[3]], id = isonRead.id+'Agree')
        isonavediff.long_name = isonRead.long_name
        isonavediff.units     = isonRead.units

        outFile_f.write(    isonave.astype('float32'))
        outFile_f.write(isonavediff.astype('float32'))
    # <--- end of loop on variables 

    outFile_f.close()
    fi.close()

def mmeAveMsk1D(listFiles, years, indDir, outDir, outFile, debug=True):
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
    - debug <optional>       - boolean value

    Notes:
    -----
    - EG 25 Nov 2014   - Initial function write
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
    fi      = cdm.open(indir+'/'+listFiles[0])
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

    print 'Number of runs:',runN

    valmask = ptopd0.missing_value

    # loop on variables

    varList = ['ptopdepth','ptopsigma','ptopso','ptopthetao','volpers','salpers','tempers']
    varFill = [0.,0.,valmask,valmask,0.,0.,0.,valmask,valmask]
    varFill = [valmask,valmask,valmask,valmask,valmask,valmask,valmask,valmask,valmask]
    varDim  = [1,1,1,1,0,0,0]
    #varList = ['ptopdepth']

    valmask = ptopd0.missing_value

    # init percent array
    percent  = npy.ma.ones([runN,timN,basN,latN], dtype='float32')*0.
    # init time axis
    time       = cdm.createAxis(npy.float32(range(timN)))
    time.id    = 'time'
    time.units = 'years since 1851'
    time.designateTime()

    axis1D = [time,axesList[1],axesList[2]]
    axis0D = [time,axesList[1]]

    # loop on 1D variables
    for iv,var in enumerate(varList):

        # Array inits
        if varDim[iv] == 1:
            isonvar = npy.ma.ones([runN,timN,basN,latN], dtype='float32')*valmask
            axisVar = axis1D
        else:
            isonvar = npy.ma.ones([runN,timN,basN], dtype='float32')*valmask
            axisVar = axis0D
        print ' Variable ',iv, var, varDim[iv]
        print isonvar.shape
        # loop over files to fill up array
        for i,file in enumerate(listFiles):
            print i, file
            ft      = cdm.open(indir+'/'+file)
            timeax  = ft.getAxis('time')
            if i == 0:
                tmax0 = timeax.shape[0]
                tmax = timeax.shape[0]
            if tmax != tmax0:
                print "wrong time axis: exiting..."
                return
            # read array
            isonRead = ft(var,time = slice(t1,t2))
            print isonRead.shape
            if varFill[iv] != valmask:
                isonvar[i,...] = isonRead.filled(varFill[iv])
            else:
                isonvar[i,...] = isonRead
            ft.close()

            # compute percentage of non-masked points accros MME
            if iv == 0:
                maskvar = mv.masked_values(isonRead.data,valmask).mask
                percent[i,...] = npy.float32(npy.equal(maskvar,0))

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

        # average
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

        outFile_f.write(isonave.astype('float32'))
    # <--- end of loop on variables 

    outFile_f.close()
    fi.close()
#
# ----------------------------------------------------------------------------
# Work

# Model ensemble mean

#twoD = False
#oneD = True
twoD = True
oneD = False
mme  = False
mm = True 

exper  = 'historical'
models = ['ACCESS1-0','ACCESS1-3','BNU-ESM','CCSM4','CESM1-BGC','EC-EARTH','FGOALS-s2','GFDL-CM2p1','GISS-E2-R','HadCM3','HadGEM2-CC','HadGEM2-ES','IPSL-CM5A-LR','IPSL-CM5A-MR','IPSL-CM5B-LR','MIROC-ESM-CHEM','MIROC-ESM']
years = [[10,156],[10,156],[10,156],[10,156],[10,156],[10,156],[10,156],[10,156],[10,156],[0,146],[0,146],[0,146],[10,156],[10,156],[10,156],[10,156],[10,156]]
#models = ['ACCESS1-3']#,'ACCESS1-3']
#models = ['EC-EARTH']
#years = [[10,156]]

# Years for difference
iniyear = 1861
per1i = (1950-iniyear)+1
per1f = (1950-iniyear)+2
per2i = (2000-iniyear)+1
per2f = (2000-iniyear)+2
timeInt=[per1i,per1f,per2i,per2f]
print timeInt
indir  = '/Users/ericg/Projets/Density_bining/Prod_density_nov14/z_individual'
outdir = '/Users/ericg/Projets/Density_bining/Prod_density_nov14/test_mme'
listens = []
listens1 = []
for i,mod in enumerate(models):
    os.chdir(indir)
    print i,mod, 'slice', years[i]
    listf  = glob.glob('cmip5.'+mod+'.*2D*')
    listf1 = glob.glob('cmip5.'+mod+'.*1D*')
    start = listf[0].find(exper)+len(exper)
    end = listf[0].find('.an.')
    rip = listf[0][start:end]
    outFile = replace(listf[0],rip,'.ensm')
    outFile1 = replace(outFile,'2D','1D')
    listens.append(outFile)
    listens1.append(outFile1)
    print outFile
    print outFile1
    if mm:
        if twoD:
            mmeAveMsk2D(listf,years[i],indir,outdir,outFile,timeInt,mme)
        if oneD:
            mmeAveMsk1D(listf1,years[i],indir,outdir,outFile1,timeInt,mme)

if mme:
    # MME
    indir  = outdir
    if twoD:
        outFile = 'cmip5.multimodel.historical.ensm.an.ocn.Omon.density_zon2D.nc'
        mmeAveMsk2D(listens,[0,146],indir,outdir,outFile,timeInt,mme)
    if oneD:
        outFile1 = 'cmip5.multimodel.historical.ensm.an.ocn.Omon.density_zon1D.nc'
        mmeAveMsk1D(listens1,[0,146],indir,outdir,outFile1,timeInt,mme)

modelsurf = ['ACCESS1-0','ACCESS1-3','CMCC-CESM','CMCC-CM','CMCC-CMS','CNRM-CM5','CSIRO-Mk3-6-0','EC-EARTH','FGOALS-s2','GFDL-ESM2G','GISS-E2-R-CC','GISS-E2-R','MIROC5','MIROC-ESM-CHEM','MIROC-ESM','MPI-ESM-LR','MPI-ESM-MR','MPI-ESM-P','NorESM1-ME','NorESM1-M']

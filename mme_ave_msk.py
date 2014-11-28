import gc,os,resource,timeit,glob,re
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

def mmeAveMsk(listFiles, years, indDir, outDir, outFile, debug=True):
    '''
    The mmeAveMsk() function averages density bined files with differing masks
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
    - EG 27 Nov 2014   - Rewrite with loop on variables
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

    # loop on variables

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

    for iv,var in enumerate(varList):

        # Array inits
        isonvar  = npy.ma.ones([runN,timN,basN,levN,latN], dtype='float32')*valmask
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
                print "wrong time axis: exiting..."
                return
            # read array
            isonRead = ft(var,time = slice(t1,t2))
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
            percentw = cdm.createVariable(percenta, axes = [time,axesList[1],axesList[2],axesList[3]], id = 'isonpercent')
            percentw._FillValue = valmask
            percentw.long_name = 'percentage of MME bin'
            percentw.units     = '%'
            outFile_f.write(percentw.astype('float32'))

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

        outFile_f.write(isonave.astype('float32'))
    # <--- end of loop on variables 

    outFile_f.close()
    fi.close()

def mmeAveMsk1D(listFiles, years, indDir, outDir, outFile, debug=True):
    '''
    The mmeAveMsk1D() function averages density bined files with differing masks
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

    latN    = ptopd0.shape[2]
    basN    = ptopd0.shape[1]
    timN    = ptopd0.shape[0]

    valmask = ptopd0.missing_value

    # Array inits
    ptopdc  = npy.ma.ones([timN,basN,latN], dtype='float32')*0.
    ptopgc,ptopsc,ptoptc,percent = [npy.ma.ones(npy.shape(ptopdc)) for _ in range(4)]
    volpc   = npy.ma.ones([timN,basN], dtype='float32')*0.
    salpc,tempc = [npy.ma.ones(npy.shape(volpc)) for _ in range(2)]

    count = 0
    # loop over files
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
        # read arrays
        ptopd = ft('ptopdepth',time=slice(t1,t2))
        ptopg = ft('ptopsigma',time=slice(t1,t2))
        ptops = ft('ptopso',time=slice(t1,t2))
        ptopt = ft('ptopthetao',time=slice(t1,t2))
        volp  = ft('volpers',time=slice(t1,t2))
        salp  = ft('salpers',time=slice(t1,t2))
        temp  = ft('tempers',time=slice(t1,t2))
        # compute percentage of non-masked points accros MME
        maskvar = mv.masked_values(ptopd.data,valmask).mask
        nomask  = npy.equal(maskvar,0)
        if i == 0:
            percent = npy.float32(nomask)
        else:
            percent = percent + npy.float32(nomask)

        # accumulate
        ptopdc  = ptopdc + ptopd#.filled(0.)
        ptopgc  = ptopgc + ptopg#.filled(0.)
        ptopsc  = ptopsc + ptops#.filled(0.)
        ptoptc  = ptoptc + ptopt#.filled(0.)
        volpc   = volpc  + volp #.filled(0.)
        salpc   = salpc  + salp #.filled(0.)
        tempc   = tempc  + temp #.filled(0.)

        ft.close()
        count = count + 1
        # <--- end file loop 


    # Average & mask

    ptopdc = ptopdc/float(count)
    ptopgc = ptopgc/float(count)
    ptopsc = ptopsc/float(count)
    ptoptc = ptoptc/float(count)

    ptopdc = mv.masked_where(ptopdc > 10000,ptopdc)
    ptopgc.mask = ptopdc.mask
    ptopsc.mask = ptopdc.mask
    ptoptc.mask = ptopdc.mask

    ptopdc = maskVal(ptopdc, valmask)
    ptopgc = maskVal(ptopgc, valmask)
    ptopsc = maskVal(ptopsc, valmask)
    ptoptc = maskVal(ptoptc, valmask)

    volpc  = volpc/float(count)
    salpc  = salpc/float(count)
    tempc  = tempc/float(count)

    percent = percent/float(count)*100.

    #print axesList
    time       = cdm.createAxis(npy.float32(range(timN)))
    time.id    = 'time'
    time.units = 'years since 1851'
    time.designateTime()

    # Write
    ptopdcw = cdm.createVariable(ptopdc, axes = [time,axesList[1],axesList[2]], id = 'ptopdepth')
    ptopdcw.long_name = ptopd.long_name
    ptopdcw.units     = ptopd.units
    ptopgcw = cdm.createVariable(ptopgc, axes = [time,axesList[1],axesList[2]], id = 'ptopsigma')
    ptopgcw.long_name = ptopg.long_name
    ptopgcw.units     = ptopg.units
    ptopscw = cdm.createVariable(ptopsc, axes = [time,axesList[1],axesList[2]], id = 'ptopso')
    ptopscw.long_name = ptops.long_name
    ptopscw.units     = ptops.units
    ptoptcw = cdm.createVariable(ptoptc, axes = [time,axesList[1],axesList[2]], id = 'ptopthetao')
    ptoptcw.long_name = ptopt.long_name
    ptoptcw.units     = ptopt.units

    volpw = cdm.createVariable(volpc, axes = [time,axesList[1]], id = 'volpers')
    volpw.long_name = volp.long_name
    volpw.units     = volp.units
    salpw = cdm.createVariable(salpc, axes = [time,axesList[1]], id = 'salpers')
    salpw.long_name = salp.long_name
    salpw.units     = salp.units
    tempw = cdm.createVariable(tempc, axes = [time,axesList[1]], id = 'tempers')
    tempw.long_name = temp.long_name
    tempw.units     = temp.units

    percenta  = cdm.createVariable(percent, axes = [time,axesList[1],axesList[2]], id = 'isonpercent')
    percenta.long_name = 'percentage of MME bin'
    percenta.units     = '%'

    outFile_f.write(ptopdcw.astype('float32'))
    outFile_f.write(ptopgcw.astype('float32'))
    outFile_f.write(ptopscw.astype('float32'))
    outFile_f.write(ptoptcw.astype('float32'))
    outFile_f.write(volpw.astype('float32'))
    outFile_f.write(salpw.astype('float32'))
    outFile_f.write(tempw.astype('float32'))
    outFile_f.write(percenta.astype('float32'))

    outFile_f.close()
    fi.close

# Work

# Model ensemble mean

#twoD = False
#oneD = True
twoD = True
oneD = False
mm  = False
mme = True 

exper  = 'historical'
models = ['ACCESS1-0','ACCESS1-3','BNU-ESM','CCSM4','CESM1-BGC','EC-EARTH','FGOALS-s2','GFDL-CM2p1','GISS-E2-R','HadCM3','HadGEM2-CC','HadGEM2-ES','IPSL-CM5A-LR','IPSL-CM5A-MR','IPSL-CM5B-LR','MIROC-ESM-CHEM','MIROC-ESM']
years = [[10,156],[10,156],[10,156],[10,156],[10,156],[10,156],[10,156],[10,156],[10,156],[0,146],[0,146],[0,146],[10,156],[10,156],[10,156],[10,156],[10,156]]
#models = ['ACCESS1-3']#,'ACCESS1-3']
#models = ['EC-EARTH']
#years = [[10,156]]
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
            mmeAveMsk  (listf,years[i],indir,outdir,outFile)
        if oneD:
            mmeAveMsk1D(listf1,years[i],indir,outdir,outFile1)

if mme:
    # MME
    indir  = outdir
    if twoD:
        outFile = 'cmip5.multimodel.historical.ensm.an.ocn.Omon.density_zon2D.nc'
        mmeAveMsk(listens,[0,146],indir,outdir,outFile)
    if oneD:
        outFile1 = 'cmip5.multimodel.historical.ensm.an.ocn.Omon.density_zon1D.nc'
        mmeAveMsk1D(listens1,[0,146],indir,outdir,outFile1)

modelsurf = ['ACCESS1-0','ACCESS1-3','CMCC-CESM','CMCC-CM','CMCC-CMS','CNRM-CM5','CSIRO-Mk3-6-0','EC-EARTH','FGOALS-s2','GFDL-ESM2G','GISS-E2-R-CC','GISS-E2-R','MIROC5','MIROC-ESM-CHEM','MIROC-ESM','MPI-ESM-LR','MPI-ESM-MR','MPI-ESM-P','NorESM1-ME','NorESM1-M']

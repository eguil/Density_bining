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

def mmeAveMsk(listFiles, indDir, outDir, outFile, debug=True):
    '''
    The mmeAveMsk() function averages density bined files with differing masks
    It ouputs the MME and a percentage of non-masked bins
    
    Created on Tue Nov 25 13:56:20 CET 2014

    Inputs:
    -------
    - listFiles(str)         - the list of files to be averaged
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
    fi      = cdm.open(indir+'/'+listFiles[0])
    isond0  = fi('isondepth') ; # Create variable handle
    # Get grid objects
    axesList = isond0.getAxisList()
    # Declare and open files for writing
    if os.path.isfile(outDir+'/'+outFile):
        os.remove(outDir+'/'+outFile)
    outFile_f = cdm.open(outDir+'/'+outFile,'w')

    latN    = isond0.shape[3]
    levN    = isond0.shape[2]
    basN    = isond0.shape[1]
    timN    = isond0.shape[0]

    valmask = isond0.missing_value

    # Array inits
    isondc  = npy.ma.ones([timN,basN,levN,latN], dtype='float32')*0.
    isonpc,isonsc,isontc,isonhc,isonvc,percent = [npy.ma.ones(npy.shape(isondc)) for _ in range(6)]

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
        isond = ft('isondepth')
        isonp = ft('isonpers')
        isons = ft('isonso')
        isont = ft('isonthetao')
        isonh = ft('isonthick')
        isonv = ft('isonvol')
        # compute percentage of non-masked points accros MME
        maskvar = mv.masked_values(isond.data,valmask).mask
        nomask  = npy.equal(maskvar,0)
        if i == 0:
            percent = npy.float32(nomask)
        else:
            percent = percent + npy.float32(nomask)

        # accumulate
        isondc  = isondc + isond.filled(0.)
        isonpc  = isonpc + isonp.filled(0.)
        isonsc  = isonsc + isons.filled(0.)
        isontc  = isontc + isont.filled(0.)
        isonhc  = isonhc + isonh.filled(0.)
        isonvc  = isonvc + isonv.filled(0.)

        ft.close()
        count = count + 1
        # <--- end file loop 


    # Average
    print 'count = ',count 

    isondc = isondc/float(count)
    isonpc = isonpc/float(count)
    isonsc = isonsc/float(count)
    isontc = isontc/float(count)
    isonhc = isonhc/float(count)
    isonvc = isonvc/float(count)
    isondc = mv.masked_where(isondc < 0.01,isondc)
    isonpc.mask = isondc.mask
    isonsc.mask = isondc.mask
    isontc.mask = isondc.mask
    isonhc.mask = isondc.mask
    isonvc.mask = isondc.mask

    isondc = maskVal(isondc, valmask)
    isonpc = maskVal(isonpc, valmask)
    isonsc = maskVal(isonsc, valmask)
    isontc = maskVal(isontc, valmask)
    isonhc = maskVal(isonhc, valmask)
    isonvc = maskVal(isonvc, valmask)

    percent = percent/float(count)*100.

    #print axesList
    time       = cdm.createAxis(npy.float32(range(timN)))
    time.id    = 'time'
    time.units = 'years since 1851'
    time.designateTime()

    # Write
    isondcw = cdm.createVariable(isondc, axes = [time,axesList[1],axesList[2],axesList[3]], id = 'isondepth')
    isondcw.long_name = isond.long_name
    isondcw.units     = isond.units
    isonpcw = cdm.createVariable(isonpc, axes = [time,axesList[1],axesList[2],axesList[3]], id = 'isonpers')
    isonpcw.long_name = isonp.long_name
    isonpcw.units     = isonp.units
    isonscw = cdm.createVariable(isonsc, axes = [time,axesList[1],axesList[2],axesList[3]], id = 'isonso')
    isonscw.long_name = isons.long_name
    isonscw.units     = isons.units
    isontcw = cdm.createVariable(isontc, axes = [time,axesList[1],axesList[2],axesList[3]], id = 'isonthetao')
    isontcw.long_name = isont.long_name
    isontcw.units     = isont.units
    isonhcw = cdm.createVariable(isonhc, axes = [time,axesList[1],axesList[2],axesList[3]], id = 'isonthick')
    isonhcw.long_name = isonh.long_name
    isonhcw.units     = isonh.units
    isonvcw = cdm.createVariable(isonvc, axes = [time,axesList[1],axesList[2],axesList[3]], id = 'isoncol')
    isonvcw.long_name = isonv.long_name
    isonvcw.units     = isonv.units

    percenta  = cdm.createVariable(percent, axes = [time,axesList[1],axesList[2],axesList[3]], id = 'isonpercent')
    percenta.long_name = 'percentage of MME bin'
    percenta.units     = '%'

    outFile_f.write(isondcw.astype('float32'))
    outFile_f.write(isonpcw.astype('float32'))
    outFile_f.write(isonscw.astype('float32'))
    outFile_f.write(isontcw.astype('float32'))
    outFile_f.write(isonhcw.astype('float32'))
    outFile_f.write(isonvcw.astype('float32'))
    outFile_f.write(percenta.astype('float32'))

    outFile_f.close()
    fi.close()


# Work

# Model ensemble mean

mm  = True
mme = False 

exper  = 'historical'
models = ['ACCESS1-0','ACCESS1-3','BNU-ESM','CCSM4','CESM1-BGC','EC-EARTH','FGOALS-s2','GFDL-CM2p1','GISS-E2-R','HadCM3','HadGEM2-CC','HadGEM2-ES','IPSL-CM5A-LR','IPSL-CM5A-MR','IPSL-CM5B-LR','MIROC-ESM-CHEM','MIROC-ESM']
indir  = '/Users/ericg/Projets/Density_bining/Prod_density_nov14/z_individual'
outdir = '/Users/ericg/Projets/Density_bining/Prod_density_nov14'
listens = []
for i,mod in enumerate(models):
    os.chdir(indir)
    print i,mod
    listf = glob.glob('cmip5.'+mod+'*2D*')
    start = listf[0].find(exper)+len(exper)
    end = listf[0].find('.an.')
    rip = listf[0][start:end]
    outFile = replace(listf[0],rip,'.ensm')
    listens.append(outFile)
    print outFile
    if mm:
        mmeAveMsk(listf,indir,outdir,outFile)

if mme:
    # MME
    indir  = outdir
    outfile = cmip5.multimodel.historical.ensm.an.ocn.Omon.density_zon2D.nc
    mmeAveMsk(listens,indir,outdir,outFile)

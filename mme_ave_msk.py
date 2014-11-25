import gc,os,resource,timeit
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

def mmeAveMsk(listFiles, indDir, outFile, debug=True):
    '''
    The mmeAveMsk() function averages density bined files with differing masks
    It ouputs the MME and a percentage of non-masked bins
    
    Created on Tue Nov 25 13:56:20 CET 2014

    Inputs:
    -------
    - listFiles(str)         - the list of files to be averaged
    - inDir(str)             - input directory where files are stored
    - outFile(str)           - output file with full path specified
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
    isond   = fi('isondepth') ; # Create variable handle
    # Get grid objects
    axesList = isond.getAxisList()
    # Declare and open files for writing
    if os.path.isfile(outFile):
        os.remove(outFile)
    outFile_f = cdm.open(outFile,'w')

    latN    = isond.shape[3]
    levN    = isond.shape[2]
    basN    = isond.shape[1]
    timN    = isond.shape[0]

    valmask = isond.missing_value
    print valmask

    # Array inits
    isondc  = npy.ma.ones([timN,basN,levN,latN], dtype='float32')*0.
    isonpc,isonsc,isontc,isonhc,isonvc = [npy.ma.ones(npy.shape(isondc)) for _ in range(5)]

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

        # accumulate

        isondc = isondc + isond

        print isondc.shape

        ft.close()
        count = count + 1
        # <--- end file loop 


    # Average
    print 'count = ',count 

    isondc = isondc/float(count)
    #isondc = mv.masked_where(isondc > 10000.,isondc)
    #isondc.data = isondc.filled(valmask)
    #isondc.mask = isond.mask
    isondc = maskVal(isondc, valmask)

    #print timeax,axesList[1],axesList[2],axesList[3]

    # Write
    isondcw = cdm.createVariable(isondc, axes = axesList, id = 'isondepth')
    isondcw.long_name = isond.long_name
    isondcw.units     = isond.units

    print isondcw.shape

    outFile_f.write(isond.astype('float32'))

    outFile_f.close()
    fi.close()



listf = ['cmip5.ACCESS1-3.historical.r1i1p1.an.ocn.Omon.density.ver-1_zon2D.nc','cmip5.ACCESS1-3.historical.r2i1p1.an.ocn.Omon.density.ver-v2_zon2D.nc','cmip5.ACCESS1-3.historical.r3i1p1.an.ocn.Omon.density.ver-v3_zon2D.nc']
listf = ['cmip5.ACCESS1-3.historical.r1i1p1.an.ocn.Omon.density.ver-1_zon2D.nc']
indir = 'Prod_density_nov14/z_individual'
outf  = 'testMME.nc'

mmeAveMsk(listf,indir,outf)

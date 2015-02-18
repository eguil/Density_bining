# -*- coding: utf-8 -*-
"""
Created on Mon May 27 16:49:02 2013

Paul J. Durack 27th May 2013

PJD  2 Jun 2013     - Added clearall and outer_locals functions
PJD 13 Jun 2013     - Updated global_att_write to globalAttWrite
PJD 13 Jun 2013     - Added writeToLog function
PJD 26 Jun 2013     - Added fixInterpAxis function
PJD 18 Jul 2013     - Sphinx docs/syntax http://pythonhosted.org/an_example_pypi_project/sphinx.html
PJD 23 Jul 2013     - Added fixVarUnits function
PJD 25 Jul 2013     - Updated fixVarUnits function to print and log changes
PJD  5 Aug 2013     - Added fitPolynomial function following Pete G's code example
PJD  9 Aug 2013     - Added writePacked function
PJD  9 Aug 2013     - Added keyboard function
PJD 22 Aug 2013     - Added setTimeBoundsYearly() to fixInterpAxis
PJD  1 Apr 2014     - Added trimModelList
PJD 20 Aug 2014     - Added mkDirNoOSErr and sysCallTimeout functions
PJD 13 Oct 2014     - Added getGitInfo function
                    - TODO: Consider implementing multivariate polynomial regression:
                      https://github.com/mrocklin/multipolyfit

This library contains all functions written to replicate matlab functionality in python

@author: durack1
"""

## Import common modules ##
import cdat_info,cdtime,code,datetime,errno,gc,inspect,os,pytz,re,string,sys,time
import cdms2 as cdm
import cdutil as cdu
#import genutil as genu
#import matplotlib as plt
import MV2 as mv
import numpy as np
import subprocess
#import scipy as sp
from numpy.core.fromnumeric import shape
from socket import gethostname
from string import replace
# Consider modules listed in /work/durack1/Shared/130103_data_SteveGriffies/130523_mplib_tips/importNPB.py
##

## Specify UVCDAT specific stuff ##
# Turn off cdat ping reporting - Does this speed up Spyder?
cdat_info.ping = False
# Set netcdf file criterion - turned on from default 0s
cdm.setCompressionWarnings(0) ; # Suppress warnings
cdm.setNetcdfShuffleFlag(0)
cdm.setNetcdfDeflateFlag(1)
cdm.setNetcdfDeflateLevelFlag(9)
# Hi compression: 1.4Gb file ; # Single salt variable
# No compression: 5.6Gb ; Standard (compression/shuffling): 1.5Gb ; Hi compression w/ shuffling: 1.5Gb
cdm.setAutoBounds(1) ; # Ensure bounds on time and depth axes are generated
##

## Define useful functions ##

def clearAll():
    """
    Documentation for clearall():
    -------
    The clearall() function purges all variables in global namespace
    
    Author: Paul J. Durack : pauldurack@llnl.gov
    
    Usage:
    ------
        >>> from durolib import clearall
        >>> clearall()
    
    Notes:
    -----
        Currently not working ...
    """
    for uniquevariable in [variable for variable in globals().copy() if variable[0] != "_" and variable != 'clearall']:
        del globals()[uniquevariable]


def environment():
    return False


def fillHoles(var):
    return var
    #http://tcl-nap.cvs.sourceforge.net/viewvc/tcl-nap/tcl-nap/library/nap_function_lib.tcl?revision=1.56&view=markup
    #http://tcl-nap.cvs.sourceforge.net/viewvc/tcl-nap/tcl-nap/library/stat.tcl?revision=1.29&view=markup
    #http://stackoverflow.com/questions/5551286/filling-gaps-in-a-numpy-array
    #http://stackoverflow.com/questions/3662361/fill-in-missing-values-with-nearest-neighbour-in-python-numpy-masked-arrays
    #https://www.google.com/search?q=python+nearest+neighbor+fill
    """
     # fill_holes --
329 	#
330 	# Replace missing values by estimates based on means of neighbours
331 	#
332 	# Usage:
333 	# fill_holes(x, max_nloops)
334 	# where:
335 	# - x is array to be filled
336 	# - max_nloops is max. no. iterations (Default is to keep going until
337 	# there are no missing values)
338 	
339 	proc fill_holes {
340 	x
341 	{max_nloops -1}
342 	} {
343 	set max_nloops [[nap "max_nloops"]]
344 	set n [$x nels]
345 	set n_present 0; # ensure at least one loop
346 	for {set nloops 0} {$n_present < $n && $nloops != $max_nloops} {incr nloops} {
347 	nap "ip = count(x, 0)"; # Is present? (0 = missing, 1 = present)
348 	set n_present [[nap "sum_elements(ip)"]]
349 	if {$n_present == 0} {
350 	error "fill_holes: All elements are missing"
351 	} elseif {$n_present < $n} {
352 	nap "x = ip ? x : moving_average(x, 3, -1)"
353 	}
354 	}
355 	nap "x"
356 	}
    """



def fitPolynomial(var,time,polyOrder):
    """
    Documentation for fitPolynomial(var):
    -------
    The fitPolynomial(var,time,polyOrder) function returns a new variable which is the polyOrder
    estimate of the variable argument
    
    Author: Paul J. Durack : pauldurack@llnl.gov
    
    Usage:
    ------
        >>> from durolib import fitPolynomial
        >>> var_cubic = fitPolynomial(var,time,polyOrder=3)
    
    Notes:
    -----
    - PJD  5 Aug 2013 - Implemented following examples from Pete G.
    - TODO: only works on 2D arrays, improve to work on 3D
    http://docs.scipy.org/doc/numpy/reference/generated/numpy.polyfit.html
    """
    if polyOrder > 3:
        print "".join(['** fitPolynomial Error: >cubic fits not supported **',])
        return
    varFitted = mv.multiply(var,0.) ; # Preallocate output array    
    coefs,residuals,rank,singularValues,rcond = np.polyfit(time,var,polyOrder,full=True)
    for timeIndex in range(len(time)):
        timeVal = time[timeIndex]
        if polyOrder == 1:
            varFitted[timeIndex] = (coefs[0]*timeVal + coefs[1])
        elif polyOrder == 2:
            varFitted[timeIndex] = (coefs[0]*(timeVal**2) + coefs[1]*timeVal + coefs[2])
        elif polyOrder == 3:
            varFitted[timeIndex] = (coefs[0]*(timeVal**3) + coefs[1]*(timeVal**2) + coefs[2]*timeVal + coefs[3])
    return varFitted
    

def fixInterpAxis(var):
    """
    Documentation for fixInterpAxis(var):
    -------
    The fixInterpAxis(var) function corrects temporal axis so that genutil.statistics.linearregression
    returns coefficients which are unscaled by the time axis
    
    Author: Paul J. Durack : pauldurack@llnl.gov
    
    Usage:
    ------
        >>> from durolib import fixInterpAxis
        >>> (slope),(slope_err) = linearregression(fixInterpAxis(var),error=1,nointercept=1)
    
    Notes:
    -----
        ...
    """
    tind = range(shape(var)[0]) ; # Assume time axis is dimension 0
    t = cdm.createAxis(tind,id='time')
    t.units = 'years since 0-01-01 0:0:0.0'
    t.calendar = var.getTime().calendar
    cdu.times.setTimeBoundsYearly(t) ; # Explicitly set time bounds to yearly
    var.setAxis(0,t)
    return var


def fixVarUnits(var,varName,report=False,logFile=None):
    """
    Documentation for fixVarUnits():
    -------
    The fixVarUnits() function corrects units of salinity and converts thetao from K to degrees_C
    
    Author: Paul J. Durack : pauldurack@llnl.gov
    
    Usage:
    ------
        >>> from durolib import fixVarUnits
        >>> [var,var_fixed] = fixVarUnits(var,'so',True,'logfile.txt')
    
    Notes:
    -----
        ...
    """
    var_fixed = False
    if varName in ['so','sos']:
        if var.max() < 1. and var.mean() < 1.:
            if report:
                print "".join(['*SO mean:     {:+06.2f}'.format(var.mean()),'; min: {:+06.2f}'.format(var.min().astype('float64')),'; max: {:+06.2f}'.format(var.max().astype('float64'))])
            if logFile is not None:
                writeToLog(logFile,"".join(['*SO mean:     {:+06.2f}'.format(var.mean()),'; min: {:+06.2f}'.format(var.min().astype('float64')),'; max: {:+06.2f}'.format(var.max().astype('float64'))]))
            var_ = var*1000
            var_.id = var.id
            var_.name = var.id
            for k in var.attributes.keys():
                setattr(var_,k,var.attributes[k])
            var = var_
            var_fixed = True
            if report:
                print "".join(['*SO mean:     {:+06.2f}'.format(var.mean()),'; min: {:+06.2f}'.format(var.min().astype('float64')),'; max: {:+06.2f}'.format(var.max().astype('float64'))])
            if logFile is not None:
                writeToLog(logFile,"".join(['*SO mean:     {:+06.2f}'.format(var.mean()),'; min: {:+06.2f}'.format(var.min().astype('float64')),'; max: {:+06.2f}'.format(var.max().astype('float64'))]))
    elif varName in 'thetao':
        if var.max() > 50. and var.mean() > 265.:
            if report:
                print "".join(['*THETAO mean: {:+06.2f}'.format(var.mean()),'; min: {:+06.2f}'.format(var.min().astype('float64')),'; max: {:+06.2f}'.format(var.max().astype('float64'))])
            if logFile is not None:
                writeToLog(logFile,"".join(['*THETAO mean: {:+06.2f}'.format(var.mean()),'; min: {:+06.2f}'.format(var.min().astype('float64')),'; max: {:+06.2f}'.format(var.max().astype('float64'))]))
            var_ = var-273.15
            var_.id = var.id
            var_.name = var.id
            for k in var.attributes.keys():
                setattr(var_,k,var.attributes[k])
            var = var_
            var_fixed = True
            if report:
                print "".join(['*THETAO mean: {:+06.2f}'.format(var.mean()),'; min: {:+06.2f}'.format(var.min().astype('float64')),'; max: {:+06.2f}'.format(var.max().astype('float64'))])
            if logFile is not None:
                writeToLog(logFile,"".join(['*THETAO mean: {:+06.2f}'.format(var.mean()),'; min: {:+06.2f}'.format(var.min().astype('float64')),'; max: {:+06.2f}'.format(var.max().astype('float64'))]))

    return var,var_fixed


def getGitInfo(filePath):
    """
    Documentation for getGitInfo():
    -------
    The getGitInfo() function retrieves latest commit info specified by filePath
    
    Author: Paul J. Durack : pauldurack@llnl.gov
    
    Returns:
    -------
           gitTag[0] - commit hash
           gitTag[1] - commit author
           gitTag[2] - commit date and time
           gitTag[3] - commit notes
           
    Usage: 
    ------
        >>> from durolib import getGitInfo
        >>> gitTag = getGitInfo(filePath)
    
    Where filePath is a file which is monitored by git
            
    Notes:
    -----
        When ...
    """
    p = subprocess.Popen(['git','log','-n1','--',filePath],stdout=subprocess.PIPE,stderr=subprocess.PIPE,cwd='/'.join(filePath.split('/')[0:-1]))
    if 'fatal: Not a git repository' in p.stderr.read():
        print 'filePath not a valid git-tracked file'
        return
    gitTagFull = p.stdout.read() ; # git full tag
    del(filePath,p)
    gitTag = []
    for count,gitStr in enumerate(gitTagFull.split('\n')):
        if gitStr == '':
            pass
        else:
            gitStr = replace(gitStr,'   ',' ') ; # Trim excess whitespace in date
            gitTag.extend(["".join(gitStr.strip())])
    
    return gitTag


def globalAttWrite(file_handle,options):
    """
    Documentation for globalAttWrite():
    -------
    The globalAttWrite() function writes standard global_attributes to an
    open netcdf specified by file_handle
    
    Author: Paul J. Durack : pauldurack@llnl.gov
    
    Returns:
    -------
           Nothing.
    Usage: 
    ------
        >>> from durolib import globalAttWrite
        >>> globalAttWrite(file_handle)
    
    Where file_handle is a handle to an open, writeable netcdf file
    
    Optional Arguments:
    -------------------
    option=optionalArguments   
    Restrictions: option has to be a string
    Default : ...
    
    You can pass option='SOMETHING', ...
    
    Examples:
    ---------
        >>> from durolib import globalAttWrite
        >>> f = cdms2.open('data_file_name','w')
        >>> globalAttWrite(f)
        # Writes standard global attributes to the netcdf file specified by file_handle
            
    Notes:
    -----
        When ...
    """
    # Create timestamp, corrected to UTC for history
    local                       = pytz.timezone("America/Los_Angeles")
    time_now                    = datetime.datetime.now();
    local_time_now              = time_now.replace(tzinfo = local)
    utc_time_now                = local_time_now.astimezone(pytz.utc)
    time_format                 = utc_time_now.strftime("%d-%m-%Y %H:%M:%S %p")
    file_handle.institution     = "Program for Climate Model Diagnosis and Intercomparison (LLNL)"
    file_handle.data_contact    = "Paul J. Durack; pauldurack@llnl.gov; +1 925 422 5208"
    file_handle.history         = "".join(['File processed: ',time_format,' UTC; San Francisco, CA, USA'])
    file_handle.host            = "".join([gethostname(),'; UVCDAT version: ',".".join(["%s" % el for el in cdat_info.version()]),
                                           '; Python version: ',replace(replace(sys.version,'\n','; '),') ;',');')])

def inpaint(array,method):
    #/work/durack1/csiro/Backup/110808/Z_dur041_linux/bin/inpaint_nans/inpaint_nans.m
    return False


def keyboard(banner=None):
    """
    Documentation for keyboard():
    -------
    The keyboard() function mimics matlab's keyboard function allowing control
    sent to the keyboard within a running script
    
    Author: Paul J. Durack : pauldurack@llnl.gov
    
    Returns:
    -------
           Nothing.
    Usage: 
    ------
        >>> from durolib import keyboard
        >>> keyboard()
    
    Examples:
    ---------
        ...
            
    Notes:
    -----
        ...
    """    
    # use exception trick to pick up the current frame
    try:
        raise None
    except:
        frame = sys.exc_info()[2].tb_frame.f_back
    print "# Use quit() to exit :) Happy debugging!"
    # evaluate commands in current namespace
    namespace = frame.f_globals.copy()
    namespace.update(frame.f_locals)
    try:
        code.interact(banner=banner,local=namespace)
    except SystemExit:
        return
        
def mkDirNoOSErr(newdir,mode=0777):
    """
    Documentation for mkDirNoOSErr(newdir,mode=0777):
    -------
    The mkDirNoOSErr() function mimics os.makedirs however does not fail if the directory already
    exists
    
    Author: Paul J. Durack : pauldurack@llnl.gov
    
    Returns:
    -------
           Nothing.
    Usage: 
    ------
        >>> from durolib import mkDirNoOSErr
        >>> mkDirNoOSErr('newPath',mode=0777)
            
    Notes:
    -----
    """
    try:
        os.makedirs(newdir,mode)
    except OSError as err:
        #Re-raise the error unless it's about an already existing directory
        if err.errno != errno.EEXIST or not os.path.isdir(newdir):
            raise
    
    
def outerLocals(depth=0):
    return inspect.getouterframes(inspect.currentframe())[depth+1][0].f_locals


def smooth(array,method):
    #/apps/MATLAB/R2011b/toolbox/matlab/specgraph/smooth3.m
    #/apps/MATLAB/R2011b/toolbox/curvefit/curvefit/smooth.m
    return False


def spyderClean():
    """
    Documentation for spyder_clean():
    -------
    The spyder_clean() function purges variables initialised upon startup
    
    Author: Paul J. Durack : pauldurack@llnl.gov
    
    Usage: 
    ------
        >>> from durolib import spyder_clean
        >>> spyder_clean()
    
    Notes:
    -----
        Currently not working ...
    """
    local_vars = outerLocals()
    if 'e' in local_vars:
        print 'yep..'
        del(e,pi,sctypeNA,typeNA)
        gc.collect()


def sysCallTimeout(cmd,timeout):
    """
    Documentation for sysCallTimeout(cmd,timeout):
    -------
    The sysCallTimeout(cmd,timeout) function attempts to execute a system call (cmd) and times out in
    a specified time
    
    Author: Paul J. Durack : pauldurack@llnl.gov
    
    Usage: 
    ------
        >>> from durolib import sysCallTimeout
        >>> sysCallTimeout(cmd,timeout)
    
    Notes:
    -----
    """
    start = time.time()
    p = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    while time.time() - start < timeout:
        if p.poll() is not None:
            return
        time.sleep(0.1)
    p.kill()
    raise OSError('sysCallTimeout: System call timed out')


def trimModelList(modelFileList):
    """
    Documentation for trimModelList(modelFileList):
    -------
    The trimModelList(modelFileList) function takes a python list of model files
    and trims these for duplicates using file creation_date attribute along with
    temporal ordering info obtained from the file version identifier
    
    Author: Paul J. Durack : pauldurack@llnl.gov
    
    Usage:
    ------
        >>> modelFileList = glob.glob(os.path.join(filePath,'*.nc')) ; # provides full directory/file path        
        >>> from durolib import trimModelList
        >>> modelFileListTrimmed = trimModelList(modelFileList)
    
    Notes:
    -----
    - PJD  1 Apr 2014 - Implement sanity checks for r1i1p1 matching for e.g.
    - PJD  1 Apr 2014 - Removed hard-coded ver- position
    - PJD  1 Apr 2014 - Added realisation test to ensure expected format
    """
    # Check for list variable
    if type(modelFileList) is not list:
        print '** Function argument not type list, exiting.. **'
        return ''
    
    # Sort list and declare output
    modelFileList.sort()
    modelFileListTmp = []
    modelFileIndex = []
    
    # Create subset modelFileList
    for file1 in modelFileList:
        file1   = file1.split('/')[-1]
        mod     = file1.split('.')[1]
        exp     = file1.split('.')[2]
        rea     = file1.split('.')[3]
        # Test rea for r1i1p111 format match
        reaTest = re.compile('^r\d{1,2}i\d{1,2}p\d{1,3}')
        if not reaTest.match(rea):
            print '** Filename format invalid - rea: ',rea,', exiting.. **'
            return ''            
        modelFileListTmp.append('.'.join([mod,exp,rea]))
        
    # Create unique list and index
    modelFileListTmpUnique = list(set(modelFileListTmp)) ; modelFileListTmpUnique.sort()
    findMatches = lambda searchList,elem: [[i for i, x in enumerate(searchList) if x == e] for e in elem]
    modelFileListTmpIndex = findMatches(modelFileListTmp,modelFileListTmpUnique)        
        
    # Loop through unique list
    for count,modelNum in enumerate(modelFileListTmpIndex):
        if len(modelFileListTmpIndex[count]) == 1: # Case single version            
            modelFileIndex.append(int(str(modelNum).strip('[]')))
        else: # Case multiple versions
            # Get version and creation_date info from file
            modelFileListVersion = [] ; modelFileListCreationDate = [] ; modelFileListIndex = []
            for index in modelFileListTmpIndex[count]:
                file1 = modelFileList[index].split('/')[-1]
                verInd = int(str([count for count,x in enumerate(file1.split('.')) if 'ver-' in x]).strip('[]'))
                ver1 = file1.split('.')[verInd].replace('ver-','')                
                f_h = cdm.open(modelFileList[index])
                CD = f_h.creation_date
                f_h.close()
                modelFileListVersion.append(ver1)
                modelFileListCreationDate.append(CD)
                modelFileListIndex.append(index)
            #print modelFileListVersion
            #print modelFileListCreationDate
            
            # Use creation_date to determine latest file
            listLen = len(modelFileListCreationDate)
            modelFileListCreationDate = map(string.replace,modelFileListCreationDate,['T',]*listLen, [' ',]*listLen)
            modelFileListCreationDate = map(cdtime.s2c,modelFileListCreationDate)
            modelFileListCreationDate = map(cdtime.c2r,modelFileListCreationDate,['days since 1-1-1',]*listLen)
            modelFileListCreationDate = [x.value for x in modelFileListCreationDate]
            maxes = [i for i,x in enumerate(modelFileListCreationDate) if x == max(modelFileListCreationDate)]
            ver = [modelFileListVersion[i] for i in maxes]
            ind = [modelFileListIndex[i] for i in maxes]
            #print modelFileListCreationDate
            #print maxes,ver,ind
            
            # If creation_dates match check version info to determine latest file
            indTest = '-' ; #verTest = '-'
            if len(maxes) > 1:
                pubTest = 0 ; dateTest = 0;
                for count,ver1 in reversed(list(enumerate(ver))):
                    # Take datestamp versioned data
                    if 'v' in ver1 and ver1 > dateTest:
                        #verTest = ver[count]
                        indTest = ind[count]
                        dateTest = ver1
                    # Use published data preferentially: 1,2,3,4, ...
                    if ver1.isdigit() and ver1 > pubTest:
                        indTest = ind[count]
                        pubTest = ver1
                modelFileIndex.append(int(str(indTest).strip('[]')))
            else:
                modelFileIndex.append(int(str(ind).strip('[]')))
            #print pubTest,dateTest

    # Trim original list with new index
    modelFileListTrimmed = [modelFileList[i] for i in modelFileIndex]
        
    #return modelFileListTrimmed,modelFileIndex,modelFileListTmp,modelFileListTmpUnique,modelFileListTmpIndex ; # Debugging
    return modelFileListTrimmed
    
    
def writeToLog(logFilePath,textToWrite):
    """
    Documentation for writeToLog(logFilePath,textToWrite):
    -------
    The writeToLog() function writes specified text to a text log file
    
    Author: Paul J. Durack : pauldurack@llnl.gov
    
    Usage: 
    ------
        >>> from durolib import writeToLog
        >>> writeToLog(~/somefile.txt,'text to write to log file')
    
    Notes:
    -----
        Current version appends a new line after each call to the function.
        File will be created if it doesn't already exist, otherwise new text
        will be appended to an existing log file.
    """
    if os.path.isfile(logFilePath):
        logHandle = open(logFilePath,'a') ; # Open to append
    else:
        logHandle = open(logFilePath,'w') ; # Open to write     
    logHandle.write("".join([textToWrite,'\n']))
    logHandle.close()
    
    
def writePacked(var,fileObject='tmp.nc'):
    """
    Documentation for writePacked(var,fileObject):
    -------
    The writePacked() function generates a 16-bit (int16) cdms2 variable and
    writes this to a netcdf file
    
    Author: Paul J. Durack : pauldurack@llnl.gov
    
    Usage: 
    ------
        >>> from durolib import writePacked
        >>> writePacked(var,'16bitPacked.nc')
    
    Notes:
    -----
        TODO: clean up fileObject existence..
        TODO: deal with incredibly slow write-times
        TODO: deal with input data precision
    """
    #varType             = var.dtype
    varMin              = var.min()
    varMax              = var.max()
    var.scale_factor    = np.float32((varMax-varMin)/(2**16)-1)
    var.add_offset      = np.float32(varMin+var.scale_factor*(2**15))
    if 'tmp.nc' in fileObject:
        fileObject = cdm.open(fileObject,'w')
    fileObject.write((var-var.add_offset)/var.scale_factor,dtype=np.int16)
    return

##

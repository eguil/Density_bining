import os,glob,sys,resource,socket
from libDensityPostpro import mmeAveMsk1D,mmeAveMsk2D, mmeAveMsk3D
from modelsDef import defModels
from correctBinFiles import correctFile
from string import replace
import warnings
import time as timc
warnings.filterwarnings("ignore")
# ----------------------------------------------------------------------------
#
# Perform model ensemble mean and other statistics for density binning output
# run with 'pythoncd mme_ave_msk.py' (cdms python)
#
# April 2016 : add ToE computation support (for 2D files only)
# May 2016   : add obs support
# Nov 2016   : add 3D files support
# Jan 2017   : add picontrol and 1pctCo2 support
#
# TODO : add arguments to proc for INIT part (exper, raw, fullTS, test, keepfiles, oneD/twoD, mm/mme, ToE...) or per step
#
# ----------------------------------------------------------------------------
tcpu0 = timc.clock()
#
#  ----------------------------
# !!! Compulsory work order !!!
#  ----------------------------
# 0) create ptopsigmaxy and correct grid interpolation issues (hist and histNat)
#   0.1) raw, oneD, mm, fullTS = T, correctF = F
#   0.2) for file in cmip5.* ; do  ncks -A -v ptopsigmaxy $file ../$file ; echo $file; done
#   0.3) raw, oneD, mm, fullTS = F, correctF = T
# 1) run oneD first (mm and mme) for historical and histNat
# 2) run twoD mm for histNat
# 3) run twoD + ToE mm for historical
# 4) run twoD mme for historical (still to implement for ToE)
#
# ===============================================================================================================
#                                        INIT - work definition
# ===============================================================================================================
#raw = True
raw = False
# fullTS = True # to compute for the full range of time (used for raw/oneD to compute ptopsigmaxy)
fullTS = False
#testOneModel = True
testOneModel = False

# Initial correction of Raw binned files (longitude interpolation and bowl issues)
correctF = False  # only active if Raw = True

# Keep existing files or replace (if True and file present, ignores the model mm or mme computation)
# Use False for testing
keepFiles = False

oneD = False
twoD = False

#oneD = True
twoD = True
mme  = False
mm = True
# experiment
#exper = 'historical'
#exper = 'historicalNat'
#exper = 'piControl'
exper = '1pctCO2'
#exper = 'obs'

# Time mean/max bowl calculation used to mask out bowl
timeBowl = 'max'

if twoD:
    correctF = False # already done for oneD
# ToE
#ToE = True
ToE = False
ToeType = 'histnat'    # working from hist and histnat
#ToeType = 'picontrol' # working from hist and picontrol
if not ToE:
    ToeType ='F'

# ===============================================================================================================

hostname = socket.gethostname()
if 'locean-ipsl.upmc.fr' in hostname:
    baseDir = '/Volumes/hciclad/data/Density_binning/'
elif 'waippo.local' in hostname or 'canalip.upmc.fr' in hostname:
    if raw:
        baseDir = '/Volumes/hciclad/data/Density_binning/'
    else:
        baseDir ='/Users/ericg/Projets/Density_bining/'
elif 'private.ipsl.fr' in hostname:
    baseDir = '/data/ericglod/Density_binning/'
elif 'crunchy.llnl.gov' in hostname:
    baseDir = '/work/guilyardi/'
else:
    print hostname
    sys.exit('Unknown hostname')

if exper <> 'obs':
    # define all models
    models = defModels()

    # Years interval for difference reference
    iniyear = 1861
    peri1 = (1861-iniyear)+1
    peri2 = (1950-iniyear)+2

    # I/O directories
    #rootDir = '/Users/ericg/Projets/Density_bining/Prod_density_april15/'
    #rootDir = '/Volumes/hciclad/data/Density_binning/Prod_density_april15/Raw/'
    #rootDir = '/data/ericglod/Density_binning/Prod_density_april15/Raw/'
    #rootdir = '/work/guilyardi/Prod_density_april15/Raw'
    if raw:
        rootDir =baseDir+'Prod_density_april15/Raw/'
    else:
        rootDir =baseDir+'Prod_density_april15/'
    histDir    = rootDir+'historical'
    histNatDir = rootDir+'historicalNat'
    piControlDir = rootDir+'piControl'
    pctCO2Dir = rootDir+'1pctCO2'
    histMMEOut = rootDir+'mme_hist'
    histNatMMEOut = rootDir+'mme_histNat'
    picMMEOut = rootDir+'mme_piControl'
    pctMMEOut = rootDir+'mme_1pctCO2'
    ToeNatOut = rootDir+'toe_histNat'

    # output name
    outroot = 'cmip5.multimodel'
    inroot  = 'cmip5'
else:

# Specific variables for observations
    obsm = {'name':'EN4'     ,'props':[1,0,0,114], 'picontrol':[0]}
    #obsm = {'name':'Ishii'   ,'props':[1,0,0,67], 'picontrol':[0]}
    models = [obsm]
    if models[0]['name'] == 'EN4': # 1900.01 - 2015.04 (115 time steps, ignore last year) Good et al.
        iniyear = 1900
        peri1 = (2014-iniyear)+1
        peri2 = (1900-iniyear)+2
        idxtime = [0,114]
    elif models[0]['name']  == 'Ishii': # 1945.01 - 2012.12 (68 time steps)
        iniyear = 1945
        peri1 = (2012-iniyear)+1
        peri2 = (1945-iniyear)+2
        idxtime = [0,67]
    #rootDir = '/Users/ericg/Projets/Density_bining/Prod_density_obs_april16/'
    rootDir ='/Volumes/hciclad/data/Density_binning/Prod_density_obs_april16/'
    ObsMMEOut = rootDir+'mme_obs'
    outroot = models[0]['name']
    inroot = 'obs'
    mm = True
    mme = False
#
nmodels = len(models)

# perform a selection of a few models (for testing or updating)?
modelSel = range(nmodels)
# modelSel = [3,10,18,19,25,27,28]
#modelSel = [22,23]
if testOneModel:
    modelSel = [0]

# Select range of MME
selMME = 'All' # select all models for MME
#selMME = 'Nat' # select only models for which there are hist AND histNat simulations
#selMME = '1pct' # select only models for which there are piControl AND 1pctCO2 simulations

if mme:
    fullTS = False
    correctF = False

if ToE:
    if ToeType == 'histnat':
        selMME = 'Nat'        # force if ToE & histnat used

if exper == 'historical':
    indir  = [histDir]
    outdir = histMMEOut
    idxtime=[0,145]
elif exper == 'historicalNat':
    indir  = [histNatDir]
    outdir = histNatMMEOut
    idxtime=[0,145]
elif exper == 'piControl':
    indir  = [piControlDir]
    outdir = picMMEOut
    idxtime=[0,-140] # last 140 years are used for mme
    selMME = '1pct' # select on runs that also have a 1pctCO2
elif exper == '1pctCO2':
    indir  = [pctCO2Dir]
    outdir = pctMMEOut
    idxtime=[0,140]
    selMME = 'piCtl' # select on runs that also have a piControl
elif exper == 'obs':
    indir  = [rootDir]
    outdir = ObsMMEOut


if ToE:
    if ToeType == 'histnat':
        indir  = [histDir, histNatMMEOut]
        outdir  = ToeNatOut
if raw:
    dim = 2
    appendDim1d='2D'
    appendDim2d='3D'
    if mme:
        if exper == 'historical':
            indir = [rootDir+'mme_hist']
            outdir = rootDir+'mme_hist'
    if mme:
        if exper == 'historicalNat':
            indir = [rootDir+'mme_histNat']
            outdir = rootDir+'mme_histNat'
else:
    dim = 1
    appendDim1d='zon1D'
    appendDim2d='zon2D'

if raw & twoD :
    outdir = outdir+'/mme'
    if mme:
        indir[0] = indir[0]+'/mme'

if mme:
    indir[0]  = outdir


timeInt=[peri1,peri2]

listens = []
listens1 = []
print
print '-----------------------------------------------------------------------------------------------'
print ' Enter mme_ave_mask.py for multi-model ensemble averaging for density bins'
print '-----------------------------------------------------------------------------------------------'
if oneD:
    print ' -> work on 1D files'
if twoD:
    print ' -> work on 2D files'
if raw:
    print ' -> work on raw 4D data'
    if correctF:
        print ' -> Correct files for longitude and bowl issues'
if ToE:
    print ' -> computing ToE for type = ',ToeType
if mm:
        print ' -> Performing ensemble(s) for',exper
        print ' -> Type of time selection on bowl (mean or max):',timeBowl
if mme:
        print ' -> Performing MME for',selMME, 'models for', exper
print
print '  --> indir = ',indir
print '  --> outdir = ',outdir
print '-----------------------------------------------------------------------------------------------'
print

os.chdir(indir[0])
for i in modelSel:
    mod = models[i]['name']
    years = [models[i]['props'][3],models[i]['props'][4]]
    if exper == 'historical':
        nens = models[i]['props'][0]
        chartest = exper
    elif exper == 'historicalNat':
        nens = models[i]['props'][1]
        chartest = exper
    elif exper == 'piControl':
        nyears = models[i]['picontrol'][0]
        nens = 1
        years=[0,nyears]
        if selMME == '1pct' and nyears < 140 and nyears > 0:
            nens = 1
            print ' TOO SHORT: IGNORE model', mod
        chartest = exper
    elif exper == '1pctCO2':
        nens = models[i]['props'][2]
        years=[0,140]
        chartest = exper
    elif exper == 'obs':
        nens = models[i]['props'][0]
        chartest = 'historical'
    if ToE:
        if ToeType == 'histnat':
            nens = models[i]['props'][1]
    if years[1] <> 0: # do not ignore model
        if nens > 0: # only if 1 member or more
            if raw:
                listf  = glob.glob(inroot+'.'+mod+'.*.nc')
                listf1 = listf
            else:
                listf  = glob.glob(inroot+'.'+mod+'.*zon2D*')
                listf1 = glob.glob(inroot+'.'+mod+'.*zon1D*')

            if len(listf) == 0:
                sys.exit('### no such file !')
            start = listf[0].find(chartest)+len(chartest)
            end = listf[0].find('.an.')
            rip = listf[0][start:end]
            if raw:
                outFile = replace(listf[0],rip,'.ensm')
                outFile1 = outFile
                if mm & correctF: # correct file in indir+'/correct' and change indir
                    idxcorr = models[i]['correctFile']
                    outDirc = indir[0]+'/correct'
                    print ' -> correct',len(listf),'files towards', outDirc
                    for filec in listf:
                        # test if file is here before creating
                        if os.path.isfile(outDirc+'/'+filec):
                            print ' -> corrected file present: ',filec
                        else:
                            print ' -> correct ',filec
                            correctFile(idxcorr, 1, filec, indir[0], filec, outDirc)
                    i#ndirnew = outDirc
            else:
                outFile = replace(listf[0],rip,'.ensm')
                outFile1 = replace(outFile,'2D','1D')
            # Create lists for mme
            if mme:
                if selMME == 'All':
                    listens.append(outFile)
                    listens1.append(outFile1)
                    print ' Add ',i,mod, '(slice', years, nens, 'members) to MME'

                if selMME == 'Nat': # only select model if histNat mm is present
                    if models[i]['props'][1] > 0:
                        listens.append(outFile)
                        listens1.append(outFile1)
                        print ' Add ',i,mod, '(slice', years, nens, 'members) to MME'
                if selMME == '1pct': # only select model if 1pctCO2 mm is present
                    if models[i]['props'][2] > 0:
                        listens.append(outFile)
                        listens1.append(outFile1)
                        print ' Add ',i,mod, '(slice', years, nens, 'members) to MME'
                if selMME == 'piCtl': # only select model if piCtl mm is present
                    if models[i]['picontrol'][0] > 0:
                        listens.append(outFile)
                        listens1.append(outFile1)
                        print ' Add ',i,mod, '(slice', years, nens, 'members) to MME'
            # Perform model ensemble
            if mm:
                if twoD:
                    if os.path.isfile(outdir+'/'+outFile) & keepFiles:
                        print ' -> File exists - IGNORE mm of',outFile,'already in',outdir
                    else:
                        print ' -> working on: ', i,mod, 'slice', years, nens, 'members'
                        if dim == 1:
                            mmeAveMsk2D(listf,years,indir,outdir,outFile,timeInt,mme,timeBowl,ToeType)
                        elif dim == 2:
                            mmeAveMsk3D(listf,years,indir,outdir,outFile,timeInt,mme,ToeType)
                        print 'Wrote ',outdir+'/'+outFile
                if oneD:
                    if os.path.isfile(outdir+'/'+outFile1) & keepFiles:
                        print ' -> File exists - IGNORE mm of',outFile1,'already in',outdir
                    else:
                        print ' -> working on: ', i,mod, 'slice', years, nens, 'member(s)'
                        mmeAveMsk1D(listf1,dim,years,indir,outdir,outFile1,timeInt,mme,ToeType,fullTS)
                        print 'Wrote ',outdir+'/'+outFile1
                    
if mme:
    # run 1D MME first
    if twoD:
        outFile = outroot+'_'+selMME+'.'+exper+'.ensm.an.ocn.Omon.density_'+appendDim2d+'.nc'
        if os.path.isfile(outdir+'/'+outFile) & keepFiles:
            print ' -> IGNORE: mme of',outFile,'already in',outdir
        else:
            if dim == 1:
                mmeAveMsk2D(listens,idxtime,indir,outdir,outFile,timeInt,mme,timeBowl,ToeType)
            elif dim ==2:
                mmeAveMsk3D(listens,idxtime,indir,outdir,outFile,timeInt,mme,ToeType)
            print 'Wrote ',outdir+'/'+outFile
    if oneD:
        outFile1 = outroot+'_'+selMME+'.'+exper+'.ensm.an.ocn.Omon.density_'+appendDim1d+'.nc'
        if os.path.isfile(outdir+'/'+outFile1) & keepFiles:
            print ' -> IGNORE: mme of',outFile1,'already in',outdir
        else:
            mmeAveMsk1D(listens1,dim,idxtime,indir,outdir,outFile1,timeInt,mme,ToeType,False)
            print 'Wrote ',outdir+'/'+outFile1
tcpu1 = timc.clock()

print ' Max memory use',resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1.e6,'GB'
print ' CPU use',tcpu1-tcpu0


# ---------------------------

#modelsurf = ['ACCESS1-0','ACCESS1-3','CMCC-CESM','CMCC-CM','CMCC-CMS','CNRM-CM5','CSIRO-Mk3-6-0','EC-EARTH','FGOALS-s2','GFDL-ESM2G','GISS-E2-R-CC','GISS-E2-R','MIROC5','MIROC-ESM-CHEM','MIROC-ESM','MPI-ESM-LR','MPI-ESM-MR','MPI-ESM-P','NorESM1-ME','NorESM1-M']



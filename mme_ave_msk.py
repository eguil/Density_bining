import os,glob,sys
from libDensity import defModels,mmeAveMsk2D,mmeAveMsk1D
from string import replace
import warnings

warnings.filterwarnings("ignore")
# ----------------------------------------------------------------------------
#
# Perform model ensemble mean and other statistics for density binning output
# run with 'pythoncd -W ignore mme_ave_msk.py' (cdms python on mac)
#
# April 2016 : add ToE computation support (for 2D files only)
# Mai 2016   : add obs support
#
# ----------------------------------------------------------------------------

#
#  ----------------------------
# !!! Compulsory work order !!!
#  ----------------------------
# 1) run oneD first (mm and mme) for historical and histNat
# 2) run twoD mm for histNat
# 3) run twoD + ToE mm for historical
# 4) run twoD mme for historical (still to implement for ToE)

oneD = False
twoD = False

oneD = True
twoD = True
mme  = False
mm = True
# experiment
exper  = 'historical'
#exper  = 'historicalNat'
exper = 'obs'

# ToE
ToE = True
ToE = False
ToeType = 'histnat'    # working from hist and histnat
#ToeType = 'picontrol' # working from hist and picontrol
if not ToE:
    ToeType ='F'
if exper <> 'obs':
    # define all models
    models = defModels()

    # Years interval for difference reference
    iniyear = 1861
    peri1 = (1861-iniyear)+1
    peri2 = (1950-iniyear)+2
    idxtime=[0,145]

    # I/O directories
    rootDir = '/Users/ericg/Projets/Density_bining/Prod_density_april15/'
    histDir    = rootDir+'historical'
    histNatDir = rootDir+'historicalNat'
    histMMEOut = rootDir+'mme_hist'
    histNatMMEOut = rootDir+'mme_histNat'
    ToeNatOut = rootDir+'toe_histNat'

    # output name
    outroot = 'cmip5.multimodel'
    inroot  = 'cmip5'
else:
# Specific variables for observations
    #obsm = {'name':'EN4'     ,'props':[1,0,0,114], 'picontrol':[0]}
    obsm = {'name':'Ishii'   ,'props':[1,0,0,67], 'picontrol':[0]}
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
    rootDir = '/Users/ericg/Projets/Density_bining/Prod_density_obs_april16/'
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

# Select range of MME
selMME = 'All' # select all models for MME
#selMME = 'Nat' # select only models for which there are hist AND histNat simulations

if ToE:
    if ToeType == 'histnat':
        selMME = 'Nat'        # force if ToE & histnat used

if exper == 'historical':
    indir  = [histDir]
    outdir = histMMEOut
elif exper == 'historicalNat':
    indir  = [histNatDir]
    outdir = histNatMMEOut
elif exper == 'obs':
    indir  = [rootDir]
    outdir = ObsMMEOut

if ToE:
    if ToeType == 'histnat':
        indir  = [histDir, histNatMMEOut]
        outdir  = ToeNatOut

timeInt=[peri1,peri2]

listens = []
listens1 = []
print
print '-------------------------------------------------------------------------'
print 'Enter mme_ave_mask.py for multi-model ensemble averaging for density bins'
if oneD:
    print ' -> work on 1D files'
if twoD:
    print ' -> work on 2D files'
if ToE:
    print ' -> computing ToE for type = ',ToeType
if mm:
        print ' -> Performing ensemble(s) for',exper
if mme:
        print ' -> Performing MME for',selMME, 'models for', exper
print

os.chdir(indir[0])
for i in modelSel:
    mod = models[i]['name']
    if exper == 'historical':
        nens = models[i]['props'][0]
        chartest = exper
    elif exper == 'historicalNat':
        nens = models[i]['props'][1]
        chartest = exper
    elif exper == 'obs':
        nens = models[i]['props'][0]
        chartest = 'historical'
    if ToE:
        if ToeType == 'histnat':
            nens = models[i]['props'][1]
    years = [models[i]['props'][2],models[i]['props'][3]]
    if years[1] <> 0: # do not ignore model
        if nens > 0: # only if 1 member or more
            listf  = glob.glob(inroot+'.'+mod+'.*zon2D*')
            listf1 = glob.glob(inroot+'.'+mod+'.*zon1D*')
            if len(listf) == 0:
                sys.exit('### no such file !')
            start = listf[0].find(chartest)+len(chartest)
            end = listf[0].find('.an.')
            rip = listf[0][start:end]
            outFile = replace(listf[0],rip,'.ensm')
            outFile1 = replace(outFile,'2D','1D')
            # Create lists for mme
            if mme:
                if selMME == 'All':
                    listens.append(outFile)
                    listens1.append(outFile1)
                    print ' Add ',i,mod, 'slice', years, nens, 'members to MME'

                if selMME == 'Nat': # only select model if histNat mm is present
                    if models[i]['props'][1] > 0:
                        listens.append(outFile)
                        listens1.append(outFile1)
                        print ' Add ',i,mod, 'slice', years, nens, 'members to MME'
            # Perform model ensemble
            if mm:
                if twoD:
                    print ' -> working on: ', i,mod, 'slice', years, nens, 'members'
                    mmeAveMsk2D(listf,years,indir,outdir,outFile,timeInt,mme,ToeType)
                    print 'Wrote ',outdir+'/'+outFile
                if oneD:
                    print ' -> working on: ', i,mod, 'slice', years, nens, 'members'
                    mmeAveMsk1D(listf1,years,indir,outdir,outFile1,timeInt,mme,ToeType)
                    print 'Wrote ',outdir+'/'+outFile1
                    
if mme:
    # run 1D MME first
    indir  = outdir
    if twoD:
        outFile = outroot+'_'+selMME+'.'+exper+'.ensm.an.ocn.Omon.density_zon2D.nc'
        mmeAveMsk2D(listens,idxtime,indir,outdir,outFile,timeInt,mme,Toetype)
        print 'Wrote ',outdir+'/'+outFile
    if oneD:
        outFile1 = outroot+'_'+selMME+'.'+exper+'.ensm.an.ocn.Omon.density_zon1D.nc'
        mmeAveMsk1D(listens1,idxtime,indir,outdir,outFile1,timeInt,mme,Toetype)
        print 'Wrote ',outdir+'/'+outFile1

# ---------------------------

#modelsurf = ['ACCESS1-0','ACCESS1-3','CMCC-CESM','CMCC-CM','CMCC-CMS','CNRM-CM5','CSIRO-Mk3-6-0','EC-EARTH','FGOALS-s2','GFDL-ESM2G','GISS-E2-R-CC','GISS-E2-R','MIROC5','MIROC-ESM-CHEM','MIROC-ESM','MPI-ESM-LR','MPI-ESM-MR','MPI-ESM-P','NorESM1-ME','NorESM1-M']



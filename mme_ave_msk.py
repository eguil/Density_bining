import os,glob
from libDensity import defModels,mmeAveMsk2D,mmeAveMsk1D
from string import replace
import warnings

warnings.filterwarnings("ignore")
# ----------------------------------------------------------------------------
#
# Perform model ensemble mean and other statistics for density binning output
# run with 'pythoncd -W ignore mme_ave_msk.py' (cdms python on mac)
#
# ----------------------------------------------------------------------------

oneD = False
twoD = False

#oneD = True
twoD = True
mme  = False
mm = True

exper  = 'historicalNat'

# define all models
models = defModels()
nmodels = len(models)

# perform a selection of a few models (for testing or updating)?
modelSel = range(nmodels)
# modelSel = [3,10,18,19,25,27,28]
# modelSel = [22,23,24]

# Select range of MME
selMME = 'All' # select all models for MME
#selMME = 'Nat' # select only models for which there are hist AND histNat simulations

# Years interval for difference reference
iniyear = 1861
peri1 = (1861-iniyear)+1
peri2 = (1950-iniyear)+2
timeInt=[peri1,peri2]
idxtime=[0,145]

if exper == 'historical':
    indir  = '/Users/ericg/Projets/Density_bining/Prod_density_april15/historical'
    outdir = '/Users/ericg/Projets/Density_bining/Prod_density_april15/mme_hist'
if exper == 'historicalNat':
    indir  = '/Users/ericg/Projets/Density_bining/Prod_density_april15/historicalNat'
    outdir = '/Users/ericg/Projets/Density_bining/Prod_density_april15/mme_histNat'
listens = []
listens1 = []
print
print 'Number of models to consider in mme_ave_mask.py:',nmodels
print ' -> work: oneD, twoD, mm, mme:', oneD, twoD
if mm:
        print ' -> Performing model ensembles for',exper
if mme:
        print ' -> Performing MME for',selMME, 'models for', exper
print
os.chdir(indir)
for i in modelSel:
    mod = models[i]['name']
    if exper == 'historical':
        nens = models[i]['props'][0]
    if exper == 'historicalNat':
        nens = models[i]['props'][1]
    years = [models[i]['props'][2],models[i]['props'][3]]
    if years[1] <> 0: # do not ignore model
        if nens > 0: # only if 1 member or more
            listf  = glob.glob('cmip5.'+mod+'.*zon2D*')
            listf1 = glob.glob('cmip5.'+mod+'.*zon1D*')
            start = listf[0].find(exper)+len(exper)
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
                    print ' -> working on: ', i,mod, 'slice', years, nens
                    mmeAveMsk2D(listf,years,indir,outdir,outFile,timeInt,mme)
                    print 'Wrote ',outdir+'/'+outFile
                if oneD:
                    print ' -> working on: ', i,mod, 'slice', years, nens
                    mmeAveMsk1D(listf1,years,indir,outdir,outFile1,timeInt,mme)
                    print 'Wrote ',outdir+'/'+outFile1
                    
if mme:
    # run 1D MME first
    indir  = outdir
    if twoD:
        outFile = 'cmip5.multimodel_'+selMME+'.historical.ensm.an.ocn.Omon.density_zon2D.nc'
        mmeAveMsk2D(listens,idxtime,indir,outdir,outFile,timeInt,mme)
        print 'Wrote ',outdir+'/'+outFile
    if oneD:
        outFile1 = 'cmip5.multimodel_'+selMME+'.historical.ensm.an.ocn.Omon.density_zon1D.nc'
        mmeAveMsk1D(listens1,idxtime,indir,outdir,outFile1,timeInt,mme)
        print 'Wrote ',outdir+'/'+outFile1

# ---------------------------

#modelsurf = ['ACCESS1-0','ACCESS1-3','CMCC-CESM','CMCC-CM','CMCC-CMS','CNRM-CM5','CSIRO-Mk3-6-0','EC-EARTH','FGOALS-s2','GFDL-ESM2G','GISS-E2-R-CC','GISS-E2-R','MIROC5','MIROC-ESM-CHEM','MIROC-ESM','MPI-ESM-LR','MPI-ESM-MR','MPI-ESM-P','NorESM1-ME','NorESM1-M']



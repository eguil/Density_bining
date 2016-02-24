import gc,os,resource,timeit,glob,re,math
from libDensity import defModels,mmeAveMsk2D,mmeAveMsk1D
from string import replace
import warnings

def fxn():
    warnings.warn("deprecated", DeprecationWarning)

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    fxn()
# ----------------------------------------------------------------------------
#
# Perform model ensemble mean and other statistics for density binning output
# run with pythoncd mme_ave_msk.py (cdms python on mac)
#
# ----------------------------------------------------------------------------

#twoD = False
#oneD = True
twoD = True
oneD = False
mm  = False
mme = True 

exper  = 'historical'
# define models
models = defModels()

# Years interval for difference reference
iniyear = 1861
peri1 = (1861-iniyear)+1
peri2 = (1950-iniyear)+2
timeInt=[peri1,peri2]

#indir  = '/Users/ericg/Projets/Density_bining/Prod_density_nov14/z_individual'
#outdir = '/Users/ericg/Projets/Density_bining/Prod_density_nov14/test_mme'
indir  = '/Users/ericg/Projets/Density_bining/Prod_density_april15/historical'
outdir = '/Users/ericg/Projets/Density_bining/Prod_density_april15/mme_hist'
listens = []
listens1 = []
nmodels = len(models)
print nmodels
os.chdir(indir)
for i in range(nmodels):
    print 
    mod = models[i]['name']
    if exper == 'historical':
        nens = models[i]['props'][0]
    if exper == 'historicalNat':
        nens = models[i]['props'][1]
    years = [models[i]['props'][2],models[i]['props'][3]]
    print i,mod, 'slice', years
    if years[1] <> 0: # do not ignore model
        listf  = glob.glob('cmip5.'+mod+'.*zon2D*')
        listf1 = glob.glob('cmip5.'+mod+'.*zon1D*')
        start = listf[0].find(exper)+len(exper)
        end = listf[0].find('.an.')
        rip = listf[0][start:end]
        outFile = replace(listf[0],rip,'.ensm')
        outFile1 = replace(outFile,'2D','1D')
        listens.append(outFile)
        listens1.append(outFile1)
        if mm:
            if twoD:
                print ' -> working on: ', outFile
                mmeAveMsk2D(listf,years,indir,outdir,outFile,timeInt,mme)
                print 'Wrote ',outdir+'/'+outFile
            if oneD:
                print ' -> working on: ', outFile1
                mmeAveMsk1D(listf1,years,indir,outdir,outFile1,timeInt,mme)
                print 'Wrote ',outdir+'/'+outFile1
                    
        if mme:
            # MME
            indir  = outdir
            if twoD:
                outFile = 'cmip5.multimodel.historical.ensm.an.ocn.Omon.density_zon2D.nc'
                mmeAveMsk2D(listens,[0,146],indir,outdir,outFile,timeInt,mme)
                print 'Wrote ',outdir+'/'+outFile
            if oneD:
                outFile1 = 'cmip5.multimodel.historical.ensm.an.ocn.Omon.density_zon1D.nc'
                mmeAveMsk1D(listens1,[0,146],indir,outdir,outFile1,timeInt,mme)
                print 'Wrote ',outdir+'/'+outFile1

# ---------------------------

#modelsurf = ['ACCESS1-0','ACCESS1-3','CMCC-CESM','CMCC-CM','CMCC-CMS','CNRM-CM5','CSIRO-Mk3-6-0','EC-EARTH','FGOALS-s2','GFDL-ESM2G','GISS-E2-R-CC','GISS-E2-R','MIROC5','MIROC-ESM-CHEM','MIROC-ESM','MPI-ESM-LR','MPI-ESM-MR','MPI-ESM-P','NorESM1-ME','NorESM1-M']



import numpy as np
from netCDF4 import Dataset as open_ncfile
import os, glob
from modelsDef import defModels
import matplotlib.pyplot as plt


def read_toe_rcp85(varread, listfiles, ignore, ndomains):
    """
    Read ToE array from different model files, save in one ouput array with one dimension being total number
    of runs
    :param varread: string, ToE variable name in files
           listfiles: list of files to read (one file per model)
           ignore: list of models to ignore
           ndomains: number of domain names
    :return: varToEA: numpy array (number of total runs, number of domains) ToE in the Atlantic
             varToEP numpy array (number of total runs, number of domains)
             varToEI numpy array (number of total runs, number of domains)
             nMembers: numpy array containing number of members per model, dimension=nb of models
    """
    nruns=0
    nrunmax = 100
    nmodels = len(listfiles)
    nMembers = np.ma.zeros(nmodels)

    # Initialize varToE containing ToE of all runs
    varToEA = np.ma.masked_all((nrunmax, ndomains))
    varToEP = np.ma.masked_all((nrunmax, ndomains))
    varToEI = np.ma.masked_all((nrunmax, ndomains))

    for i in range(nmodels):
        file_toe = listfiles[i]
        ftoe = open_ncfile(file_toe, 'r')
        name = os.path.basename(file_toe).split('.')[1]

        if name not in ignore:
            # Read ToE (members, basin, domain)
            toeread = ftoe.variables[varread][:]
            nMembers[i] = int(toeread.shape[0])
            print('- Reading ToE of %s with %d members'%(name,nMembers[i]))
            nruns1 = int(nruns + nMembers[i])
            
            run_labels = ftoe.variables['run_label'][:]
            #print('   ',run_labels)

            # Save ToE
            varToEA[nruns:nruns1,:] = toeread[:,1,:]
            varToEP[nruns:nruns1,:] = toeread[:,2,:]
            varToEI[nruns:nruns1,:] = toeread[:,3,:]

            nruns = nruns1

    print('Total number of runs:', nruns)
    varToEA = varToEA[0:nruns,:]
    varToEP = varToEP[0:nruns,:]
    varToEI = varToEI[0:nruns,:]

    return varToEA, varToEP, varToEI, nMembers


def read_gsat_rcp85(indir,listfiles,ignore):

    """
    Read GSAT array from different model files, save anomaly in one ouput array with one dimension being total number
    of runs, the other being the time
    :param indir: directories where GSAT files are stored
           listfiles: list of ToE files to read (one file per model) to keep same order as when reading ToE
           ignore: list of models to ignore
    :return: GSAT anomaly array (years, nb of runs)

    """

    nruns=0
    nrunmax = 100
    nmodels = len(listfiles)
    nMembers = np.ma.empty(nmodels)
    timN=240

    # -- Initialize GSAT
    gsat_anom = np.ma.masked_all((timN,nrunmax))

    for i in range(nmodels):
        file_toe = listfiles[i] # Read toe file in the same order to retrieve model name
        name = os.path.basename(file_toe).split('.')[1]
        ftoe = open_ncfile(file_toe,'r')
        run_labels = ftoe.variables['run_label'][:] # Read run labels of model i
        if name != 'HadGEM2-ES':
            iystart = 11
        else:
            iystart = 2

        if name not in ignore:

            # Read GSAT
            file_gsat = glob.glob(indir+ 'GSAT.*'+name+'*.nc')[0] # GSAT File
            fgsat = open_ncfile(file_gsat, 'r')
            gsatread = fgsat.variables['GSAT'][:]
            run_labels_GSAT = fgsat.variables['members_name'][:]
            nMembers[i] = int(gsatread.shape[1])
            print('- Reading GSAT of %s with %d members'%(name,nMembers[i]))
            # print('  gsatread shape : ',gsatread.shape)
            nruns1 = int(nruns + nMembers[i])
            
            # Re-organize order of members so that it's the same as ToE array
            gsatread_cor = np.ma.masked_all_like(gsatread)
            for k in range(int(nMembers[i])):
                idx_cor = list(run_labels).index(run_labels_GSAT[k])
                gsatread_cor[:,idx_cor] = gsatread[:,k]
                #print('  ',k,run_labels_GSAT[k], run_labels[k], run_labels[idx_cor])
            
            # Save GSAT
            gsatread_anom = gsatread_cor - np.ma.average(gsatread_cor[0:50,:],axis=0) # Anomaly relative to first 50 years
            gsat_anom[:,nruns:nruns1] = gsatread_anom[iystart:,:] # Keep 1861-2005
            # print('  gsatanom shape : ',gsat_anom.shape)
            nruns = nruns1

    print('Total number of runs:', nruns)
    gsat_anom = gsat_anom[:,0:nruns]

    return gsat_anom

def read_toe_1pctCO2(varread, indir, models, ignore, ndomains):
    """
    Read ToE array from different model files, save in one ouput array with one dimension being total number
    of models
    :param varread: string, ToE variable name in files
           indir: directory of files to read (one file per model)
           models: dict containing models' characteristics
           ignore: list of models to ignore
           ndomains: number of domain names
    :return: varToEA: numpy array (number of total runs, number of domains) ToE in the Atlantic
             varToEP numpy array (number of total runs, number of domains)
             varToEI numpy array (number of total runs, number of domains)
             nMembers: numpy array containing number of members per model, dimension=nb of models
    """

    nmodels=len(models)

    # -- Initialize varToE containing ToE
    varToEA = np.ma.masked_all((nmodels,ndomains))
    varToEP = np.ma.masked_all((nmodels,ndomains))
    varToEI = np.ma.masked_all((nmodels,ndomains))

    nMembers=np.ma.zeros(nmodels)

    for i, model in enumerate(models):

        if model['name'] not in ignore:

            # Read file
            file_CO2piC = glob.glob(indir+'/'+ '*'+model['name']+'*.nc')[0]
            fpiC = open_ncfile(file_CO2piC, 'r')

            # Read ToE (basin, domain)
            toeread = fpiC.variables[varread][:]
            print('- Reading ToE of ' + model['name'])

            # Save ToE
            varToEA[i,:] = toeread[1,:]
            varToEP[i,:] = toeread[2,:]
            varToEI[i,:] = toeread[3,:]

            nMembers[i] = 1

    print('Total number of models:', np.sum(nMembers))

    # Remove rows with masked values (due to ignored model(s))
    if len(ignore) != 0:
        idx = np.argwhere(nMembers==0)[0]
        if len(idx) !=0:
            idx = idx[0]
            varToEA_bis = np.delete(varToEA,idx,0) # Delete row at index idx
            varToEA = np.ma.masked_greater(varToEA_bis,200) # Now turn back into a masked array
            varToEP_bis = np.delete(varToEP,idx,0)
            varToEP = np.ma.masked_greater(varToEP_bis,200)
            varToEI_bis = np.delete(varToEI,idx,0)
            varToEI = np.ma.masked_greater(varToEI_bis,200)

    return varToEA, varToEP, varToEI, nMembers


def read_gsat_1pctCO2(indir,models,ignore):

    """
    Read GSAT array from different model files, save anomaly in one ouput array with one dimension being total number
    of runs, the other being the time
    :param indir: directories where GSAT files are stored
           models: dict containing models' characteristics
           ignore: list of models to ignore
    :return: GSAT anomaly array (years, nb of runs)

    """

    nmodels = len(models)
    nMembers = np.ma.zeros(nmodels)
    timN=140

    indir_1pctCO2 = indir + '1pctCO2/'
    indir_piC = indir + 'piControl/'

    # -- Initialize GSAT
    gsat_anom = np.ma.masked_all((timN,nmodels-len(ignore)))
    imod = 0

    for i,model in enumerate(models):

        if model['name'] not in ignore:

            # Read GSAT 1pctCO2
            file_gsat_CO2 = glob.glob(indir_1pctCO2+ '*'+model['name']+'*.nc')[0] # GSAT File
            fgsatCO2 = open_ncfile(file_gsat_CO2, 'r')
            gsatread_CO2 = fgsatCO2.variables['GSAT'][0:timN]
            print('- Reading GSAT of %s'%(model['name'],))

            # Read GSAT PiControl
            file_gsat_piC = glob.glob(indir_piC+ '*'+model['name']+'*.nc')[0] # GSAT File
            fgsatpiC = open_ncfile(file_gsat_piC, 'r')
            gsatread_piC = fgsatpiC.variables['GSAT'][:]

            # Compute and Save GSAT anomaly
            gsat_anom[:,imod] = gsatread_CO2 - np.ma.average(gsatread_piC,axis=0)
            imod = imod+1

            nMembers[i]=1

    print('Total number of models:', np.sum(nMembers))

    return gsat_anom
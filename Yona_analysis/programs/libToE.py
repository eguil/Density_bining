import numpy as np
import time as timc

def findToE(signal, noise, mult):
    '''
    define Time of Emergence (ToE) from last time index at which signal is larger than mult*noise
        signal is [time,space]
        noise is [space]
        mult is float
    '''
    #tcpu0 = timc.clock()
    timN = signal.shape[0]
    toe_wrk = np.ma.ones(signal.shape)*1. # init toe_wrk array to 1
    signaltile = np.ma.reshape(np.tile(noise,timN),signal.shape) # repeat noise timN
    toe_idx = np.argwhere(abs(signal) >= mult*signaltile) # find indices of points where signal > noise
    if signal.size > timN: # if there are at least 2 dimensions
        toe_wrk[toe_idx[:,0],toe_idx[:,1]] = 0. # set corresponding points in toe_wrk to zero
    else: # if there is only the time dimension
        toe_wrk[toe_idx[:,0]] = 0
    toe = timN-np.flipud(toe_wrk).argmax(axis=0) # compute ToE as last index when signal > noise
    #tcpu1 = timc.clock()
    # perf
    #print ' ToE CPU = ',tcpu1-tcpu0

    return toe

def findToE_2thresholds(signal, noise1, noise2, tidx, mult):
    '''
    define Time of Emergence (ToE) from last time index at which signal is larger than mult*noise
        signal is [time]
        noise1 and noise2 are single-values
        mult is float
        tidx is an integer
        use noise1 over first part of the signal (until tidx), noise2 over second part
    1D case for now
    '''
    timN = signal.shape[0]
    toe_wrk = np.ma.ones(signal.shape)*1. # init toe_wrk array to 1
    signaltile1 = np.tile(noise1,tidx) # repeat noise1
    signaltile2 = np.tile(noise2,timN-tidx)
    signaltile = np.concatenate((signaltile1, signaltile2))
    toe_idx = np.argwhere(abs(signal) >= mult*signaltile) # find indices of points where signal > noise
    # if signal.size > timN: # if there are at least 2 dimensions
    #     toe_wrk[toe_idx[:,0],toe_idx[:,1]] = 0. # set corresponding points in toe_wrk to zero
    # else: # if there is only the time dimension
    toe_wrk[toe_idx[:,0]] = 0
    toe = timN-np.flipud(toe_wrk).argmax(axis=0) # compute ToE as last index when signal > noise

    return toe



def ToEdomainhistvshistNat(model_name, domain_name):

    '''
    Save domain boxes for hist vs. histNat ensemble means (salinity/temperature)

    :param model_name: name of model
    :param domain_name: Southern ST, Southern Ocean, etc...
    :return: box(es) of the specified model and domain

    '''
    if domain_name == 'Southern ST':

        # Southern Subtropics (cooling/freshening in all three basins)
        domains = [
            {'name':'bcc-csm1-1'    , 'Atlantic': [-30,-25,25.4,26.2],  'Pacific': [-20,-12,24.4,25.4], 'Indian': [-30,-20,25.4,26.1]}, # 0
            {'name':'CanESM2'       , 'Atlantic': [-18,-10,25.4,26.1],  'Pacific': [-30,-10,24.8,26],   'Indian': [-38,-20,26,26.4]}, # 1
            {'name':'CCSM4'         , 'Atlantic': [-35,-15,25,26],      'Pacific': [-40,-15,24.8,26.4], 'Indian': [-37,-20,25.6,26.2]}, # 2
            {'name':'CESM1-CAM5'    , 'Atlantic': [-20,-7,25,26.4],     'Pacific': [-40,-20,25.2,26.25],'Indian': [-35,-20,25.6,26.2]}, # 3
            {'name':'CNRM-CM5'      , 'Atlantic': [-40,-15,25.2,26.3],  'Pacific': [-25,-10,24.6,25.6], 'Indian': [-35,-15,25.4,26.2]}, # 4
            {'name':'CSIRO-Mk3-6-0' , 'Atlantic': [-18,-10,26,26.3],    'Pacific': [-30,-15,24.6,25.6], 'Indian': [-40,-15,26.3,26.6]}, #5
            {'name':'FGOALS-g2'     , 'Atlantic': [-15,-10,26.1,26.6],  'Pacific': [-25,-15,24.4,25.4], 'Indian': [-35,-20,25.4,26]}, # 6
            {'name':'GFDL-CM3'      , 'Atlantic': [-33,-25,25.8,26.4],  'Pacific': [-35,-10,24.6,26.2], 'Indian': [-37,-30,26.2,26.5]}, # 7
            {'name':'GFDL-ESM2M'    , 'Atlantic': [-30,-20,25.4,26],    'Pacific': [-28,-15,26.1,26.4], 'Indian': [-35,-20,25.6,26.2]}, # 8
            {'name':'GISS-E2-R'     , 'Atlantic': [-35,-25,25.6,26.5],  'Pacific': [-25,-10,24,25.8],   'Indian': [-30,-15,25.4,26.1]}, # 9
            {'name':'HadGEM2-ES'    , 'Atlantic': [-23,-5,26.1,27],     'Pacific': [-25,-15,24,25.2],   'Indian': [-30,-20,25.6,26]}, # 10
            {'name':'IPSL-CM5A-LR'  , 'Atlantic': [-43,-35,27.2,27.5],  'Pacific': [-25,-15,25.8,26.6], 'Indian': [-35,-20,26.5,26.9]}, # 11
            {'name':'IPSL-CM5A-MR'  , 'Atlantic': [-40,-30,26.9,27.2],  'Pacific': [-30,-15,25.4,26.6], 'Indian': [-40,-25,26.6,26.9]}, # 12
            {'name':'MIROC-ESM-CHEM', 'Atlantic': [-25,-10,25.2,26.1],  'Pacific': [-20,-8,24.6,26],    'Indian': [-37,-20,25.6,26.3]}, # 13
            {'name':'MIROC-ESM'     , 'Atlantic': [-35,-20,25.8,26.4],  'Pacific': [-35,-25,25.8,26.2], 'Indian': [-35,-20,25.6,26.4]}, # 14
            {'name':'MME'     ,       'Atlantic': [-30,-10,25.8,26.5],  'Pacific': [-25,-10,24.6,26.2], 'Indian': [-40,-15,25.8,26.5]}, # 15
        ]
        domain_char = {'nb_basins': 3, 'Atlantic': True, 'Pacific': True, 'Indian': True}

    if domain_name == 'SO':

        # Southern Ocean (warmer/saltier in all three basins)
        domains = [
            {'name':'bcc-csm1-1'    , 'Atlantic': [-60,-47,27,27.5],    'Pacific': [-63,-57,27.1,27.4],  'Indian': [-60,-50,27.1,27.5]}, # 0
            {'name':'CanESM2'       , 'Atlantic': [-48,-40,26.4,27.4],  'Pacific': [-60,-50,27.1,27.5],  'Indian': [-55,-45,27.1,27.6]}, # 1
            {'name':'CCSM4'         , 'Atlantic': [-57,-50,27.5,27.9],  'Pacific': [-57,-50,27,27.6],    'Indian': [-55,-45,27,27.6]}, # 2
            {'name':'CESM1-CAM5'    , 'Atlantic': [-55,-47,26.9,27.7],  'Pacific': [-60,-50,26.8,27.5],  'Indian': [-55,-50,26.9,27.5]}, # 3
            {'name':'CNRM-CM5'      , 'Atlantic': [-55,-45,26.5,27.1],  'Pacific': [-55,-50,26.4,27],    'Indian': [-55,-47,26.5,27.1]}, # 4
            {'name':'CSIRO-Mk3-6-0' , 'Atlantic': [-58,-40,27.2,27.9],  'Pacific': None,                 'Indian': [-60,-50,27.6,28]}, #5
            {'name':'FGOALS-g2'     , 'Atlantic': [-65,-55,26.2,27.3],  'Pacific': [-65,-60,27.6,28],    'Indian': [-60,-45,27.1,27.4]}, # 6
            {'name':'GFDL-CM3'      , 'Atlantic': [-65,-50,27.6,28],    'Pacific': [-67,-60,27.9,28],    'Indian': [-70,-60,27.7,28]}, # 7
            {'name':'GFDL-ESM2M'    , 'Atlantic': [-60,-50,27,27.4],    'Pacific': [-70,-65,27.6,27.9],  'Indian': [-60,-50,27.7,27.9]}, # 8
            {'name':'GISS-E2-R'     , 'Atlantic': [-55,-45,27,27.5],    'Pacific': [-57,-50,26.9,27.4],  'Indian': [-53,-45,27,27.3]}, # 9
            {'name':'HadGEM2-ES'    , 'Atlantic': [-55,-47,27,27.3],    'Pacific': [-60,-58,27.2,27.6],  'Indian': [-55,-50,27,27.3]}, # 10
            {'name':'IPSL-CM5A-LR'  , 'Atlantic': [-60,-50,27.7,27.9],  'Pacific': [-65,-60,27.7,27.9],  'Indian': [-55,-50,27.7,27.8]}, # 11
            {'name':'IPSL-CM5A-MR'  , 'Atlantic': [-55,-50,27.6,27.9],  'Pacific': [-65,-60,27.5,27.8],  'Indian': [-55,-50,27.6,27.8]}, # 12
            {'name':'MIROC-ESM-CHEM', 'Atlantic': [-55,-50,27.2,27.8],  'Pacific': [-60,-55,27.4,27.6],  'Indian': [-55,-50,27.3,27.7]}, # 13
            {'name':'MIROC-ESM'     , 'Atlantic': [-55,-50,27.4,27.9],  'Pacific': None,                 'Indian': [-55,-50,27.3,27.8]}, # 14
            {'name':'MME'           , 'Atlantic': [-55,-45,27,27.5],    'Pacific': [-70,-60,27.7,27.9],  'Indian': [-70,-55,27.6,27.9]} # 15
        ]

        domain_char = {'nb_basins': 3, 'Atlantic': True, 'Pacific': True, 'Indian': True}

    if domain_name == 'North Atlantic':

        # North Atlantic (warmer/saltier)
        domains = [
            {'name':'bcc-csm1-1'    , 'Atlantic': [45,50,26.3,27.2],    'Pacific': None, 'Indian': None}, # 0
            {'name':'CanESM2'       , 'Atlantic': [30,45,26,26.8],    'Pacific': None, 'Indian': None}, # 1
            {'name':'CCSM4'         , 'Atlantic': [25,40,26,27],    'Pacific': None, 'Indian': None}, # 2
            {'name':'CESM1-CAM5'    , 'Atlantic': [20,40,26.,26.9],    'Pacific': None, 'Indian': None}, # 3
            {'name':'CNRM-CM5'      , 'Atlantic': [40,50,26.2,26.9],    'Pacific': None, 'Indian': None}, # 4
            {'name':'CSIRO-Mk3-6-0' , 'Atlantic': [30,50,26.1,27.2],    'Pacific': None, 'Indian': None}, #5
            {'name':'FGOALS-g2'     , 'Atlantic': [25,35,26.1,27.2],    'Pacific': None, 'Indian': None}, # 6
            {'name':'GFDL-CM3'      , 'Atlantic': [45,60,26.4,27.1],    'Pacific': None, 'Indian': None}, # 7
            {'name':'GFDL-ESM2M'    , 'Atlantic': [20,35,26,27],    'Pacific': None, 'Indian': None}, # 8
            {'name':'GISS-E2-R'     , 'Atlantic': [45,65,26.5,27.8],    'Pacific': None, 'Indian': None}, # 9
            {'name':'HadGEM2-ES'    , 'Atlantic': [20,40,26,26.8],    'Pacific': None, 'Indian': None}, # 10
            {'name':'IPSL-CM5A-LR'  , 'Atlantic': [30,45,26.4,27],    'Pacific': None, 'Indian': None}, # 11
            {'name':'IPSL-CM5A-MR'  , 'Atlantic': [30,45,26.5,27.1],    'Pacific': None, 'Indian': None}, # 12
            {'name':'MIROC-ESM-CHEM', 'Atlantic': [35,45,26.4,27.2],    'Pacific': None, 'Indian': None}, # 13
            {'name':'MIROC-ESM'     , 'Atlantic': [45,60,26.6,27.6],    'Pacific': None, 'Indian': None}, # 14
            {'name':'MME'           , 'Atlantic': [35,70,26.3,27],    'Pacific': None, 'Indian': None}, # 15
        ]

        domain_char = {'nb_basins': 1, 'Atlantic': True, 'Pacific': False, 'Indian': False}

    if domain_name == 'Northern ST':

        # Northern Subtropics (cooling/freshening in the Pacific and Indian oceans)
        domains = [
            {'name':'bcc-csm1-1'    , 'Atlantic': None,    'Pacific': [30,35,25,25.2],    'Indian': [22,24,25.5,26.1]}, # 0
            {'name':'CanESM2'       , 'Atlantic': None,    'Pacific': [20,35,24.8,25.2],  'Indian': [20,25,25.4,26.5]}, # 1
            {'name':'CCSM4'         , 'Atlantic': None,    'Pacific': [15,30,24.2,25.2],  'Indian': [20,25,25.8,26.3]}, # 2
            {'name':'CESM1-CAM5'    , 'Atlantic': None,    'Pacific': None,               'Indian': [20,25,25.8,26.5]}, # 3
            {'name':'CNRM-CM5'      , 'Atlantic': None,    'Pacific': [25,30,24.8,25.2],  'Indian': [20,25,25.8,26.5]}, # 4
            {'name':'CSIRO-Mk3-6-0' , 'Atlantic': None,    'Pacific': [10,20,23.4,25],    'Indian': None}, #5
            {'name':'FGOALS-g2'     , 'Atlantic': None,    'Pacific': [15,17,24,25],      'Indian': [23,25,26.6,27]}, # 6
            {'name':'GFDL-CM3'      , 'Atlantic': None,    'Pacific': None,               'Indian': None}, # 7
            {'name':'GFDL-ESM2M'    , 'Atlantic': None,    'Pacific': [15,28,25.4,26],    'Indian': [20,25,25.8,27]}, # 8
            {'name':'GISS-E2-R'     , 'Atlantic': None,    'Pacific': [10,30,24,25.6],    'Indian': [15,25,25.8,26.7]}, # 9
            {'name':'HadGEM2-ES'    , 'Atlantic': None,    'Pacific': [22,32,25.8,26.5],  'Indian': [5,10,23.4,24.4]}, # 10
            {'name':'IPSL-CM5A-LR'  , 'Atlantic': None,    'Pacific': [20,30,25,25.8],    'Indian': [18,22,26.1,26.7]}, # 11
            {'name':'IPSL-CM5A-MR'  , 'Atlantic': None,    'Pacific': [20,35,25,26],      'Indian': [15,20,25.8,26.7]}, # 12
            {'name':'MIROC-ESM-CHEM', 'Atlantic': None,    'Pacific': [25,35,25.2,25.8],  'Indian': [15,25,26,26.7]}, # 13
            {'name':'MIROC-ESM'     , 'Atlantic': None,    'Pacific': [15,25,26,26.3],    'Indian': [15,20,22.8,24]}, # 14
            {'name':'MME'           , 'Atlantic': None,    'Pacific': [20,30,24.8,25.8],  'Indian': [20,25,25.8,26.9]}, # 15
        ]

        domain_char = {'nb_basins': 2, 'Atlantic': False, 'Pacific': True, 'Indian': True}

    if domain_name == 'North Pacific':

        # North Pacific (warmer/saltier)
        domains = [
            {'name':'bcc-csm1-1'    , 'Atlantic': None,    'Pacific': [40,55,26,26.4],   'Indian': None}, # 0
            {'name':'CanESM2'       , 'Atlantic': None,    'Pacific': [45,60,26.5,26.8], 'Indian': None}, # 1
            {'name':'CCSM4'         , 'Atlantic': None,    'Pacific': [40,60,26,26.7],   'Indian': None}, # 2
            {'name':'CESM1-CAM5'    , 'Atlantic': None,    'Pacific': [60,65,26,26.7],   'Indian': None}, # 3
            {'name':'CNRM-CM5'      , 'Atlantic': None,    'Pacific': [50,65,25.8,26.4], 'Indian': None}, # 4
            {'name':'CSIRO-Mk3-6-0' , 'Atlantic': None,    'Pacific': [60,65,25.4,26.3], 'Indian': None}, #5
            {'name':'FGOALS-g2'     , 'Atlantic': None,    'Pacific': [40,60,26.6,26.9], 'Indian': None}, # 6
            {'name':'GFDL-CM3'      , 'Atlantic': None,    'Pacific': [45,60,26.5,26.9], 'Indian': None}, # 7
            {'name':'GFDL-ESM2M'    , 'Atlantic': None,    'Pacific': [38,50,26,26.5],   'Indian': None}, # 8
            {'name':'GISS-E2-R'     , 'Atlantic': None,    'Pacific': [40,60,26.7,27.1], 'Indian': None}, # 9
            {'name':'HadGEM2-ES'    , 'Atlantic': None,    'Pacific': [45,65,26.4,26.9], 'Indian': None}, # 10
            {'name':'IPSL-CM5A-LR'  , 'Atlantic': None,    'Pacific': [45,60,26.8,27.1], 'Indian': None}, # 11
            {'name':'IPSL-CM5A-MR'  , 'Atlantic': None,    'Pacific': [45,60,26.8,27],   'Indian': None}, # 12
            {'name':'MIROC-ESM-CHEM', 'Atlantic': None,    'Pacific': [45,60,27,27.4],   'Indian': None}, # 13
            {'name':'MIROC-ESM'     , 'Atlantic': None,    'Pacific': [50,62,27,27.3],   'Indian': None}, # 14
            {'name':'MME'           , 'Atlantic': None,    'Pacific': [55,65,25.8,26.7],   'Indian': None}, # 15
        ]

        domain_char = {'nb_basins': 1, 'Atlantic': False, 'Pacific': True, 'Indian': False}


    for imodel in range(len(domains)):
        if domains[imodel]['name'] == model_name :
            varout = domains[imodel]

    return varout, domain_char



def ToEdomainrcp85vshistNat(model_name, domain_name):

    '''
    Save domain boxes for rcp8.5 vs. histNat ensemble means (salinity/temperature)
    November-December 2018

    :param model_name: name of model
    :param domain_name: Southern ST, Southern Ocean, etc...
    :return: box(es) of the specified model and domain

    EDIT : March 30th 2020 - some corrections to the coordinates especially in the subtropical North Atlantic 
    to stay within about 20-40ºN (some were betwwen 40-60)
    '''
    if domain_name == 'Southern ST':

        # Southern Subtropics (cooling/freshening in all three basins)
        domains = [
            {'name':'CanESM2'       , 'Atlantic': [-40,-25,25,26.4],    'Pacific': [-40,-10,24.8,26.3], 'Indian': [-40,-20,25.4,26.5]},
            {'name':'CCSM4'         , 'Atlantic': [-40,-15,24.8,26.5],  'Pacific': [-40,-15,24.8,26.4], 'Indian': [-40,-20,25.4,26.3]},
            {'name':'CESM1-CAM5'    , 'Atlantic': [-30,-10,24.8,26.3],  'Pacific': [-40,-15,25.2,26.3], 'Indian': [-35,-20,25.4,26.2]},
            {'name':'CSIRO-Mk3-6-0' , 'Atlantic': [-20,-13,25.6,26.3],  'Pacific': [-35,-20,25.6,26.3], 'Indian': [-40,-20,26,26.6]},
            {'name':'FGOALS-g2'     , 'Atlantic': [-20,-10,25.8,26.6],  'Pacific': [-28,-15,24.8,25.8], 'Indian': [-35,-15,24.6,26.4]},
            {'name':'GISS-E2-R'     , 'Atlantic': [-35,-27,25.6,26.5],  'Pacific': [-30,-15,24,25.8],     'Indian': [-35,-15,25.2,26.1]},
            {'name':'HadGEM2-ES'    , 'Atlantic': [-38,-25,25.8,26.3],  'Pacific': [-30,-15,24.2,25.8], 'Indian': [-35,-15,25.,26]},
            {'name':'IPSL-CM5A-LR'  , 'Atlantic': [-40,-25,26.8,27.4],  'Pacific': [-30,-15,25.8,26.6], 'Indian': [-35,-20,26.3,27]},
            {'name':'IPSL-CM5A-MR'  , 'Atlantic': [-30,-20,25.8,26.5],  'Pacific': [-30,-15,25.4,26.6], 'Indian': [-40,-20,26,27]},
            {'name':'MIROC-ESM-CHEM', 'Atlantic': [-30,-10,25.2,26.1],  'Pacific': [-30,-20,25.6,26.1], 'Indian': [-37,-20,25.4,26.3]},
            {'name':'MIROC-ESM'     , 'Atlantic': [-30,-10,25.2,26.1],  'Pacific': [-30,-15,25.6,26.2], 'Indian': [-37,-20,25.2,26.3]}
        ]
        domain_char = {'nb_basins': 3, 'Atlantic': True, 'Pacific': True, 'Indian': True}

    if domain_name == 'SO':

        # Southern Ocean (warmer/saltier in all three basins)
        domains = [
            {'name':'CanESM2'       , 'Atlantic': [-50,-40,26.8,27.5],  'Pacific': [-60,-50,26.8,27.5],  'Indian': [-60,-45,27,27.7]},
            {'name':'CCSM4'         , 'Atlantic': [-60,-45,27.1,27.9],  'Pacific': [-60,-50,27,27.6],    'Indian': [-60,-45,26.9,27.6]},
            {'name':'CESM1-CAM5'    , 'Atlantic': [-60,-45,26.9,27.8],  'Pacific': [-65,-50,26.8,27.5],  'Indian': [-65,-50,26.8,27.6]},
            {'name':'CSIRO-Mk3-6-0' , 'Atlantic': [-50,-35,27.1,27.8],  'Pacific': [-65,-50,27,27.9],    'Indian': [-60,-50,27.2,28]},
            {'name':'FGOALS-g2'     , 'Atlantic': [-65,-45,26.75,27.5],  'Pacific': [-67,-60,27.2,28],    'Indian': [-60,-45,26.7,27.7]},
            {'name':'GISS-E2-R'     , 'Atlantic': [-55,-45,27,27.5],    'Pacific': [-60,-50,26.5,27.5],  'Indian': [-55,-45,27,27.3]},
            {'name':'HadGEM2-ES'    , 'Atlantic': [-55,-45,26.7,27.3],  'Pacific': [-63,-55,26.5,27.5],  'Indian': [-55,-45,26,27.3]},
            {'name':'IPSL-CM5A-LR'  , 'Atlantic': [-60,-50,27.5,28],    'Pacific': [-65,-60,27.7,27.9],  'Indian': [-55,-50,27.5,27.8]},
            {'name':'IPSL-CM5A-MR'  , 'Atlantic': [-55,-50,27.4,27.9],  'Pacific': [-65,-60,27.5,27.8],  'Indian': [-55,-50,27.6,27.8]},
            {'name':'MIROC-ESM-CHEM', 'Atlantic': [-55,-45,27.2,27.8],  'Pacific': [-65,-50,27,27.8],    'Indian': [-60,-45,27.3,27.9]},
            {'name':'MIROC-ESM'     , 'Atlantic': [-55,-45,27,27.8],    'Pacific': [-65,-50,27.2,27.8],  'Indian': [-62,-50,27.3,27.8]}
        ]
        domain_char = {'nb_basins': 3, 'Atlantic': True, 'Pacific': True, 'Indian': True}

    if domain_name == 'North Atlantic':

        # North Atlantic (warmer/saltier)
        domains = [
            {'name':'CanESM2'       , 'Atlantic': [25,45,26,27],    'Pacific': None, 'Indian': None},
            {'name':'CCSM4'         , 'Atlantic': [25,40,26,27],    'Pacific': None, 'Indian': None},
            {'name':'CESM1-CAM5'    , 'Atlantic': [20,40,26,26.9],  'Pacific': None, 'Indian': None},
            {'name':'CSIRO-Mk3-6-0' , 'Atlantic': [25,45,26.1,27.2],'Pacific': None, 'Indian': None},
            {'name':'FGOALS-g2'     , 'Atlantic': [25,40,26.1,27.3],'Pacific': None, 'Indian': None},
            {'name':'GISS-E2-R'     , 'Atlantic': [25,45,26,27],    'Pacific': None, 'Indian': None},
            {'name':'HadGEM2-ES'    , 'Atlantic': [25,45,25.2,26.6],'Pacific': None, 'Indian': None},
            {'name':'IPSL-CM5A-LR'  , 'Atlantic': [25,45,25.8,27.2],'Pacific': None, 'Indian': None},
            {'name':'IPSL-CM5A-MR'  , 'Atlantic': [25,45,26.1,27.3],'Pacific': None, 'Indian': None},
            {'name':'MIROC-ESM-CHEM', 'Atlantic': [25,45,25.6,26.7],'Pacific': None, 'Indian': None},
            {'name':'MIROC-ESM'     , 'Atlantic': [25,45,25.6,26.7],'Pacific': None, 'Indian': None}
        ]
        domain_char = {'nb_basins': 1, 'Atlantic': True, 'Pacific': False, 'Indian': False}

    if domain_name == 'Northern ST':

        # Northern Subtropics (cooling/freshening in the Pacific and Indian oceans)
        domains = [
            {'name':'CanESM2'       , 'Atlantic': None,    'Pacific': [15,35,23.6,25.2],  'Indian': [15,25,25,26.5]},
            {'name':'CCSM4'         , 'Atlantic': None,    'Pacific': [15,35,24,25.2],  'Indian': [20,25,25.8,26.3]},
            {'name':'CESM1-CAM5'    , 'Atlantic': None,    'Pacific': [15,30,23.8,25.2],  'Indian': [20,25,25.6,26.5]},
            {'name':'CSIRO-Mk3-6-0' , 'Atlantic': None,    'Pacific': [10,30,23.4,25],    'Indian': [20,25,26,26.8]},
            {'name':'FGOALS-g2'     , 'Atlantic': None,    'Pacific': [10,25,23.6,25],    'Indian': [20,25,26,27]},
            {'name':'GISS-E2-R'     , 'Atlantic': None,    'Pacific': [15,30,23.8,25.6],    'Indian': [15,25,25.8,26.7]},
            {'name':'HadGEM2-ES'    , 'Atlantic': None,    'Pacific': [12,30,23.6,25],    'Indian': [20,25,26,26.6]},
            {'name':'IPSL-CM5A-LR'  , 'Atlantic': None,    'Pacific': [15,35,24.5,26],    'Indian': [15,22,26.1,26.7]},
            {'name':'IPSL-CM5A-MR'  , 'Atlantic': None,    'Pacific': [18,35,24,25.8],    'Indian': [15,22,25.8,26.7]},
            {'name':'MIROC-ESM-CHEM', 'Atlantic': None,    'Pacific': [15,35,24.4,26.1],  'Indian': [15,20,26,26.6]},
            {'name':'MIROC-ESM'     , 'Atlantic': None,    'Pacific': [15,35,24.5,26],    'Indian': [15,25,25,26.5]}
        ]

        domain_char = {'nb_basins': 2, 'Atlantic': False, 'Pacific': True, 'Indian': True}

    if domain_name == 'North Pacific':

        # North Pacific (warmer/saltier)
        domains = [
            {'name':'CanESM2'       , 'Atlantic': None,    'Pacific': [45,60,26,26.8],     'Indian': None},
            {'name':'CCSM4'         , 'Atlantic': None,    'Pacific': [45,60,25.8,26.7],   'Indian': None},
            {'name':'CESM1-CAM5'    , 'Atlantic': None,    'Pacific': [45,65,25.8,26.6],   'Indian': None},
            {'name':'CSIRO-Mk3-6-0' , 'Atlantic': None,    'Pacific': [45,60,25.8,26.6],   'Indian': None},
            {'name':'FGOALS-g2'     , 'Atlantic': None,    'Pacific': [40,60,26.,27.3],    'Indian': None},
            {'name':'GISS-E2-R'     , 'Atlantic': None,    'Pacific': [40,60,26.3,27.1],   'Indian': None},
            {'name':'HadGEM2-ES'    , 'Atlantic': None,    'Pacific': [45,60,26,26.9],     'Indian': None},
            {'name':'IPSL-CM5A-LR'  , 'Atlantic': None,    'Pacific': [45,60,26.5,27.1],   'Indian': None},
            {'name':'IPSL-CM5A-MR'  , 'Atlantic': None,    'Pacific': [45,60,26.4,27],     'Indian': None},
            {'name':'MIROC-ESM-CHEM', 'Atlantic': None,    'Pacific': [45,60,26.3,27.25],  'Indian': None},
            {'name':'MIROC-ESM'     , 'Atlantic': None,    'Pacific': [45,60,26.3,27.25],  'Indian': None}
        ]

        domain_char = {'nb_basins': 1, 'Atlantic': False, 'Pacific': True, 'Indian': False}


    for imodel in range(len(domains)):
        if domains[imodel]['name'] == model_name :
            varout = domains[imodel]

    return varout, domain_char


def ToEdomain1pctCO2vsPiC(model_name, domain_name):

    '''
    Save domain boxes for 1pctCO2 vs. pre-industrial control runs of each model (salinity/temperature)

    :param model_name: name of model
    :param domain_name: Southern ST, Southern Ocean, etc...
    :return: box(es) of the specified model and domain

    '''
    if domain_name == 'Southern ST':

        # Southern Subtropics (cooling/freshening in all three basins)
        domains = [
            {'name':'ACCESS1-0'   , 'Atlantic': [-35,-10,26,26.6],    'Pacific': [-40,-25,25.4,26.4], 'Indian': [-40,-25,25.8,26.5]}, # 0
            {'name':'ACCESS1-3'   , 'Atlantic': [-20,-10,25.8,26.5],  'Pacific': [-40,-20,25.2,26.3], 'Indian': [-40,-20,25.4,26.2]}, # 1
            {'name':'BNU-ESM'     , 'Atlantic': [-30,-25,25.4,26],    'Pacific': [-35,-20,24.4,25.8], 'Indian': [-35,-20,24.8,25.6]}, # 2
            {'name':'CCSM4'       , 'Atlantic': [-40,-15,24.8,26.5],  'Pacific': [-40,-15,24.8,26.4], 'Indian': [-40,-20,25.4,26.3]}, # 3
            {'name':'CESM1-BGC'   , 'Atlantic': [-45,-20,24.8,26.5],  'Pacific': [-40,-20,25.2,26.3], 'Indian': [-38,-20,25.6,26.3]}, # 4
            {'name':'CESM1-CAM5'  , 'Atlantic': [-30,-10,24.8,26.3],  'Pacific': [-40,-15,25.2,26.3], 'Indian': [-35,-20,25.4,26.2]}, #5
            {'name':'CNRM-CM5'    , 'Atlantic': None,                 'Pacific': None,                'Indian': None}, # 6
            {'name':'CNRM-CM5-2'  , 'Atlantic': [-35,-15,26.2,26.6],  'Pacific': [-35,-15,25,26.3],   'Indian': [-40,-20,25.8,26.4]}, # 7
            {'name':'CSIRO-Mk3-6-0','Atlantic': [-20,-13,25.6,26.3],  'Pacific': [-35,-20,25.6,26.3], 'Indian': [-40,-20,26,26.6]}, # 8
            {'name':'FGOALS-g2'   , 'Atlantic': [-15,-10,26,27],      'Pacific': [-30,-15,25,25.8],   'Indian': [-35,-20,25.4,26.3]}, # 9
            {'name':'GFDL-ESM2G'  , 'Atlantic': [-53,-48,26.6,27],    'Pacific': [-35,-20,25.4,26.3], 'Indian': [-25,-15,25.2,25.8]}, # 10
            {'name':'GFDL-ESM2M'  , 'Atlantic': [-32,-26,25.4,26.4],  'Pacific': [-35,-20,25.8,26.3], 'Indian': [-40,-20,25.6,26.2]}, # 11
            {'name':'HadGEM2-ES'  , 'Atlantic': [-38,-25,25.8,26.3],  'Pacific': [-30,-15,24.2,25.8], 'Indian': [-30,-20,25.,26]}, # 12
            {'name':'IPSL-CM5A-LR', 'Atlantic': [-40,-25,26.8,27.4],  'Pacific': [-30,-15,25.8,26.6], 'Indian': [-35,-20,26.3,27]}, # 13
            {'name':'IPSL-CM5A-MR', 'Atlantic': [-40,-25,26.4,27.2],  'Pacific': [-30,-15,25.4,26.6], 'Indian': [-40,-20,26,27]}, # 14
            {'name':'IPSL-CM5B-LR', 'Atlantic': [-40,-10,26,27],      'Pacific': [-40,-15,26,26.8],   'Indian': [-40,-20,26,27]} # 15
        ]
        domain_char = {'nb_basins': 3, 'Atlantic': True, 'Pacific': True, 'Indian': True}

    if domain_name == 'SO':

        # Southern Ocean (warmer/saltier in all three basins)
        domains = [
            {'name':'ACCESS1-0'   , 'Atlantic': [-55,-45,26.8,27.3],    'Pacific': [-60,-55,26.8,27.4], 'Indian': [-55,-45,26.5,27.3]}, # 0
            {'name':'ACCESS1-3'   , 'Atlantic': [-55,-45,27,27.5],      'Pacific': [-60,-50,26.7,27.5], 'Indian': [-55,-45,26.6,27.4]}, # 1
            {'name':'BNU-ESM'     , 'Atlantic': [-65,-50,26.8,27.6],    'Pacific': [-65,-50,26.3,27.3], 'Indian': [-65,-45,26.8,27.7]}, # 2
            {'name':'CCSM4'       , 'Atlantic': [-60,-45,27.1,27.9],    'Pacific': [-60,-50,27,27.6],   'Indian': [-60,-45,26.9,27.6]},# 3
            {'name':'CESM1-BGC'   , 'Atlantic': [-55,-45,27.1,27.8],    'Pacific': [-60,-53,27,27.6],   'Indian': [-55,-45,26.9,27.6]}, # 4
            {'name':'CESM1-CAM5'  , 'Atlantic': [-60,-45,26.9,27.8],    'Pacific': [-60,-50,26.8,27.5], 'Indian': [-65,-50,26.8,27.6]}, #5
            {'name':'CNRM-CM5'    , 'Atlantic': None,                   'Pacific': None,                'Indian': None}, # 6
            {'name':'CNRM-CM5-2'  , 'Atlantic': [-57,-45,27,27.4],      'Pacific': [-60,-45,26.5,27.2], 'Indian': [-55,-45,26.5,27.2]}, # 7
            {'name':'CSIRO-Mk3-6-0','Atlantic': [-50,-35,27.1,27.8],    'Pacific': [-65,-50,27,27.9],   'Indian': [-60,-50,27.2,28]}, # 8
            {'name':'FGOALS-g2'   , 'Atlantic': [-40,-30,26.5,26.9],    'Pacific': [-65,-55,27.1,28],   'Indian': [-53,-45,27.1,27.6]}, # 9
            {'name':'GFDL-ESM2G'  , 'Atlantic': [-55,-50,27.3,27.6],    'Pacific': None,                'Indian': [-55,-45,27,27.5]}, # 10
            {'name':'GFDL-ESM2M' ,  'Atlantic': [-55,-35,27.1,27.5],    'Pacific': [-70,-55,27.5,27.9], 'Indian': [-53,-45,27.2,27.5]}, # 11
            {'name':'HadGEM2-ES'  , 'Atlantic': [-55,-45,26.7,27.3],    'Pacific': [-63,-55,26.5,27.5], 'Indian': [-55,-45,26,27.3]}, # 12
            {'name':'IPSL-CM5A-LR', 'Atlantic': [-60,-50,27.5,28],      'Pacific': [-65,-60,27.7,27.9], 'Indian': [-55,-50,27.5,27.8]}, # 13
            {'name':'IPSL-CM5A-MR', 'Atlantic': [-55,-50,27.4,27.9],    'Pacific': [-65,-60,27.5,27.8], 'Indian': [-55,-50,27.6,27.8]}, # 14
            {'name':'IPSL-CM5B-LR', 'Atlantic': [-65,-55,27.5,28],      'Pacific': [-70,-60,27.8,28],   'Indian': [-70,-56,27.7,28]} # 15
        ]

        domain_char = {'nb_basins': 3, 'Atlantic': True, 'Pacific': True, 'Indian': True}

    if domain_name == 'North Atlantic':

        # North Atlantic (warmer/saltier)
        domains = [
            {'name':'ACCESS1-0'   , 'Atlantic': [20,35,24.8,26.2], 'Pacific': None, 'Indian': None}, # 0
            {'name':'ACCESS1-3'   , 'Atlantic': [20,40,25,26.5], 'Pacific': None, 'Indian': None}, # 1
            {'name':'BNU-ESM'     , 'Atlantic': [25,45,25.8,27.2], 'Pacific': None, 'Indian': None}, # 2
            {'name':'CCSM4'       , 'Atlantic': [25,40,26,27], 'Pacific': None, 'Indian': None}, # 3
            {'name':'CESM1-BGC'   , 'Atlantic': [20,40,26,27.1], 'Pacific': None, 'Indian': None}, # 4
            {'name':'CESM1-CAM5'  , 'Atlantic': [20,40,26,26.9], 'Pacific': None, 'Indian': None}, #5
            {'name':'CNRM-CM5'    , 'Atlantic': None,          'Pacific': None, 'Indian': None}, # 6
            {'name':'CNRM-CM5-2'  , 'Atlantic': [25,45,26,27], 'Pacific': None, 'Indian': None}, # 7
            {'name':'CSIRO-Mk3-6-0','Atlantic': [25,45,26.1,27.2], 'Pacific': None, 'Indian': None}, # 8
            {'name':'FGOALS-g2'   , 'Atlantic': [25,35,26.1,27.2], 'Pacific': None, 'Indian': None}, # 9
            {'name':'GFDL-ESM2G'  , 'Atlantic': [20,30,25.6,27], 'Pacific': None, 'Indian': None}, # 10
            {'name':'GFDL-ESM2M'  , 'Atlantic': [20,40,26,27], 'Pacific': None, 'Indian': None}, # 11
            {'name':'HadGEM2-ES'  , 'Atlantic': [25,45,25.2,26.6], 'Pacific': None, 'Indian': None}, # 12
            {'name':'IPSL-CM5A-LR', 'Atlantic': [25,45,25.8,27.2],'Pacific': None, 'Indian': None}, # 13
            {'name':'IPSL-CM5A-MR', 'Atlantic': [30,45,26.3,27.3],'Pacific': None, 'Indian': None}, # 14
            {'name':'IPSL-CM5B-LR', 'Atlantic': [35,45,25.5,26.5], 'Pacific': None, 'Indian': None} # 15
        ]

        domain_char = {'nb_basins': 1, 'Atlantic': True, 'Pacific': False, 'Indian': False}

    if domain_name == 'Northern ST':

        # Northern Subtropics (cooling/freshening in the Pacific and Indian oceans)
        domains = [
            {'name':'ACCESS1-0'   , 'Pacific': [20,35,24.2,25.8],  'Indian': [10,20,25.6,26.6], 'Atlantic': None}, # 0
            {'name':'ACCESS1-3'   , 'Pacific': [15,35,24,25.6],    'Indian': [20,25,25.8,27],   'Atlantic': None}, # 1
            {'name':'BNU-ESM'     , 'Pacific': [15,35,24,25.4],    'Indian': [15,25,26,26.8],   'Atlantic': None}, # 2
            {'name':'CCSM4'       , 'Pacific': [15,30,23.8,25.2],  'Indian': [20,25,25.6,26.5], 'Atlantic': None}, # 3
            {'name':'CESM1-BGC'   , 'Pacific': [20,35,24,25.2],    'Indian': [20,25,25.6,26.4], 'Atlantic': None}, # 4
            {'name':'CESM1-CAM5'  , 'Pacific': [15,30,24,25.4],    'Indian': [20,25,25.6,26.5], 'Atlantic': None}, #5
            {'name':'CNRM-CM5'    , 'Pacific': None,               'Indian': None,              'Atlantic': None}, # 6
            {'name':'CNRM-CM5-2'  , 'Pacific': [15,30,24.4,25.8],  'Indian': [15,25,26,27],     'Atlantic': None}, # 7
            {'name':'CSIRO-Mk3-6-0','Pacific': [10,30,23.4,25],    'Indian': [20,25,26,26.8],   'Atlantic': None}, # 8
            {'name':'FGOALS-g2'   , 'Pacific': [13,20,24,24.8],    'Indian': [22,25,26.6,27.1], 'Atlantic': None}, # 9
            {'name':'GFDL-ESM2G'  , 'Pacific': [17,30,25.2,25.8],  'Indian': [20,25,25.4,26.9], 'Atlantic': None}, # 10
            {'name':'GFDL-ESM2M'  , 'Pacific': [15,30,25.2,26.1],  'Indian': [20,25,25.8,27],   'Atlantic': None}, # 11
            {'name':'HadGEM2-ES'  , 'Pacific': [20,35,23.8,25],    'Indian': [20,25,26,26.6],   'Atlantic': None}, # 12
            {'name':'IPSL-CM5A-LR', 'Pacific': [20,30,24.8,26],    'Indian': [15,22,26.1,26.7], 'Atlantic': None}, # 13
            {'name':'IPSL-CM5A-MR', 'Pacific': [20,35,24.4,26],    'Indian': [15,22,25.8,26.7], 'Atlantic': None}, # 14
            {'name':'IPSL-CM5B-LR', 'Pacific': [15,35,25,26.1],    'Indian': [10,20,26.3,27],   'Atlantic': None} # 15
        ]

        domain_char = {'nb_basins': 2, 'Atlantic': False, 'Pacific': True, 'Indian': True}

    if domain_name == 'North Pacific':

        # North Pacific (warmer/saltier)
        domains = [
            {'name':'ACCESS1-0'   , 'Pacific': [43,60,26.2,27],  'Atlantic': None, 'Indian': None}, # 0
            {'name':'ACCESS1-3'   , 'Pacific': [45,60,26.1,27],  'Atlantic': None, 'Indian': None}, # 1
            {'name':'BNU-ESM'     , 'Pacific': [45,60,26,27.2],  'Atlantic': None, 'Indian': None}, # 2
            {'name':'CCSM4'       , 'Pacific': [45,65,25.8,26.6],'Atlantic': None, 'Indian': None}, # 3
            {'name':'CESM1-BGC'   , 'Pacific': [45,60,26,27],    'Atlantic': None, 'Indian': None}, # 4
            {'name':'CESM1-CAM5'  , 'Pacific': [45,65,25.8,26.6], 'Atlantic': None, 'Indian': None}, #5
            {'name':'CNRM-CM5'    , 'Pacific': None,              'Atlantic': None, 'Indian': None}, # 6
            {'name':'CNRM-CM5-2'  , 'Pacific': [45,60,26,27],     'Atlantic': None, 'Indian': None}, # 7
            {'name':'CSIRO-Mk3-6-0','Pacific': [45,60,25.8,26.6], 'Atlantic': None, 'Indian': None}, # 8
            {'name':'FGOALS-g2'   , 'Pacific': [40,60,26.6,26.9], 'Atlantic': None, 'Indian': None}, # 9
            {'name':'GFDL-ESM2G'  , 'Pacific': [45,60,26.4,26.7], 'Atlantic': None, 'Indian': None}, # 10
            {'name':'GFDL-ESM2M'  , 'Pacific': [45,60,26,27    ], 'Atlantic': None, 'Indian': None}, # 11
            {'name':'HadGEM2-ES'  , 'Pacific': [45,60,26,26.9],   'Atlantic': None, 'Indian': None}, # 12
            {'name':'IPSL-CM5A-LR' ,'Pacific': [45,60,26.5,27.1], 'Atlantic': None, 'Indian': None}, # 13
            {'name':'IPSL-CM5A-MR', 'Pacific': [45,60,26.4,27],   'Atlantic': None, 'Indian': None}, # 14
            {'name':'IPSL-CM5B-LR', 'Pacific': [40,60,26.4,26.7], 'Atlantic': None, 'Indian': None} # 15
        ]

        domain_char = {'nb_basins': 1, 'Atlantic': False, 'Pacific': True, 'Indian': False}


    for imodel in range(len(domains)):
        if domains[imodel]['name'] == model_name :
            varout = domains[imodel]

    return varout, domain_char






def ToEdomain1pctCO2(model_name, domain_name):

    '''
    Save domain boxes for 1pctCO2 runs of each model (salinity/temperature)

    :param model_name: name of model
    :param domain_name: Southern ST, Southern Ocean, etc...
    :return: box(es) of the specified model and domain

    '''
    if domain_name == 'Southern ST':

        # Southern Subtropics (cooling/freshening in all three basins)
        domains = [
            {'name':'ACCESS1-0'   , 'Atlantic': [-20,-10,26,26.7],    'Pacific': [-37,-25,25.4,26.4], 'Indian': [-37,-25,25.8,26.5]}, # 0
            {'name':'ACCESS1-3'   , 'Atlantic': [-20,-10,25.8,26.5],  'Pacific': [-35,-20,25.2,26.3], 'Indian': [-40,-20,25.4,26.2]}, # 1
            {'name':'BNU-ESM'     , 'Atlantic': [-15,-7,24.8,25.8],   'Pacific': [-35,-20,25,25.8],   'Indian': [-30,-20,25,25.4]}, # 2
            {'name':'CCSM4'       , 'Atlantic': [-35,-15,25,26],      'Pacific': [-40,-15,24.8,26.4], 'Indian': [-37,-20,25.6,26.2]}, # 3
            {'name':'CESM1-BGC'   , 'Atlantic': [-40,-30,25.6,26.5],  'Pacific': [-40,-20,25.2,26.4], 'Indian': [-38,-20,25.6,26.3]}, # 4
            {'name':'CESM1-CAM5'  , 'Atlantic': [-20,-7,25,26.4],     'Pacific': [-40,-20,25.2,26.25],'Indian': [-35,-20,25.6,26.2]}, #5
            {'name':'CNRM-CM5'    , 'Atlantic': [-20,-10,24.4,26.2],  'Pacific': [-25,-10,24.6,25.6], 'Indian': [-35,-15,25.4,26.2]}, # 6
            {'name':'CNRM-CM5-2'  , 'Atlantic': [-35,-20,26.2,26.6],  'Pacific': [-40,-20,25.2,26.3], 'Indian': [-40,-15,25.6,26.5]}, # 7
            {'name':'CSIRO-Mk3-6-0','Atlantic': [-20,-12,26,26.3],    'Pacific': [-35,-15,25.2,26.3], 'Indian': [-40,-20,26,26.6]}, # 8
            {'name':'FGOALS-g2'   , 'Atlantic': [-15,-10,26,26.6],    'Pacific': [-25,-15,24.6,25.6], 'Indian': [-35,-20,25.4,26.3]}, # 9
            {'name':'GFDL-ESM2G'  , 'Atlantic': [-53,-48,26.6,27],    'Pacific': [-35,-20,25.4,26.3], 'Indian': [-25,-15,25.2,25.8]}, # 10
            {'name':'GFDL-ESM2M'  , 'Atlantic': [-27,-20,25.6,26],    'Pacific': [-33,-20,25.8,26.3], 'Indian': [-35,-20,25.6,26.2]}, # 11
            #{'name':'GISS-E2-H'   , 'Atlantic': [-30,-20,25.6,26.1],  'Pacific': [-25,-17,25.4,26.2], 'Indian': [-40,-20,25.8,26.4]}, # 12
            {'name':'HadGEM2-ES'  , 'Atlantic': [-38,-33,26,26.5],    'Pacific': [-35,-10,24.8,25.8], 'Indian': [-35,-15,25.4,26]}, # 13
            {'name':'IPSL-CM5A-LR', 'Atlantic': [-43,-35,27.2,27.5],  'Pacific': [-25,-15,25.8,26.6], 'Indian': [-35,-20,26.5,26.9]}, # 14
            {'name':'IPSL-CM5A-MR', 'Atlantic': [-30,-25,26.2,26.6],  'Pacific': [-30,-15,25.4,26.6], 'Indian': [-35,-20,26.3,26.8]}, # 15
            {'name':'IPSL-CM5B-LR', 'Atlantic': [-37,-25,26.4,27.1],  'Pacific': [-25,-10,25.2,26.1], 'Indian': [-35,-20,26.1,26.5]} # 16
        ]
        domain_char = {'nb_basins': 3, 'Atlantic': True, 'Pacific': True, 'Indian': True}

    if domain_name == 'SO':

        # Southern Ocean (warmer/saltier in all three basins)
        domains = [
            {'name':'ACCESS1-0'   , 'Atlantic': [-55,-40,26.7,27.4],    'Pacific': [-62,-55,27.1,27.7], 'Indian': [-55,-45,26.7,27.6]}, # 0
            {'name':'ACCESS1-3'   , 'Atlantic': [-55,-50,27,27.4],      'Pacific': [-60,-50,26.7,27.3], 'Indian': [-55,-45,26.6,27.4]}, # 1
            {'name':'BNU-ESM'     , 'Atlantic': [-65,-50,27.4,27.8],    'Pacific': [-65,-55,26.9,27.5], 'Indian': [-65,-50,27.2,27.8]}, # 2
            {'name':'CCSM4'       , 'Atlantic': [-57,-50,27.5,27.9],    'Pacific': [-57,-50,27,27.6],   'Indian': [-55,-45,27,27.6]}, # 3
            {'name':'CESM1-BGC'   , 'Atlantic': [-55,-50,27.4,27.8],    'Pacific': [-60,-53,27,27.6],   'Indian': [-55,-45,26.9,27.6]}, # 4
            {'name':'CESM1-CAM5'  , 'Atlantic': [-55,-47,26.9,27.7],    'Pacific': [-60,-50,26.8,27.5], 'Indian': [-55,-50,26.9,27.5]}, #5
            {'name':'CNRM-CM5'    , 'Atlantic': [-52,-45,26.5,27.1],    'Pacific': [-55,-50,26.4,27.1], 'Indian': [-52,-47,26.5,27.1]}, # 6
            {'name':'CNRM-CM5-2'  , 'Atlantic': [-50,-40,27.2,27.7],    'Pacific': [-55,-45,26.6,27.2], 'Indian': [-52,-45,26.6,27.3]}, # 7
            {'name':'CSIRO-Mk3-6-0','Atlantic': [-45,-40,27.1,27.6],    'Pacific': [-60,-50,27.3,27.9], 'Indian': [-55,-50,27.3,27.9]}, # 8
            {'name':'FGOALS-g2'   , 'Atlantic': [-40,-30,26.5,26.9],    'Pacific': [-70,-55,27.1,27.8], 'Indian': [-53,-45,27.1,27.6]}, # 9
            {'name':'GFDL-ESM2G'  , 'Atlantic': [-55,-50,27.3,27.6],    'Pacific': [-60,-50,27,27.7],   'Indian': [-55,-45,27,27.5]}, # 10
            {'name':'GFDL-ESM2M' ,  'Atlantic': [-55,-50,27.2,27.6],    'Pacific': [-70,-65,27.7,27.9], 'Indian': [-53,-45,27.2,27.5]}, # 11
            #{'name':'GISS-E2-H'   , 'Atlantic': [-65,-55,27.5,27.9],    'Pacific': [-65,-60,27.7,27.9], 'Indian': [-58,-50,27.6,27.9]}, # 12
            {'name':'HadGEM2-ES'  , 'Atlantic': [-55,-45,26.8,27.2],    'Pacific': [-55,-45,26.3,26.9], 'Indian': [-55,-45,26.6,27.2]}, # 13
            {'name':'IPSL-CM5A-LR', 'Atlantic': [-50,-40,27.7,27.9],    'Pacific': [-60,-55,27.7,27.9], 'Indian': None}, # 14
            {'name':'IPSL-CM5A-MR', 'Atlantic': [-53,-50,27.6,27.9],    'Pacific': [-60,-55,27.6,27.9], 'Indian': [-55,-45,27.8,27.9]}, # 15
            {'name':'IPSL-CM5B-LR', 'Atlantic': [-55,-50,27.3,27.7],    'Pacific': [-60,-55,27.4,27.7], 'Indian': [-60,-55,27.6,27.8]} # 16
        ]

        domain_char = {'nb_basins': 3, 'Atlantic': True, 'Pacific': True, 'Indian': True}

    if domain_name == 'North Atlantic':

        # North Atlantic (warmer/saltier)
        domains = [
            {'name':'ACCESS1-0'   , 'Atlantic': [40,55,25.8,26.9], 'Pacific': None, 'Indian': None}, # 0
            {'name':'ACCESS1-3'   , 'Atlantic': [30,50,26,26.5], 'Pacific': None, 'Indian': None}, # 1
            {'name':'BNU-ESM'     , 'Atlantic': [30,45,26,27.2], 'Pacific': None, 'Indian': None}, # 2
            {'name':'CCSM4'       , 'Atlantic': [25,40,26,27], 'Pacific': None, 'Indian': None}, # 3
            {'name':'CESM1-BGC'   , 'Atlantic': [25,40,26,27.1], 'Pacific': None, 'Indian': None}, # 4
            {'name':'CESM1-CAM5'  , 'Atlantic': [20,40,26,27], 'Pacific': None, 'Indian': None}, #5
            {'name':'CNRM-CM5'    , 'Atlantic': [40,50,26.2,26.9], 'Pacific': None, 'Indian': None}, # 6
            {'name':'CNRM-CM5-2'  , 'Atlantic': [30,45,26.4,27], 'Pacific': None, 'Indian': None}, # 7
            {'name':'CSIRO-Mk3-6-0','Atlantic': [20,40,26,26.8], 'Pacific': None, 'Indian': None}, # 8
            {'name':'FGOALS-g2'   , 'Atlantic': [25,35,26.1,27.2], 'Pacific': None, 'Indian': None}, # 9
            {'name':'GFDL-ESM2G'  , 'Atlantic': [10,30,26,27], 'Pacific': None, 'Indian': None}, # 10
            {'name':'GFDL-ESM2M'  , 'Atlantic': [20,35,26,27], 'Pacific': None, 'Indian': None}, # 11
            #{'name':'GISS-E2-H'   , 'Atlantic': [20,40,26,27], 'Pacific': None, 'Indian': None}, # 12
            {'name':'HadGEM2-ES'  , 'Atlantic': [20,45,26,26.7], 'Pacific': None, 'Indian': None}, # 13
            {'name':'IPSL-CM5A-LR', 'Atlantic': [30,45,26.6,27.4], 'Pacific': None, 'Indian': None}, # 14
            {'name':'IPSL-CM5A-MR', 'Atlantic': [40,45,26.8,27.3], 'Pacific': None, 'Indian': None}, # 15
            {'name':'IPSL-CM5B-LR', 'Atlantic': [40,50,26,26.9], 'Pacific': None, 'Indian': None} # 16
        ]

        domain_char = {'nb_basins': 1, 'Atlantic': True, 'Pacific': False, 'Indian': False}

    if domain_name == 'Northern ST':

        # Northern Subtropics (cooling/freshening in the Pacific and Indian oceans)
        domains = [
            {'name':'ACCESS1-0'   , 'Pacific': [25,32.5,24.6,25.4],'Indian': [5,15,25.6,26.5], 'Atlantic': None}, # 0
            {'name':'ACCESS1-3'   , 'Pacific': [12,30,24.4,25.6],  'Indian': [20,25,26.2,26.7], 'Atlantic': None}, # 1
            {'name':'BNU-ESM'     , 'Pacific': [15,30,24.2,25.4],  'Indian': [15,25,26,26.8], 'Atlantic': None}, # 2
            {'name':'CCSM4'       , 'Pacific': [15,35,24.2,25.6],  'Indian': [20,25,25.8,26.3], 'Atlantic': None}, # 3
            {'name':'CESM1-BGC'   , 'Pacific': [15,35,24.2,25.4],  'Indian': [20,25,25.4,26.2], 'Atlantic': None}, # 4
            {'name':'CESM1-CAM5'  , 'Pacific': [15,30,24.2,25.2],  'Indian': [10,25,25.6,26.3], 'Atlantic': None}, #5
            {'name':'CNRM-CM5'    , 'Pacific': [15,30,24.8,25.6],  'Indian': [10,25,25.8,26.5], 'Atlantic': None}, # 6
            {'name':'CNRM-CM5-2'  , 'Pacific': [15,32,24.6,25.8],  'Indian': [15,25,26,26.6], 'Atlantic': None}, # 7
            {'name':'CSIRO-Mk3-6-0','Pacific': [15,35,24.4,25.4],  'Indian': [15,25,26.2,26.9], 'Atlantic': None}, # 8
            {'name':'FGOALS-g2'   , 'Pacific': [10,25,24,25.2],    'Indian': [5,20,26,26.5], 'Atlantic': None}, # 9
            {'name':'GFDL-ESM2G'  , 'Pacific': [17,30,25.2,25.8],  'Indian': [20,25,25.4,26.4], 'Atlantic': None}, # 10
            {'name':'GFDL-ESM2M'  , 'Pacific': [15,30,25.2,26.1],  'Indian': [20,25,25.8,27], 'Atlantic': None}, # 11
            #{'name':'GISS-E2-H'   , 'Pacific': [15,35,25.2,26.1],  'Indian': [15,25,26,27.2], 'Atlantic': None}, # 12
            {'name':'HadGEM2-ES'  , 'Pacific': [20,35,24.6,25.4],  'Indian': [5,25,25.6,26.3], 'Atlantic': None}, # 13
            {'name':'IPSL-CM5A-LR', 'Pacific': [15,35,24.8,26.2],  'Indian': [15,20,26,26.9], 'Atlantic': None}, # 14
            {'name':'IPSL-CM5A-MR', 'Pacific': [20,35,25,26],      'Indian': [15,20,25.8,26.7], 'Atlantic': None}, # 15
            {'name':'IPSL-CM5B-LR', 'Pacific': [15,30,25,25.8],    'Indian': [15,20,26.5,26.9], 'Atlantic': None} # 16
        ]

        domain_char = {'nb_basins': 2, 'Atlantic': False, 'Pacific': True, 'Indian': True}

    if domain_name == 'North Pacific':

        # North Pacific (warmer/saltier)
        domains = [
            {'name':'ACCESS1-0'   , 'Pacific': [43,60,26.2,27], 'Atlantic': None, 'Indian': None}, # 0
            {'name':'ACCESS1-3'   , 'Pacific': [45,60,26.1,27], 'Atlantic': None, 'Indian': None}, # 1
            {'name':'BNU-ESM'     , 'Pacific': [50,60,26.2,27.3], 'Atlantic': None, 'Indian': None}, # 2
            {'name':'CCSM4'       , 'Pacific': [40,60,26,26.7], 'Atlantic': None, 'Indian': None}, # 3
            {'name':'CESM1-BGC'   , 'Pacific': [40,60,26,26.8], 'Atlantic': None, 'Indian': None}, # 4
            {'name':'CESM1-CAM5'  , 'Pacific': [45,60,26,26.8], 'Atlantic': None, 'Indian': None}, #5
            {'name':'CNRM-CM5'    , 'Pacific': [45,60,26.5,26.9], 'Atlantic': None, 'Indian': None}, # 6
            {'name':'CNRM-CM5-2'  , 'Pacific': [45,60,26.7,27.1], 'Atlantic': None, 'Indian': None}, # 7
            {'name':'CSIRO-Mk3-6-0','Pacific': [40,60,26.1,26.6], 'Atlantic': None, 'Indian': None}, # 8
            {'name':'FGOALS-g2'   , 'Pacific': [40,60,26.6,26.8], 'Atlantic': None, 'Indian': None}, # 9
            {'name':'GFDL-ESM2G'  , 'Pacific': [45,60,26.4,26.7], 'Atlantic': None, 'Indian': None}, # 10
            {'name':'GFDL-ESM2M'  , 'Pacific': [52,60,26.3,26.9], 'Atlantic': None, 'Indian': None}, # 11
            #{'name':'GISS-E2-H'   , 'Pacific': [55,65,27.5,27.8], 'Atlantic': None, 'Indian': None}, # 12
            {'name':'HadGEM2-ES'  , 'Pacific': [45,65,26.5,27], 'Atlantic': None, 'Indian': None}, # 13
            {'name':'IPSL-CM5A-LR', 'Pacific': [40,60,26.8,27.1], 'Atlantic': None, 'Indian': None}, # 14
            {'name':'IPSL-CM5A-MR', 'Pacific': [40,60,26.7,27], 'Atlantic': None, 'Indian': None}, # 15
            {'name':'IPSL-CM5B-LR', 'Pacific': [45,60,26.4,26.7], 'Atlantic': None, 'Indian': None} # 16
        ]

        domain_char = {'nb_basins': 1, 'Atlantic': False, 'Pacific': True, 'Indian': False}


    for imodel in range(len(domains)):
        if domains[imodel]['name'] == model_name :
            varout = domains[imodel]

    return varout, domain_char
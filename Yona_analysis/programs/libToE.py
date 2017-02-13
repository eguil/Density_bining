import numpy as np
import time as timc

def findToE(signal, noise, mult):
    '''
    define Time of Emergence (ToE) from last time index at which signal is larger than mult*noise
        signal is [time,space]
        noise is [space]
        mult is float
    TODO: add valmask where ToE not reached
    '''
    #tcpu0 = timc.clock()
    timN = signal.shape[0]
    toe_wrk = np.ma.ones(signal.shape)*1. # init toe_wrk array to 1
    signaltile = np.reshape(np.tile(noise,timN),signal.shape) # repeat noise timN
    toe_idx = np.argwhere(abs(signal) >= mult*signaltile) # find indices of points where signal > noise
    toe_wrk[toe_idx[:,0],toe_idx[:,1]] = 0. # set corresponding points in toe_wrk to zero
    toe = timN-np.flipud(toe_wrk).argmax(axis=0) # compute ToE as last index when signal > noise
    #tcpu1 = timc.clock()
    # perf
    #print ' ToE CPU = ',tcpu1-tcpu0


    return toe


def ToEdomain(model_name, domain_name):

    if domain_name == 'Southern ST':

        # Southern Subtropics
        domains = [
            {'name':'ACCESS1-0'   , 'Atlantic': [-23,-10,26,26.7],    'Pacific': [-40,-23,25.1,26.4], 'Indian': [-40,-25,25.5,26.5]}, # 0
            {'name':'ACCESS1-3'   , 'Atlantic': [-25,-10,25.75,26.6], 'Pacific': [-40,-20,25,26.3],   'Indian': [-40,-20,25.1,26.25]}, # 1
            {'name':'BNU-ESM'     , 'Atlantic': [-18,-8,25,26],       'Pacific': [-35,-10,23,26],     'Indian': [-30,-20,24.8,25.5]}, # 2
            {'name':'CCSM4'       , 'Atlantic': [-40,-15,24.8,25.8],  'Pacific': [-45,-15,24.75,26.3],'Indian': [-40,-20,25.5,26.25]}, # 3
            {'name':'CESM1-BGC'   , 'Atlantic': [-42,-20,25.25,26],   'Pacific': [-45,-20,25,26.5],   'Indian': [-40,-20,25.5,26.3]}, # 4
            {'name':'CESM1-CAM5'  , 'Atlantic': [-20,-5,24.8,26.5],   'Pacific': [-40,-15,25,26.3],   'Indian': [-40,-20,25.5,26.8]}, #5
            {'name':'CNRM-CM5'    , 'Atlantic': [-35,-5,24.5,26.3],   'Pacific': [-35,-10,24.25,26],  'Indian': [-40,-10,25,26.3]}, # 6
            {'name':'CNRM-CM5-2'  , 'Atlantic': [-35,-10,26.1,26.7],  'Pacific': [-40,-15,24.8,26.5], 'Indian': [-40,-15,25.5,26.6]}, # 7
            {'name':'CSIRO-Mk3-6-0','Atlantic': [-20,-10,25.8,26.5],  'Pacific': [-42,-10,24.8,26.5], 'Indian': [-40,-20,25.8,267]}, # 8
            {'name':'FGOALS-g2'   , 'Atlantic': [-53,-45,26.5,26.7],  'Pacific': [-35,-15,24.3,26],   'Indian': [-40,-18,25,26.5]}, # 9
            {'name':'GFDL-ESM2G'  , 'Atlantic': [-65,-50,26.25,27.25],'Pacific': [-40,-20,25,26.3],   'Indian': [-30,-20,25,26]}, # 10
            {'name':'GFDL-ESM2M'  , 'Atlantic': [-30,-20,25.2,26.2],  'Pacific': [-40,-15,25.5,26.5], 'Indian': [-40,-20,25.5,26.2]}, # 11
            {'name':'GISS-E2-H'   , 'Atlantic': [-35,-20,25.5,26.3],  'Pacific': [-35,-15,25.1,26.5], 'Indian': [-40,-20,25.5,26.7]}, # 12
            {'name':'HadGEM2-ES'  , 'Atlantic': [-40,-30,26,26.7],    'Pacific': [-35,-10,24.5,26],   'Indian': [-38,-15,25.2,26.1]}, # 13
            {'name':'IPSL-CM5A-LR', 'Atlantic': [-45,-25,27.1,27.5],  'Pacific': [-40,-10,25.2,27.2], 'Indian': [-40,-20,26.25,27.3]}, # 14
            {'name':'IPSL-CM5A-MR', 'Atlantic': [-35,-15,26.25,27.25],'Pacific': [-35,-10,24.7,27],   'Indian': [-40,-15,26,27.2]}, # 15
            {'name':'IPSL-CM5B-LR', 'Atlantic': [-35,-20,25.8,27.2],  'Pacific': [-35,-10,24.7,26.6], 'Indian': [-40,-20,25.8,26.8]} # 16
        ]


    if domain_name == 'SO':

        # Southern Ocean
        domains = [
            {'name':'ACCESS1-0'   , 'Atlantic': [-57,-37,26.6,27.6],    'Pacific': [-62,-50,26.75,27.8], 'Indian': [-55,-45,26.25,27.7]}, # 0
            {'name':'ACCESS1-3'   , 'Atlantic': [-57,-45,26.8,27.6],    'Pacific': [-62,-45,26.2,27.8], 'Indian': [-55,-45,26.2,27.8]}, # 1
            {'name':'BNU-ESM'     , 'Atlantic': [-65,-30,26.25,27.1],    'Pacific': [-65,-45,26.2,27.5], 'Indian': [-70,-40,25.5,27.7]}, # 2
            {'name':'CCSM4'       , 'Atlantic': [-60,-40,27,27.9],    'Pacific': [-62,-47,26.75,27.75], 'Indian': [-58,-45,26.3,27.75]}, # 3
            {'name':'CESM1-BGC'   , 'Atlantic': [-60,-45,26.8,27.9],    'Pacific': [-62,-50,26.75,27.8], 'Indian': [-60,-45,26.3,27.9]}, # 4
            {'name':'CESM1-CAM5'  , 'Atlantic': [-62,-45,26.7,27.75],    'Pacific': [-65,-45,26.6,27.75], 'Indian': [-67,-50,26.5,28]}, #5
            {'name':'CNRM-CM5'    , 'Atlantic': [-60,-40,26.4,27.7],    'Pacific': [-65,-45,26,27.6], 'Indian': [-62,-45,26,27.7]}, # 6
            {'name':'CNRM-CM5-2'  , 'Atlantic': [-50,-40,25.7,26.7],    'Pacific': [-62,-45,26.4,27.7], 'Indian': [-55,-45,26.1,27.5]}, # 7
            {'name':'CSIRO-Mk3-6-0','Atlantic': [-50,-35,27,27.8],    'Pacific': [-65,-50,26.8,27.9], 'Indian': [-62,-45,27.2,28]}, # 8
            {'name':'FGOALS-g2'   , 'Atlantic': [-42,-30,26.4,27.25],    'Pacific': [-68,-50,26.75,27.7], 'Indian': [-65,-45,26.7,27.7]}, # 9
            {'name':'GFDL-ESM2G'  , 'Atlantic': [-58,-40,27.3,27.6],    'Pacific': [-58,-48,26.4,27.6], 'Indian': [-55,-40,26.9,27.7]}, # 10
            {'name':'GFDL-ESM2M'  , 'Atlantic': [-55,-50,27.2,28],    'Pacific': [-70,-62,27.2,27.8], 'Indian': [-65,-40,27.1,27.6]}, # 11
            {'name':'GISS-E2-H'   , 'Atlantic': [-65,-53,27.3,27.9],    'Pacific': [-65,-55,27.4,27.8], 'Indian': [-55,-45,26.7,27.75]}, # 12
            {'name':'HadGEM2-ES'  , 'Atlantic': [-60,-40,26.6,27.5],    'Pacific': [-65,-57,26.2,27.7], 'Indian': [-60,-45,25.9,27.5]}, # 13
            {'name':'IPSL-CM5A-LR', 'Atlantic': [-55,-45,27.1,27.8],    'Pacific': [-65,-55,27.5,28], 'Indian': None}, # 14
            {'name':'IPSL-CM5A-MR', 'Atlantic': [-55,-50,27.3,28],    'Pacific': [-63,-55,27.4,28], 'Indian': [-70,-50,27.2,28]}, # 15
            {'name':'IPSL-CM5B-LR', 'Atlantic': [-65,-50,27.2,27.9],    'Pacific': [-70,-53,27,27.8], 'Indian': [-68,-58,27.4,27.9]} # 16
        ]


    for imodel in range(len(domains)):
        if domains[imodel]['name'] == model_name :
            varout = domains[imodel]

    return varout
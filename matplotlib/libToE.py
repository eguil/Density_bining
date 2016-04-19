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
    #print ' ToE CPU = ',tcpu1-tcpu0

    return toe

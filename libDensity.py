'''
 libDensity.py contains all the functions needed by binDensity.py
'''


import cdms2 as cdm
import cdutil as cdu
import MV2 as mv
import numpy as npy
from string import replace
from genutil import statistics
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
    - valmask    - 1D scalar of mask value
    
    Output:
    - field     - 1D/2D/3D masked array
    
    Usage:
    ------
    >> from libDensity import maskVal
    >> maskedVariable = maskVal(unMaskedVariable,valmask)

    Notes:
    -----
    - PJD 15 Sep 2014 - 
    '''
    field [npy.isnan(field.data)] = valmask
    field._FillValue = valmask
    field = mv.masked_where(field > valmask/10, field)
    return field


def binIndices(szm,s_s):
    '''
    binIndices returns indices of z points in density bins
    Created to replace original interp loop in binDensity.py


    Author:    Eric Guilyardi : Eric.Guilyardi@locean-ipsl.upmc.fr
    Co-author:

    Created on Tue Sep  4 14:49:45 CEST 2018

    Inputs:
    ------
    - szm     - density as f(z,dim2) - 2D masked array (dim2 can be space or space x time)
    - s_s     - 1D target density grid

    Output:
    -------
    - binIndex - 3D masked array (density, depth, dim2)

    Usage:
    ------
    >> from libDensity import binIndices
    >> binIndices = binIndices(szm,s_s)

    Notes:
    -----
    TODO: Needs to fill in missing point within density grid
    '''
    print 'Enter binIndices'
    tcpu = timc.clock()
    N_z = szm.shape[0]
    dim2   = szm.shape[1]
    N_s   = len(s_s)
    print '  -> number of vertical z levels ', N_z
    print '  -> number of density s levels  ', N_s
    print '  -> dim2 ',dim2

    binIndex = npy.ma.ones([N_s, N_z, dim2], dtype='int')*0

    # Compute signs of [szm(k)-s_s(l)] and [szm(k)-s_s(l+1)]
    # if product of two is negative then szm(k) is in bin [s_s(l)-s_s(l+1)]

    # replicate s_s to have dim2
    print npy.tile(s_s, dim2).shape
    print npy.tile(s_s, dim2)
    ssr = npy.tile(s_s, dim2).reshape(dim2,N_s).transpose()
    print '  -> size of ssr', ssr.shape
    print ssr
    # Loop on k (s index)
    for l in range(N_s-1):
        print ' ====== l ======', l
        # replicate s_s(l) and s_s(l+1) on all k levels
        rszk   = npy.tile(ssr[l  ,:], N_z).reshape(N_z,dim2)
        rszkp1 = npy.tile(ssr[l+1,:], N_z).reshape(N_z,dim2)
        #print '  -> size of rszk, rszkp1',rszk.shape, rszkp1.shape
        print rszk
        print rszkp1
        product = npy.sign((rszk - szm)*(rszkp1 - szm))
        print product
        indneg = npy.argwhere (product <= 0).transpose()
        binIndex[l,indneg[0], indneg[1]] = 1

    return binIndex

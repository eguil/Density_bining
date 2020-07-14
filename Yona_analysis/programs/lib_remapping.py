#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from scipy.interpolate import griddata

def remaptoz(fieldr,depthr,targetz):
    """
    From field(sigma) and pseudo-z(sigma), build field(pseudo-z)
    
    input: fieldr - (basin,density,latitude)
           depthr - reference pseudo-z (basin,density,lat)
           targetz - target 1D depth grid for remapping
           
    output: fieldz - (basin,pseudo-z,latitude)
            zbowl - (basin,latitude)
    
    June 2019: correction for values to correspond to correct z levels
    March 2020: adding zbowl output
    """
    
    basinN = fieldr.shape[0]
    densityN = fieldr.shape[1]
    latN = fieldr.shape[2]
    
    fieldz = np.ma.masked_all((basinN,len(targetz),latN))
    zbowl = np.zeros((basinN,latN))
    
    for ibasin in range(basinN):
        for ilat in range(latN):
            
            zsig = depthr[ibasin,:,ilat] # Read pseudo-depth (function of sigma) of water column
            field_sig = fieldr[ibasin,:,ilat] # Read field values of the water column

            field_sort = np.ma.compressed(field_sig) # Remove masked values
            zsort = zsig[np.where(field_sig!=np.ma.masked)]

            if len(zsort) > 1:
                zroll = np.roll(zsort,1)
                zmean=(zsort+zroll)/2
                zmean[0] = zsort[0]/2
                fieldz[ibasin,:,ilat] = griddata(zmean,field_sort,targetz) # Grid field with target pressure grid at correct z levels
                zbowl[ibasin,ilat] = zmean[0]
            else:
                fieldz[ibasin,:,ilat] = np.ma.masked
                zbowl[ibasin,ilat] = np.ma.masked

    # Mask nans
    fieldz[np.isnan(fieldz)] = np.ma.masked
    
    return fieldz, zbowl


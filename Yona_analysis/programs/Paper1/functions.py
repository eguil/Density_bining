import xarray as xr
import numpy as np
from scipy.interpolate import griddata

def remaptoz_xr(fieldr,depthr,targetz):

    basinN = fieldr.shape[0]
    densityN = fieldr.shape[1]
    latN = fieldr.shape[2]

    fieldz = xr.DataArray(np.zeros((basinN,len(targetz),latN)),dims=['basin','pseudo-depth','latitude'],
                          coords=[fieldr.basin,targetz,fieldr.latitude])
    zbowl = np.zeros((basinN,latN))

    for ibasin in range(basinN):
        for ilat in range(latN):

            zsig = depthr[ibasin,:,ilat] # Read pseudo-depth (function of sigma) of water column
            field_sig = fieldr[ibasin,:,ilat] # Read field values of the water column

            field_sort = field_sig[~np.isnan(field_sig)] # Remove nans
            zsort = zsig[~np.isnan(field_sig)]

            if len(zsort) > 1:
                zroll = np.roll(zsort,1)
                zmean=(zsort+zroll)/2
                zmean[0] = zsort[0]/2
                fieldz[ibasin,:,ilat] = griddata(zmean,field_sort,targetz) # Grid field with target pressure grid at correct z levels
                # Convert bowl density to pseudo-depth = first non-nan element
                zbowl[ibasin,ilat] = zmean[0]
            else:
                fieldz[ibasin,:,ilat] = np.nan
                zbowl[ibasin,ilat] = np.nan

    zbowl = xr.DataArray(zbowl,dims=['basin','latitude'],
                          coords=[fieldr.basin,fieldr.latitude])
    return fieldz, zbowl


def lag_linregress_3D(x, y, lagx=0, lagy=0):
    """
    Input: Two xr.Datarrays of any dimensions with the first dim being time. 
    Thus the input data could be a 1D time series, or for example, have three 
    dimensions (time,lat,lon). 
    Datasets can be provided in any order, but note that the regression slope 
    and intercept will be calculated for y with respect to x.
    Output: Covariance, correlation, regression slope and intercept, p-value, 
    and standard error on regression between the two datasets along their 
    aligned time dimension.  
    Lag values can be assigned to either of the data, with lagx shifting x, and
    lagy shifting y, with the specified lag amount. 
    """ 
    #1. Ensure that the data are properly alinged to each other. 
    x,y = xr.align(x,y)

    #2. Add lag information if any, and shift the data accordingly
    if lagx!=0:

        # If x lags y by 1, x must be shifted 1 step backwards. 
        # But as the 'zero-th' value is nonexistant, xr assigns it as invalid 
        # (nan). Hence it needs to be dropped
        x   = x.shift(time = -lagx).dropna(dim='time')

        # Next important step is to re-align the two datasets so that y adjusts
        # to the changed coordinates of x
        x,y = xr.align(x,y)

    if lagy!=0:
        y   = y.shift(time = -lagy).dropna(dim='time')
        x,y = xr.align(x,y)

    #3. Compute data length, mean and standard deviation along time axis: 
    n = y.notnull().sum(dim='time')
    xmean = x.mean(axis=0)
    ymean = y.mean(axis=0)
    xstd  = x.std(axis=0)
    ystd  = y.std(axis=0)

    #4. Compute covariance along time axis
    cov   =  np.sum((x - xmean)*(y - ymean), axis=0)/(n)

    #5. Compute correlation along time axis
    cor   = cov/(xstd*ystd)

    #6. Compute regression slope and intercept:
    slope     = cov/(xstd**2)
    intercept = ymean - xmean*slope  

    #7. Compute P-value and standard error
    #Compute t-statistics
    tstats = cor*np.sqrt(n-2)/np.sqrt(1-cor**2)
    stderr = slope/tstats

    from scipy.stats import t
    pval   = t.sf(tstats, n-2)*2
    pval   = xr.DataArray(pval, dims=cor.dims, coords=cor.coords)

    return cov,cor,slope,intercept,pval,stderr
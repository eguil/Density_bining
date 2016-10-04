"""
Densitlib for matplotlib for density plots
(c) Eric Guilyardi
Feb 2016

"""

import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import MaxNLocator


# --------------------------------------------
#  Build zonal mean with zoom on density axis
# --------------------------------------------

def zon_2dom(plt, ax0, ax1, lat, lev, varBasin, varSigma, unit, minmax, clevsm, cmap, domrho, agreelev, agreeON, noax,
             labBowl):
    # -- variables
    var = varBasin['diffBowl']
    varm = varBasin['meanBowl']
    varag = varBasin['agree']

    # -- title and bowl labels
    title = varBasin['name']
    label1 = labBowl[0]
    label2 = labBowl[1]

    # -- contour levels
    rhomin = domrho[0]
    rhomid = domrho[1]
    rhomax = domrho[2]

    # -- Latmin/max
    latmin = -80.
    latmax = 80.
    deltalat = 20
    labels = ['', '60S', '40S', '20S', 'Eq', '20N', '40N', '60N', '']
    #
    # ====  Upper panel  ===================================================
    #
    ax0.axis([latmin, latmax, rhomin, rhomid])
    ax0.invert_yaxis()
    ax0.xaxis.set_ticks(np.arange(latmin, latmax, deltalat))
    ax0.set_xticklabels(labels)
    if noax == 'T':
        ax0.set_yticklabels([])
        ax1.set_yticklabels([])
    if noax == 'R':
        ax0.yaxis.tick_right()
        ax1.yaxis.tick_right()

    # -- levels for shade plot
    levels = MaxNLocator(nbins=minmax[2]).tick_values(minmax[0], minmax[1])

    # -- Format for contour labels
    levfmt='%.0f'
    if abs(clevsm[1]-clevsm[0]) < 1:
        levfmt='%.1f'
    if abs(clevsm[1]-clevsm[0]) < 0.1:
        levfmt='%.2f'

    # -- draw filled contours of period diff
    cnplot = ax0.contourf(lat, lev, var, cmap=cmap, levels=levels)

    if agreeON:
        # -- draw agreement contour > agreement level (agreelev)
        cmapbl = LinearSegmentedColormap('cmapbl', bluecol())
        chplot = ax0.contourf(lat, lev, varag, levels=[-agreelev, agreelev], hatches=['..'], colors='none')
        cpplot = ax0.contour(lat, lev, varag, [agreelev - .0001, agreelev + 0.00001], cmap=cmapbl, linewidths=2)
        cpplot = ax0.contour(lat, lev, varag, [-agreelev - .0001, -agreelev + 0.00001], cmap=cmapbl, linewidths=2)

    # -- draw mean contours
    cmapb = LinearSegmentedColormap('cmapb', blkcol())
    cpplot = ax0.contour(lat, lev, varm, clevsm, cmap=cmapb)
    ax0.clabel(cpplot, inline=1, fontsize=10, fmt=levfmt)

    # -- draw ptopsigma for 2 periods (yr1 = ref, yr2 = end of serie)
    lnplot1 = ax0.plot(lat, varSigma['yr1'], linestyle='--', color='black', linewidth=2)
    lnplot2 = ax0.plot(lat, varSigma['yr2'], linestyle='-', color='black', linewidth=2)

    #
    # ====  Lower panel   ===================================================
    #
    ax1.axis([latmin, latmax, rhomid, rhomax])
    ax1.invert_yaxis()
    ax1.xaxis.set_ticks(np.arange(latmin, latmax, deltalat))
    ax1.set_xticklabels(labels)

    # -- draw filled contours
    cnplot = ax1.contourf(lat, lev, var, cmap=cmap, levels=levels)

    if agreeON:
        # -- draw agreement contour > agreement level (agreelev)
        cmapbl = LinearSegmentedColormap('cmapbl', bluecol())
        chplot = ax1.contourf(lat, lev, varag, levels=[-agreelev, agreelev], hatches=['..'], colors='none')
        cpplot = ax1.contour(lat, lev, varag, [agreelev - .0001, agreelev + 0.00001], cmap=cmapbl, linewidths=2)
        cpplot = ax1.contour(lat, lev, varag, [-agreelev - .0001, -agreelev + 0.00001], cmap=cmapbl, linewidths=2)

    # -- draw mean contours
    cmapb = LinearSegmentedColormap('cmapb', blkcol())
    cpplot = ax1.contour(lat, lev, varm, clevsm, cmap=cmapb)
    ax1.clabel(cpplot, inline=1, fontsize=10, fmt=levfmt)

    # -- draw ptopsigma for 2 periods (yr1 = ref, yr2 = end of serie)
    lnplot1b = ax1.plot(lat, varSigma['yr1'], linestyle='--', color='black', linewidth=2, label=label1)
    lnplot2b = ax1.plot(lat, varSigma['yr2'], linestyle='-', color='black', linewidth=2, label=label2)

    # -- Add legend for bowl position
    plt.legend(loc='upper right', title='Bowl')

    # -- add plot title
    ax0.text(latmin + 10, rhomin + 1, title, fontsize=14, fontweight='bold')

    return [cnplot, lnplot1, lnplot2]


# --------------------------------------------
#  Build zonal mean with zoom on depth axis
# --------------------------------------------

def zon_2domz(plt, ax0, ax1, lat, lev, var, title, minmax, clevsm, cmap, domzed, noax):
    # -- variables

    # -- title and bowl labels

    # -- contour levels
    zedmin = domzed[0]
    zedmid = domzed[1]
    zedmax = domzed[2]

    # -- Latmin/max
    latmin = -80.
    latmax = 80.
    deltalat = 20
    labels = ['', '60S', '40S', '20S', 'Eq', '20N', '40N', '60N', '']
    #
    # ====  Upper panel  ===================================================
    #
    ax0.axis([latmin, latmax, zedmin, zedmid])
    ax0.invert_yaxis()
    ax0.xaxis.set_ticks(np.arange(latmin, latmax, deltalat))
    ax0.set_xticklabels(labels)
    if noax == 'T':
        ax0.set_yticklabels([])
        ax1.set_yticklabels([])
    if noax == 'R':
        ax0.yaxis.tick_right()
        ax1.yaxis.tick_right()

    # -- levels for shade plot
    levels = MaxNLocator(nbins=minmax[2]).tick_values(minmax[0], minmax[1])

    # -- Format for contour labels
    levfmt='%.0f'
    if abs(clevsm[1]-clevsm[0]) < 1:
        levfmt='%.1f'
    if abs(clevsm[1]-clevsm[0]) < 0.1:
        levfmt='%.2f'

    # -- draw filled contours of period diff
    cnplot = ax0.contourf(lat, lev, var, cmap=cmap, levels=levels)

    # -- draw mean contours
    cmapb = LinearSegmentedColormap('cmapb', blkcol())
    cpplot = ax0.contour(lat, lev, var, clevsm, cmap=cmapb)
    ax0.clabel(cpplot, inline=1, fontsize=10, fmt=levfmt)

    #
    # ====  Lower panel   ===================================================
    #
    ax1.axis([latmin, latmax, zedmid, zedmax])
    ax1.invert_yaxis()
    ax1.xaxis.set_ticks(np.arange(latmin, latmax, deltalat))
    ax1.set_xticklabels(labels)

    # -- draw filled contours
    cnplot = ax1.contourf(lat, lev, var, cmap=cmap, levels=levels)

    # -- draw mean contours
    cmapb = LinearSegmentedColormap('cmapb', blkcol())
    cpplot = ax1.contour(lat, lev, var, clevsm, cmap=cmapb)
    ax1.clabel(cpplot, inline=1, fontsize=10, fmt=levfmt)

    # -- Add legend for bowl position
    plt.legend(loc='upper right', title='Bowl')

    # -- add plot title
    ax0.text(latmin + 10, zedmin + 1, title, fontsize=14, fontweight='bold')

    return [cnplot]


# --------------------------------
#   List variables and properties
# --------------------------------

def defVar(longName):
    salinity = {
        'longN': 'salinity',  # long name
        'var': 'isonso',  # variable name
        'minmax': [-.2, .2, 16],  # for diff shading + number of color interval
        'clevsm': np.arange(30, 40, .2),  # for mean contours
        'clevsmdif': np.arange(-.2, .2, .025),  # for mean contours
        'clevsmstd': np.arange(0., .2, .005),  # for stddev contours
        '1dminmax': [-.1, .1], # for 1D ToE plots
        'legVar': "Salinity",  # Legend name
        'unit': "PSU",  # TODO: could be read from file
    }

    temp = {'var': 'isonthetao', 'minmax': [-.6, .6, 16], 'clevsm': np.arange(-2, 30, 1),
            'clevsmstd': np.arange(0, 2., .01), '1dminmax': [-.4, .4],'clevsmdif': np.arange(-.4, .4, .05),
            'legVar': "Temperature", 'unit': "C", 'longN': 'temp',
            }
    depth = {'var': 'isondepth', 'minmax': [-75., 75., 10], 'clevsm': np.arange(0, 2000, 100),
             'clevsmstd': np.arange(0, 20, 5),'1dminmax': [-10, 50],'clevsmdif': np.arange(-75, 75, 10),
             'legVar': "Depth", 'unit': "m", 'longN': 'depth',
             }
    volume = {'var': 'isonvol', 'minmax': [-20., 20., 20], 'clevsm': np.arange(0, 200, 20),
              'clevsmstd': np.arange(0, 20, 1),'1dminmax': [-5, 5],'clevsmdif': np.arange(-20, 20, 5),
              'legVar': "Volume", 'unit': "1.e12 m^3", 'longN': 'volume',
              }
    persist = {'var': 'isonpers', 'minmax': [-10., 10., 20], 'clevsm': np.arange(0, 90, 10),
               'clevsmstd': np.arange(0, 3., .5),'1dminmax': [-50, 50],'clevsmdif': np.arange(-10, 10, 2),
               'legVar': "Persistence", 'unit': "% of time", 'longN': 'persist'
               }
    heatcontent = {'var': 'isonhtc', 'minmax': [-10., 10., 20], 'clevsm': np.arange(0, 90, 10),
               'clevsmstd': np.arange(0, 3., .5),'1dminmax': [-50, 50],'clevsmdif': np.arange(-10, 10, 2),
               'legVar': "Heat content", 'unit': "10^XX J", 'longN': 'persist'
               }

    vars = [salinity, temp, depth, volume, persist, heatcontent]

    varout = 'None'
    for ivar in range(len(vars)):
        if vars[ivar]['longN'] == longName:
            varout = vars[ivar]

    return varout


# - def black color bar
def blkcol():
    cdict = {'red': ((0.0, 0.0, 0.0),
                     (1.0, 1.0, 1.0)),
             'green': ((0.0, 0.0, 0.0),
                       (1.0, 1.0, 1.0)),
             'blue': ((0.0, 0.0, 0.0),
                      (1.0, 1.0, 1.0))}
    return (cdict)


# - def blues color bar
def bluecol():
    cdict = {'red': ((0.0, 0.0, 0.0),
                     (1.0, 0.0, 0.0)),
             'green': ((0.0, 0.0, 0.0),
                       (1.0, 0.0, 0.0)),
             'blue': ((0.0, 0.0, 1.0),
                      (1.0, 1.0, 1.0))}
    return (cdict)

# - average in lat/rho domain

def averageDom(field, dim, domain, lat, rho):

    latidx = np.argwhere((lat >= domain[0]) & (lat <= domain[1])).transpose()
    rhoidx = np.argwhere((rho >= domain[2]) & (rho <= domain[3])).transpose()
    lidx1 = latidx[0][0];
    lidx2 = latidx[0][-1]
    ridx1 = rhoidx[0][0];
    ridx2 = rhoidx[0][-1]
    if dim == 3:
        vara = np.ma.average(field[:, ridx1:ridx2, :], axis=1)
        var_ave = np.ma.average(vara[:, lidx1:lidx2], axis=1)
    else:
        vara = np.ma.average(field[ridx1:ridx2, :], axis=0)
        var_ave = np.ma.average(vara[lidx1:lidx2], axis=0)

    return var_ave

def remapToZ(fieldr,depthr,volumr, valmask, targetz):
    '''
    The remaToZ() function remaps a density bined zonal field back to z
    It starts from the surface and computes the mean field for each z level, using the zonal volume of isopycnals for weighting

    Author:    Eric Guilyardi : Eric.Guilyardi@locean-ipsl.upmc.fr

    Created on Thu Aug 18 09:50:29 CEST 2016

    Inputs:
    ------
    - fieldr     - field (T or S)       - 3D density,latitude,time array
    - depthr     - depth of isopycnals  - 3D density,latitude,time array
    - volumr     - volume of isopycnals - 3D density,latitude,time array
    - valmaks    - value of masked points
    - targetz    - target z grid        - 1D

    Output:
    - fieldz    - 3D depth,latitude,time array

    Usage:
    ------

    Notes:
    -----
    - EG 18 Aug 2016   - Initial function write
   '''

    latN = fieldr.shape[3]
    rhoN = fieldr.shape[2]
    basN = fieldr.shape[1]
    timN = fieldr.shape[0]

    levN = len(targetz)

    # init fieldx array
    fieldz  = np.ma.ones([timN, basN, levN, latN], dtype='float32')*valmask
    fieldz_vol = np.ma.ones([levN], dtype='float32')*valmask
    volumez = np.ma.ones([levN], dtype='float32')*valmask
    # Method 1: loop on vertical levels of target grid - may lead to levels missing data
    method = 1
    if method == 1:
        #for t in range(timN):
        for t in range(1):
            for j in range(latN):
                for k in range(levN-1):
                    field_int = 0.
                    volum_int = 0.
                    for r in range(rhoN):
                        if volumr[t,0,r,j] <> 0:
                            if depthr[t,0,r,j] >= targetz[k] and depthr[t,0,r,j] < targetz[k+1]:
                                field_int = field_int + fieldr[t,0,r,j]*volumr[t,0,r,j]
                                volum_int = volum_int + volumr[t,0,r,j]
                                depthr[t,0,r,j] = -100. # to speed up search for next depths
                    if volum_int <> 0.:
                        fieldz[t,0,k,j] = field_int / volum_int
    # Method 2: interpolation (reverse of what was done in binDensity.py:)
    # interpolate depth(z) (=z_zt=zzm) to depth(s) at s_s densities (=z_s) using density(z) (=s_z=szm)
    # for i in range(lonN*latN):
    #     if nomask[i]:
    #         z_s [0:N_s,i] = npy.interp(s_s[:,i], szm[:,i], zzm[:,i], right = valmask) ; # depth - consider spline
    #         c1_s[0:N_s,i] = npy.interp(z_s[0:N_s,i], zzm[:,i], c1m[:,i], right = valmask) ; # thetao
    #         c2_s[0:N_s,i] = npy.interp(z_s[0:N_s,i], zzm[:,i], c2m[:,i], right = valmask) ; # so
    if method == 2:
        #for t in range(timN):
        for t in range(1):
            for j in range(latN-1):
                print 'j=',j
                nomask = np.argwhere(volumr[t,0,:,j] < valmask/10.)
                print nomask
                if len(nomask) >= 2:
                    print depthr[t,0,nomask,j], fieldr[t,0,nomask,j]
                    print 'vol:',volumr[t,0,nomask,j]
                    depr = depthr[t,0,nomask,j]
                    fxvol = fieldr[t,0,nomask,j]*volumr[t,0,nomask,j]
                    volr = volumr[t,0,nomask,j]
                    print len(fxvol), len(depr), len(targetz)
                    fieldz_vol[:] = np.interp(targetz,depr,fxvol, right=valmask)
                    volumez[:]    = np.interp(targetz,depr,volr, right=valmask)
                    print 'fld_*vol_interp:',fieldz_vol
                    print 'volumez:', volumez
                    for k in range(levN):
                        if volumez[k] <> 0.:
                            fieldz[t,0,k,j] = fieldz_vol[k]/volumez[k]
                    print fieldz[t,0,:,j]






    return fieldz

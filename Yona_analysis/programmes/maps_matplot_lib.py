#!/bin/env python
# -*- coding: utf-8 -*-


import numpy as np
import numpy.ma as ma
from mpl_toolkits.basemap import Basemap
#from matplotlib.ticker import MaxNLocator
from netCDF4 import Dataset as open_ncfile
from matplotlib import cbook
from matplotlib.colors import Normalize

# --------------------------------
#   List variables and properties
# --------------------------------

def defVar(longName):
    salinity = {
        'longN': 'salinity',  # long name
        'var': 'sog',  # variable name
        'minmax': [-0.3, 0.3, 16],  # for diff shading + number of color interval
        'clevsm': np.arange(30, 40, .4),  # for mean contours
        'clevsm_bold': np.arange(0, 30, 5),
        'clevsmdif': np.arange(-.2, .2, .025),  # for mean contours
        'clevsmstd': np.arange(0., .2, .005),  # for stddev contours
        '1dminmax': [-.1, .1], # for 1D ToE plots
        'legVar': "Salinity",  # Legend name
        'unit': "PSU",  # TODO: could be read from file
    }

    temp = {'var': 'thetaog', 'minmax': [-0.65, 0.65, 14], 'clevsm': np.arange(-2, 30, 1), 'clevsm_bold': np.arange(0,30,5),
            'clevsmstd': np.arange(0, 2., .01), '1dminmax': [-.4, .4],'clevsmdif': np.arange(-.4, .4, .05),
            'legVar': "Temperature", 'unit': "C", 'longN': 'temp',
            }
    depth = {'var': 'isondepthg', 'minmax': [-75., 75., 10], 'clevsm': np.arange(0, 2000, 100),
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
               'legVar': "Heat content", 'unit': "10^XX J", 'longN': 'heatcontent'
               }

    vars = [salinity, temp, depth, volume, persist, heatcontent]

    varout = 'None'
    for ivar in range(len(vars)):
        if vars[ivar]['longN'] == longName:
            varout = vars[ivar]

    return varout


def defVarDurack(longName):
    salinity = {'var_change': 'salinity_change', 'var_change_er':'thetao_change_error',
                'var_mean': 'salinity_mean',
                'var_mean_zonal': 'salinity_mean_basin_zonal',
                'var_change_zonal': 'salinity_change_basin_zonal', 'var_change_zonal_er' : 'salinity_change_error_basin_zonal',
                'minmax': [-0.3, 0.3, 16],
                'minmax_zonal': [-0.2, 0.2, 16],
                'clevsm': np.arange(30, 40, .25),
                'clevsm_zonal': np.arange(30, 40, .1),
                'clevsm_bold': np.arange(30, 40, .5),
                'legVar': "Salinity", 'unit': "PSU", 'longN': 'salinity'}

    temp = {'var_change':'thetao_change', 'var_change_er':'thetao_change_error',
            'var_mean':'thetao_mean',
            'var_mean_zonal': 'thetao_mean_basin_zonal',
            'var_change_zonal': 'thetao_change_basin_zonal', 'var_change_zonal_er' : 'thetao_change_error_basin_zonal',
            'minmax': [-0.65, 0.65, 14],
            'minmax_zonal' : [-0.4,0.4,16],
            'clevsm': np.arange(-2, 30, 1),
            'clevsm_zonal': np.arange(-2, 30, 1),
            'clevsm_bold': np.arange(-2,30,2),
            'legVar': "Temperature", 'unit': "C", 'longN': 'temp'}

    vars = [salinity,temp]
    for ivar in range(len(vars)):
        if vars[ivar]['longN'] == longName:
            varout = vars[ivar]

    return varout


def defVarmme(longName):
    salinity = {'var_zonal': 'isonsoBowl',
                'var_global': 'sogBowl', 'var_global_std':'sogBowlStd',
                'minmax': [-0.3, 0.3, 16],
                'minmax_zonal': [-0.2, 0.2, 16],
                'clevsm': np.arange(30, 40, .25),
                'clevsm_zonal': np.arange(30, 40, .1),
                'clevsm_bold': np.arange(30, 40, .5),
                'legVar': "Salinity", 'unit': "PSU", 'longN': 'salinity'}

    temp = {'var_zonal':'isonthetaoBowl',
            'var_global': 'thetaogBowl', 'var_global_std':'thetaogBowlStd',
            'minmax': [-0.65, 0.65, 14],
            'minmax_zonal' : [-0.4,0.4,16],
            'clevsm': np.arange(-2, 30, 1),
            'clevsm_zonal': np.arange(-2, 30, 1),
            'clevsm_bold': np.arange(-2,30,2),
            'legVar': "Temperature", 'unit': "C", 'longN': 'temp'}

    vars = [salinity,temp]
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


# Customed colormap
def custom_div_cmap(numcolors=17, name='custom_div_cmap',
                    mincol='blue', midcol='white', maxcol='red'):
    """ Create a custom diverging colormap with three colors

    Default is blue to white to red with 17 colors.  Colors can be specified
    in any way understandable by matplotlib.colors.ColorConverter.to_rgb()
    """

    from matplotlib.colors import LinearSegmentedColormap

    cmap = LinearSegmentedColormap.from_list(name=name,
                                             colors=[mincol, midcol, maxcol],
                                             N=numcolors)
    return cmap


# ----------------------------------------------------
#   Build zonal latitude/density plot
# ----------------------------------------------------

def zonal_2D(plt, action, ax0, ax1, ticks, lat, density, varBasin, minmax, domrho, clevsm=None, clevsm_bold=None):

    # Colormap
    cmap = custom_div_cmap()

    # latitude domain
    domlat = [-70, 70]


    if action == 'total' :
        var = varBasin['var_change']
        var_mean = varBasin['var_mean']
        var_er = varBasin['var_error']
        levels = np.linspace(minmax[0], minmax[1], minmax[2])

        # # -- Error field
        # var_er = var_er * 1.1  # to account for a potential underestimation of the error determined by a bootstrap analysis
        # var_er = var_er * 1.64  # 90% confidence level
        # not_signif_change = np.where(np.absolute(var) < var_er, 1, 0)
        #
        # -- Format for contour labels
        levfmt = '%.0f'
        if abs(clevsm[1] - clevsm[0]) < 1:
            levfmt = '%.1f'
        if abs(clevsm[1] - clevsm[0]) < 0.1:
            levfmt = '%.2f'

    if action == 'total_mme':
        var = varBasin['var_change']
        var_mean = varBasin['var_mean']
        levels = np.linspace(minmax[0], minmax[1], minmax[2])

        # -- Format for contour labels
        levfmt = '%.0f'
        if abs(clevsm[1] - clevsm[0]) < 1:
            levfmt = '%.1f'
        if abs(clevsm[1] - clevsm[0]) < 0.1:
            levfmt = '%.2f'

    if action == 'isopyc_mig':
        var = varBasin['dvar_dsig']*varBasin['dsig_dt'] + varBasin['dvar_dy']*varBasin['dy_dt']
        levels = np.linspace(minmax[0], minmax[1], minmax[2])

    if action == 'isopyc_mig_lat':
        var = varBasin['dvar_dy']*varBasin['dy_dt']
        levels = np.linspace(minmax[0], minmax[1], minmax[2])
        cmap = custom_div_cmap()

    if action == 'isopyc_mig_sig':
        var = varBasin['dvar_dsig']*varBasin['dsig_dt']
        levels = np.linspace(minmax[0], minmax[1], minmax[2])
        cmap = custom_div_cmap()

    if action == 'dvar_dsig':
        var = varBasin['dvar_dsig']
        levels = np.arange(-2,2.01,0.25)
        cmap = plt.get_cmap('bwr')

    if action == 'dvar_dy':
        var = varBasin['dvar_dy']
        levels = np.arange(-0.1,0.101,0.01)
        cmap = plt.get_cmap('bwr')

    if action == 'dsig_dt':
        var = varBasin['dsig_dt']
        levels = np.linspace(-0.5,0.5,16)
        cmap = custom_div_cmap()

    if action == 'dy_dt':
        var = varBasin['dy_dt']
        levels = np.linspace(-5, 5, 16)
        cmap = custom_div_cmap()

    if action == 'residual':
        var = varBasin['var_change_res']
        levels = np.linspace(minmax[0], minmax[1], minmax[2])

    if action == 'var_2000_hr':
        var = varBasin['var_2000_hr']
        levels = np.arange(33.75,36,0.25)

    if action == 'var_2000_sig_hr':
        var = varBasin['var_2000_sig_hr']
        levels = None

    if action == 'min_dist_lat':
        var = varBasin['min_dist_lat']
        levels = np.arange(-10,11,2)
        cmap = plt.get_cmap('bwr')

    if action == 'min_dist_sig':
        var = varBasin['min_dist_sig']
        levels = np.arange(-1,1.01,0.25)
        cmap = plt.get_cmap('bwr')

    if action == 'delta_var_min':
        var = varBasin['delta_var_min']
        levels = np.arange(0,0.0301,0.005)
        cmap = plt.get_cmap('jet')

    if action == 'mean_fields':
        var_1950 = varBasin['var_1950']
        var_2000 = varBasin['var_2000']
        levels = clevsm

        # -- Format for contour labels
        levfmt = '%.0f'
        if abs(clevsm[1] - clevsm[0]) < 1:
            levfmt = '%.1f'
        if abs(clevsm[1] - clevsm[0]) < 0.1:
            levfmt = '%.2f'

    if action != 'total' and action != 'total_mme':
        bowl = varBasin['bowl']

    # Create meshgrid
    lat2d, density2d = np.meshgrid(lat, density)

    # ==== Upper panel ====

    if action == 'var_2000_hr' or action == 'var_2000_sig_hr':
        cnplot1 = ax0.contourf(lat2d, density2d, var, levels, cmap=plt.get_cmap('jet'), extend='both')
    elif action == 'mean_fields':
        ax0.contour(lat2d,density2d,var_1950, levels=levels, colors='black', linewidths=0.5)
        cpplot11 = ax0.contour(lat2d,density2d,var_1950, levels=clevsm_bold, colors='black', linewidths=1.5)
        ax0.clabel(cpplot11, inline=1, fontsize=11, fmt=levfmt)
        ax0.contour(lat2d,density2d,var_2000, levels=levels, colors='black', linewidths=0.5, linestyles='dashed')
        cpplot12 = ax0.contour(lat2d,density2d,var_2000, levels=clevsm_bold, colors='black', linewidths=1.5, linestyles='dashed')
        #ax0.clabel(cpplot12, inline=1, fontsize=11, fmt=levfmt)
    else:
        cnplot1 = ax0.contourf(lat2d, density2d, var, cmap=cmap, levels=levels, extend='both')

    # if (action == 'total_mme' and var_mean != None) or action == 'total' :
    #     cpplot11 = ax0.contour(lat2d, density2d, var_mean, clevsm, colors='black', linewidths=0.5)
    #     #ax0.clabel(cpplot11, inline=1, fontsize=10, fmt=levfmt)
    #     cpplot12 = ax0.contour(lat2d, density2d, var_mean, clevsm_bold, colors='black', linewidths=2)
    #     ax0.clabel(cpplot12, inline=1, fontsize=12, fontweight='bold', fmt=levfmt)

        # error_plot = ax0.contourf(lat2d, density2d, not_signif_change, levels=[0.25,0.5,1.5], colors='None',
        #                            hatches=['','....'], edgecolor='0.3', linewidth=0.0)

    if action != 'var_2000_hr' and action != 'var_2000_sig_hr' and action!='total' and action!='total_mme':
        ax0.plot(lat, bowl, color='black')

    ax0.set_ylim([domrho[0], domrho[1]])
    ax0.set_xlim([domlat[0], domlat[1]])
    ax0.invert_yaxis()
    ax0.tick_params(
        axis='x',  # changes apply to the x axis
        which='both',  # both major and minor ticks are affected
        bottom='off',  # ticks along the bottom edge are off
        labelbottom='off',
        top='off')

    if ticks != 'left':
        ax0.tick_params(axis='y', labelleft='off')
    if ticks == 'right':
        ax0.tick_params(axis='y', labelright='on')

    ax0.axvline(x=0, color='black', ls='--')


    # === Lower panel ====

    if action == 'var_2000_hr' or action == 'var_2000_sig_hr':
        cnplot2 = ax1.contourf(lat2d, density2d, var, levels, cmap=plt.get_cmap('jet'), extend='both')
    elif action == 'mean_fields':
        ax1.contour(lat2d,density2d,var_1950, levels=levels, colors='black', linewidths=0.5)
        cpplot21 = ax1.contour(lat2d,density2d,var_1950, levels=clevsm_bold, colors='black', linewidths=1.5)
        ax1.clabel(cpplot21, inline=1, fontsize=11, fmt=levfmt)
        ax1.contour(lat2d,density2d,var_2000, levels=levels, colors='black', linewidths=0.5, linestyles='dashed')
        cpplot22 = ax1.contour(lat2d,density2d,var_2000, levels=clevsm_bold, colors='black', linewidths=1.5, linestyles='dashed')
        #ax1.clabel(cpplot22, inline=1, fontsize=11, fmt=levfmt)
        cnplot2 = cpplot22
    else:
        cnplot2 = ax1.contourf(lat2d, density2d, var, cmap=cmap, levels=levels, extend='both')

    # if (action == 'total_mme'  and var_mean != None) or action == 'total':
    #     cpplot21 = ax1.contour(lat2d, density2d, var_mean, clevsm, colors='black', linewidths=0.5)
    #     #ax1.clabel(cpplot21, inline=1, fontsize=10, fmt=levfmt)
    #     cpplot22 = ax1.contour(lat2d, density2d, var_mean, clevsm_bold, colors='black', linewidths=2)
    #     ax1.clabel(cpplot22, inline=1, fontsize=12, fontweight='bold', fmt=levfmt)

    #     # error_plot = ax1.contourf(lat2d, density2d, not_signif_change, levels=[0.25, 0.5, 1.5], colors='None',
    #     #                           hatches=['', '....'], edgecolor='0.6', linewidth=0.0)

    if action != 'var_2000_hr' and action != 'var_2000_sig_hr' and action!='total' and action!='total_mme':
        ax1.plot(lat, bowl, color='black')


    ax1.set_ylim([domrho[1], domrho[2]])
    ax1.set_xlim([domlat[0], domlat[1]])
    ax1.invert_yaxis()
    ax1.tick_params(
        axis='x',  # changes apply to the x axis
        which='both',  # both major and minor ticks are affected
        top='off')  # ticks along the bottom edge are off

    if ticks != 'left':
        ax1.tick_params(axis='y', labelleft='off')
    if ticks == 'right':
        ax1.tick_params(axis='y', labelright='on')


    ax1.axvline(x=0, color='black', ls='--')

    # Re-label x-axis
    xlabels = ['', '60S', '40S', '20S', '0', '20N', '40N', '60N']
    ax1.set_xticklabels(xlabels)

    # Remove intersecting tick at rhomid
    yticks = ax1.yaxis.get_major_ticks()
    if ticks == 'left':
        yticks[0].label1.set_visible(False)
    if ticks == 'right':
        yticks[0].label2.set_visible(False)


    # -- add plot title
    ax0.text(-60, 22, varBasin['name'], fontsize=14, fontweight='bold')


    return cnplot2, levels



# ----------------------------------------------------
#   Build Basemap of changes along a chosen isopycnal
# ----------------------------------------------------

def proj_map(kind, plt, ax, minmax, clevsm, clevsm_bold, lat, lon, cmap, isopyc, sliced_density, var1,
             var2 = None, var3 = None):

    isopyc_idx = np.argmin(np.abs(sliced_density - isopyc))

    if kind == 'hist-histNat':
        var_hist = np.squeeze(var1[:,isopyc_idx,:])
        var_histNat = np.squeeze(var2[:,isopyc_idx,:])
        # Difference
        var_diff = np.ma.average(var_hist[-5:, :], axis=0) - np.ma.average(var_histNat, axis=0)
        # Climatology
        var_mean = np.ma.average(var_hist, axis=0)

    elif kind == 'Durack':
        var_diff = np.squeeze(var1[isopyc_idx,:])
        var_mean = np.squeeze(var2[isopyc_idx,:])
        #Error field
        var_diff_er = np.squeeze(var3[isopyc_idx,:])
        var_diff_er = var_diff_er*1.1 # to account for a potential underestimation of the error determined by a bootstrap analysis
        var_diff_er = var_diff_er*2.58 # 99% level
        not_signif_change = np.where(np.absolute(var_diff)<var_diff_er, 1, 0)

    elif kind == 'hist':
        var = np.squeeze(var1[:,isopyc_idx,:,:])
        # look at difference between 1950 and end to compare with obs
        var_diff = np.ma.average(var[-5:,:], axis=0) - np.ma.average(var[0:5,:], axis=0)
        # Climatology
        var_mean = np.ma.average(var, axis=0)

    else :
        var_obs = np.squeeze(var1[:,isopyc_idx,:,:])
        var_diff = np.ma.average(var_obs[-5:,:], axis=0) - np.ma.average(var_obs[0:5,:], axis=0)
        # Climatology
        var_mean = np.ma.average(var_obs, axis=0)

    # Create meshgrid
    lon2d, lat2d = np.meshgrid(lon, lat)

    # Levels for shade plot
    #levels = MaxNLocator(nbins=minmax[2]).tick_values(minmax[0], minmax[1])
    levels = np.linspace(minmax[0],minmax[1],minmax[2])

    # Format for contour labels
    levfmt = '%.0f'
    if abs(clevsm[1] - clevsm[0]) < 1:
        levfmt = '%.1f'
    if abs(clevsm[1] - clevsm[0]) < 0.1:
        levfmt = '%.2f'

    # Read grid for coloring continents
    f2 = open_ncfile('/home/ysilvy/Density_bining/Yona_analysis/data/140807_WOD13_masks.nc', 'r')
    landsea = f2.variables['landsea'][:]

    # Create mask
    sea_mask = landsea != 1
    landsea = np.ma.array(landsea, mask=sea_mask)
    landsea[landsea == 1] = 0.2

    # Basemap
    map = Basemap(projection='cyl', llcrnrlon=20, llcrnrlat=-70., urcrnrlon=380, urcrnrlat=70, ax=ax)
    map.drawmapboundary(fill_color='1')
    map.drawparallels(np.arange(-60, 61, 20.), labels=[1, 0, 1, 0], linewidth=0.5)
    map.drawmeridians(np.arange(-180, 180, 60), labels=[0, 0, 0, 1], linewidth=0.5)

    # Draw filled contours of diff
    cnplot = map.contourf(lon2d, lat2d, var_diff, cmap=cmap, levels = levels, latlon=True, extend='both')

    if kind == 'Durack':
        map.fillcontinents(color='black')
        error_plot = map.contourf(lon2d, lat2d, not_signif_change, levels=[0.25,0.5,1.5], colors='None',
                               hatches=['','....'], edgecolor='0.6', linewidth=0.0, latlon=True)

    else :
        # Draw continents
        pc2 = map.pcolormesh(lon2d, lat2d, landsea, shading='flat', cmap=plt.cm.gray, latlon=True)

    # Draw mean contours
    cpplot1 = map.contour(lon2d, lat2d, var_mean, clevsm, colors = 'black', linewidths=0.5, latlon=True)
    #plt.clabel(cpplot1, inline=1, fontsize=10, fmt=levfmt)
    cpplot2 = map.contour(lon2d, lat2d, var_mean, clevsm_bold, colors='black', linewidths=2, latlon=True)
    plt.clabel(cpplot2, inline=1, fontsize=12, fontweight='bold', fmt=levfmt)

    # Indicate which isopycnal we're projecting on
    ax.text(0.1, 0.85, '$\gamma = %.1f$' %(isopyc,), transform=ax.transAxes, fontweight='bold',color='w', fontsize=17)

    return cnplot, levels




# -------------------------------------------------------------------------------------
#   Same function, but removing the zonal mean as well, and plotting it next to the map
# -------------------------------------------------------------------------------------

def proj_map_zonal_changes(kind, zonal_change, plt, ax1, ax2, minmax, clevsm, lat, lon, cmap, isopyc,
                           sliced_density, var1, var2 = None):


    isopyc_idx = np.argmin(np.abs(sliced_density - isopyc))

    if kind == 'model':
        var_hist = np.squeeze(var1[:,isopyc_idx,:])
        var_histNat = np.squeeze(var2[:,isopyc_idx,:])
        # Difference
        var_diff = np.ma.average(var_hist[-6:, :], axis=0) - np.ma.average(var_histNat[-6,:, :], axis=0)
        # Climatology
        var_mean = np.ma.average(var_hist, axis=0)

    elif kind == 'Durack':
        var_diff = np.squeeze(var1[isopyc_idx,:])
        var_mean = np.squeeze(var2[isopyc_idx,:])

    else:
        var_obs = np.squeeze(var1[:,isopyc_idx,:])
        var_diff = np.ma.average(var_obs[-6:,:], axis=0) - np.ma.average(var_obs[0:5,:], axis=0)
        # Climatology
        var_mean = np.ma.average(var_obs, axis=0)

    # Build zonal mean(s) and remove from var_diff
    zonal_mean = np.ma.average(var_diff, axis=1)
    zonal_mean_2D = np.tile(zonal_mean,(len(lon),1))
    zonal_mean_2D = np.transpose(zonal_mean_2D)
    if zonal_change == 'global':
        var_diff = var_diff - zonal_mean_2D
    if zonal_change == 'basin' and kind == 'Durack' :
        # Read basin mask
        f3 = open_ncfile('/home/ysilvy/data/DurackandWijffels_GlobalOceanSurfaceChanges_1950-2000_mask.nc', 'r')
        basin_mask = f3.variables['basin_mask'][:]
        # Build zonal means for each basin
        mask_p = basin_mask != 1
        mask_a = basin_mask != 2
        mask_i = basin_mask != 3
        var_diff_p = np.ma.array(var_diff, mask=mask_p)
        var_diff_a = np.ma.array(var_diff, mask=mask_a)
        var_diff_i = np.ma.array(var_diff, mask=mask_i)
        zonal_mean_p = np.ma.average(var_diff_p, axis=1)
        zonal_mean_a = np.ma.average(var_diff_a, axis=1)
        zonal_mean_i = np.ma.average(var_diff_i, axis=1)
        zonal_mean_p = np.tile(zonal_mean_p, (len(lon),1)); zonal_mean_p = np.transpose(zonal_mean_p)
        zonal_mean_a = np.tile(zonal_mean_a, (len(lon), 1)); zonal_mean_a = np.transpose(zonal_mean_a)
        zonal_mean_i = np.tile(zonal_mean_i, (len(lon), 1)); zonal_mean_i = np.transpose(zonal_mean_i)
        var_diff[basin_mask==1] = var_diff[basin_mask==1] - zonal_mean_p[basin_mask==1]
        var_diff[basin_mask==2] = var_diff[basin_mask==2] - zonal_mean_a[basin_mask==2]
        var_diff[basin_mask==3] = var_diff[basin_mask==3] - zonal_mean_i[basin_mask==3]


    # Create meshgrid
    lon2d, lat2d = np.meshgrid(lon, lat)

    # Levels for shade plot
    #levels = MaxNLocator(nbins=minmax[2]).tick_values(minmax[0], minmax[1])
    levels = np.linspace(minmax[0], minmax[1], minmax[2])

    # Read grid for coloring continents
    f2 = open_ncfile('/home/ysilvy/data/140807_WOD13_masks.nc', 'r')
    landsea = f2.variables['landsea'][:]

    # Create mask
    sea_mask = landsea != 1
    landsea = np.ma.array(landsea, mask=sea_mask)
    landsea[landsea == 1] = 0.2

    # Basemap
    map = Basemap(projection='cyl', llcrnrlon=20, llcrnrlat=-70, urcrnrlon=380, urcrnrlat=70, ax=ax1)
    map.drawmapboundary(fill_color='1')
    map.drawparallels(np.arange(-60, 61, 20.), labels=[1, 1, 0, 0], linewidth=0.5)
    map.drawmeridians(np.arange(-180, 180, 60), labels=[0, 0, 0, 1], linewidth=0.5)

    # Draw filled contours of diff
    cnplot = map.contourf(lon2d, lat2d, var_diff, cmap=cmap, levels = levels, latlon=True, extend='both')

    if kind == 'Durack':
        map.fillcontinents(color='black')
    else :
        # Draw continents
        pc2 = map.pcolormesh(lon2d, lat2d, landsea, shading='flat', cmap=plt.cm.gray, latlon=True)

    # Indicate which isopycnal we're projecting on
    ax1.text(0.1, 0.85, '$\gamma = %.1f$' %(isopyc,), transform=ax1.transAxes, fontweight='bold',color='w', fontsize=17)


    # Plot zonal mean
    ax2.plot(zonal_mean, lat, color='purple', linewidth=1.5)

    # Set axis parameters
    ax2.set_xlim(minmax[0], minmax[1])
    ax2.set_ylim(-70,70)
    ax2.tick_params(
        axis='both',  # changes apply to both axes
        which='both',  # both major and minor ticks are affected
        top='off',  # ticks along the top edge are off
        labelleft='off',
        #left='off',
        right='off',
        labelright='off')
    ax2.spines['left'].set_position('zero')
    ax2.spines['right'].set_color('none')
    ax2.spines['top'].set_color('none')
    ax2.set_yticks([-60,-40,-20,0,20,40,60])
    labels = ax2.get_xticklabels()
    plt.setp(labels, rotation=45, fontsize=10)
    ax2.grid(True)

    return cnplot, levels



# -----------------------------------------------
#   Class to define the midpoint of a colorbar
# -----------------------------------------------

class MidPointNorm(Normalize):
    def __init__(self, midpoint=0, vmin=None, vmax=None, clip=False):
        Normalize.__init__(self,vmin, vmax, clip)
        self.midpoint = midpoint

    def __call__(self, value, clip=None):
        if clip is None:
            clip = self.clip

        result, is_scalar = self.process_value(value)

        self.autoscale_None(result)
        vmin, vmax, midpoint = self.vmin, self.vmax, self.midpoint

        if not (vmin < midpoint < vmax):
            raise ValueError("midpoint must be between maxvalue and minvalue.")
        elif vmin == vmax:
            result.fill(0) # Or should it be all masked? Or 0.5?
        elif vmin > vmax:
            raise ValueError("maxvalue must be bigger than minvalue")
        else:
            vmin = float(vmin)
            vmax = float(vmax)
            if clip:
                mask = ma.getmask(result)
                result = ma.array(np.clip(result.filled(vmax), vmin, vmax),
                                  mask=mask)

            # ma division is very slow; we can take a shortcut
            resdat = result.data

            #First scale to -1 to 1 range, than to from 0 to 1.
            resdat -= midpoint
            resdat[resdat>0] /= abs(vmax - midpoint)
            resdat[resdat<0] /= abs(vmin - midpoint)

            resdat /= 2.
            resdat += 0.5
            result = ma.array(resdat, mask=result.mask, copy=False)

        if is_scalar:
            result = result[0]
        return result

    def inverse(self, value):
        if not self.scaled():
            raise ValueError("Not invertible until scaled")
        vmin, vmax, midpoint = self.vmin, self.vmax, self.midpoint

        if cbook.iterable(value):
            val = ma.asarray(value)
            val = 2 * (val-0.5)
            val[val>0]  *= abs(vmax - midpoint)
            val[val<0] *= abs(vmin - midpoint)
            val += midpoint
            return val
        else:
            val = 2 * (val - 0.5)
            if val < 0:
                return  val*abs(vmin-midpoint) + midpoint
            else:
                return  val*abs(vmax-midpoint) + midpoint


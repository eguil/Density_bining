#!/bin/env python
# -*- coding: utf-8 -*-


import numpy as np
import numpy.ma as ma
from mpl_toolkits.basemap import Basemap
#from matplotlib.ticker import MaxNLocator
from netCDF4 import Dataset as open_ncfile
from matplotlib.ticker import AutoMinorLocator
from scipy.interpolate import interp1d, InterpolatedUnivariateSpline

# --------------------------------
#       Variable properties
# --------------------------------

def defVar(longName):
    salinity = {
        'longN': 'salinity',  # long name
        'var': 'sog',  # variable name
        'minmax': [-0.3, 0.3, 16],  # for diff shading + number of color interval
        'clevsm': np.arange(30, 40, .25),  # for mean contours
        'clevsm_bold': np.arange(30, 40, 0.5),
        'clevsmdif': np.arange(-.2, .2, .025),  # for mean contours
        'clevsmstd': np.arange(0., .2, .005),  # for stddev contours
        '1dminmax': [-.1, .1], # for 1D ToE plots
        'legVar': "Salinity",  # Legend name
        'unit': "PSU",  # TODO: could be read from file
    }

    temp = {'var': 'thetaog', 'minmax': [-0.65, 0.65, 14], 'clevsm': np.arange(-2, 30, 1), 'clevsm_bold': np.arange(-2,30,2),
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
    salinity = {'var_zonal': 'isonsoBowl', 'var_zonal_w/bowl': 'isonso',
                'var_global': 'sogBowl', 'var_global_std':'sogBowlStd',
                'minmax': [-0.3, 0.3, 16],
                'minmax_zonal': [-0.2, 0.2, 16],
                'clevsm': np.arange(30, 40, .25),
                'clevsm_zonal': np.arange(30, 40, .1),
                'clevsm_bold': np.arange(30, 40, .5),
                'legVar': "Salinity", 'unit': "PSU", 'longN': 'salinity'}

    temp = {'var_zonal':'isonthetaoBowl', 'var_zonal_w/bowl': 'isonthetao',
            'var_global': 'thetaogBowl', 'var_global_std':'thetaogBowlStd',
            'minmax': [-0.65, 0.65, 14],
            'minmax_zonal' : [-0.4,0.4,16],
            'clevsm': np.arange(-2, 30, 1),
            'clevsm_zonal': np.arange(-2, 30, 1),
            'clevsm_bold': np.arange(-2,30,2),
            'legVar': "Temperature", 'unit': "C", 'longN': 'temp'}

    depth = {'var_zonal':'isondepthBowl',
             'clevsm_zonal': np.arange(0, 2000, 100),
             'clevsm_bold' : np.arange(0,2000,500),
             'minmax_zonal': [-50, 50, 16],
             'legVar': "Depth", 'unit': "m", 'longN': 'depth'}

    volume = {'var_zonal': 'isonvolBowl', 'var_zonal_w/bowl': 'isonvolBowl',
              'minmax_zonal': [-20., 20., 20],
              'clevsm_zonal': np.arange(0, 500, 50),
              'legVar': "Volume", 'unit': "1.e12 m^3", 'longN': 'volume'
              }

    vars = [salinity,temp,depth,volume]
    for ivar in range(len(vars)):
        if vars[ivar]['longN'] == longName:
            varout = vars[ivar]

    return varout


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

def zonal_2D(plt, action, ax0, ax1, ticks, lat, density, varBasin, domrho, cmap, levels, clevsm=None, clevsm_bold=None):

    # latitude domain
    domlat = [-70, 70]

    if action == 'total' :
        var = varBasin['var_change']
        var_mean = varBasin['var_mean']
        var_er = varBasin['var_error']

        # -- Error field for Durack & Wijffels data
        var_er = var_er * 1.1  # to account for a potential underestimation of the error determined by a bootstrap analysis
        var_er = var_er * 1.64  # 90% confidence level
        not_signif_change = np.where(np.absolute(var) < var_er, 1, 0)

        # -- Format for contour labels
        levfmt = '%.0f'
        if abs(clevsm[1] - clevsm[0]) < 1:
            levfmt = '%.1f'
        if abs(clevsm[1] - clevsm[0]) < 0.1:
            levfmt = '%.2f'

    elif action == 'total_mme':
        var = varBasin['var_change']
        var_mean = varBasin['var_mean']

        # -- Format for contour labels
        levfmt = '%.0f'
        if abs(clevsm[1] - clevsm[0]) < 1:
            levfmt = '%.1f'
        if abs(clevsm[1] - clevsm[0]) < 0.1:
            levfmt = '%.2f'

    elif action == 'isopyc_mig_sig':
        var = varBasin['dvar_dsig']*varBasin['dsig_dt']

    elif action == 'isopyc_mig_lat':
        var = varBasin['dvar_dy']*varBasin['dy_dt']

    elif action == 'mean_fields':
        var_1950 = varBasin['var_1950']
        var_2000 = varBasin['var_2000']

        # -- Format for contour labels
        levfmt = '%.0f'
        if abs(clevsm[1] - clevsm[0]) < 1:
            levfmt = '%.1f'
        if abs(clevsm[1] - clevsm[0]) < 0.1:
            levfmt = '%.2f'

    elif action == 'spiciness_fields':
        var_1950 = varBasin['dvar_dy_1950']
        var_2000 = varBasin['dvar_dy_2000']
        # -- Format for contour labels
        levfmt = '%.0f'
        if abs(levels[1] - levels[0]) < 1:
            levfmt = '%.1f'
        if abs(levels[1] - levels[0]) < 0.1:
            levfmt = '%.2f'

    else :
        var = varBasin[action]


    if action != 'total' and action != 'total_mme' and action != 'ToE':
        bowl = varBasin['bowl']
    if action == 'total_mme' or action == 'ToE' :
        bowl2 = varBasin['bowl2']
        bowl1 = varBasin['bowl1']
        label1 = varBasin['labBowl'][0]
        label2 = varBasin['labBowl'][1]

    # levels and color map
    levels = levels
    cmap = cmap

    # Create meshgrid
    lat2d, density2d = np.meshgrid(lat, density)

    # ==== Upper panel ====

    if action == 'mean_fields':
        ax0.contour(lat2d,density2d,var_1950, levels=levels, colors='black', linewidths=0.5)
        cpplot11 = ax0.contour(lat2d,density2d,var_1950, levels=clevsm_bold, colors='black', linewidths=1.5)
        ax0.clabel(cpplot11, inline=1, fontsize=11, fmt=levfmt)
        ax0.contour(lat2d,density2d,var_2000, levels=levels, colors='black', linewidths=0.5, linestyles='dashed')
        cpplot12 = ax0.contour(lat2d,density2d,var_2000, levels=clevsm_bold, colors='black', linewidths=1.5, linestyles='dashed')
        #ax0.clabel(cpplot12, inline=1, fontsize=11, fmt=levfmt)
    elif action == 'spiciness_fields':
        cnt_1950 = ax0.contour(lat2d,density2d,var_1950, levels=levels, colors='black')
        ax0.clabel(cnt_1950, inline=1, fontsize=11, fmt=levfmt)
        cnt_2000 = ax0.contour(lat2d,density2d,var_2000, levels=levels, colors='black', linestyles='dashed')
        ax0.clabel(cnt_2000, inline=1, fontsize=11, fmt=levfmt)
    else:
        cnplot1 = ax0.contourf(lat2d, density2d, var, levels=levels, cmap=cmap, extend='both')

    # # -- Draw mean contours
    # if (action == 'total_mme' and var_mean != None) or action == 'total' :
    #     cpplot11 = ax0.contour(lat2d, density2d, var_mean, clevsm, colors='black', linewidths=0.5)
    #     cpplot12 = ax0.contour(lat2d, density2d, var_mean, clevsm_bold, colors='black', linewidths=2)
    #     ax0.clabel(cpplot12, inline=1, fontsize=12, fontweight='bold', fmt=levfmt)

    # -- Draw areas where signal is not significant for D&W
    if action == 'total' :
        error_plot = ax0.contourf(lat2d, density2d, not_signif_change, levels=[0.25,0.5,1.5], colors='None',
                                   hatches=['','....'], edgecolor='0.3', linewidth=0.0)

    if action != 'var_2000_hr' and action != 'var_2000_sig_hr' and action!='total' and action!='total_mme' and action != 'ToE':
        ax0.plot(lat, bowl, color='black')

    if action == 'ToE' or action == 'total_mme':
        ax0.plot(lat, bowl2, linestyle = '--', linewidth=2, color='black',label=label2)
        ax0.plot(lat, bowl1, linewidth=2, color='black',label=label1)
        # -- Add legend for bowl position
        if varBasin['name'] == 'Indian':
            ax0.legend(loc='upper right', title='Bowl', fontsize=12)

    ax0.set_ylim([domrho[0], domrho[1]])
    ax0.set_xlim([domlat[0], domlat[1]])
    ax0.invert_yaxis()
    ax0.tick_params(
        axis='x',  # changes apply to the x axis
        which='both',  # both major and minor ticks are affected
        bottom='off',  # ticks along the bottom edge are off
        labelbottom='off',
        top='off')

    # # For selecting specific boxes, make a fine grid
    # xminorLocator = AutoMinorLocator(4)
    # yminorLocator = AutoMinorLocator(5)
    # ax0.xaxis.set_minor_locator(xminorLocator)
    # ax0.yaxis.set_minor_locator(yminorLocator)
    # ax0.grid(True, which='minor')
    # ax0.grid(True, which='major', ls='-')

    if ticks != 'left':
        ax0.tick_params(axis='y', labelleft='off')
    if ticks == 'right':
        ax0.tick_params(axis='y', labelright='on')

    ax0.axvline(x=0, color='black', ls='--')


    # === Lower panel ====

    if action == 'mean_fields':
        ax1.contour(lat2d,density2d,var_1950, levels=levels, colors='black', linewidths=0.5)
        cpplot21 = ax1.contour(lat2d,density2d,var_1950, levels=clevsm_bold, colors='black', linewidths=1.5)
        ax1.clabel(cpplot21, inline=1, fontsize=11, fmt=levfmt)
        ax1.contour(lat2d,density2d,var_2000, levels=levels, colors='black', linewidths=0.5, linestyles='dashed')
        cpplot22 = ax1.contour(lat2d,density2d,var_2000, levels=clevsm_bold, colors='black', linewidths=1.5, linestyles='dashed')
        #ax1.clabel(cpplot22, inline=1, fontsize=11, fmt=levfmt)
        cnplot2 = cpplot22
    elif action == 'spiciness_fields':
        cnt_1950 = ax1.contour(lat2d,density2d,var_1950, levels=levels, colors='black')
        ax1.clabel(cnt_1950, inline=1, fontsize=11, fmt=levfmt)
        cnt_2000 = ax1.contour(lat2d,density2d,var_2000, levels=levels, colors='black', linestyles='dashed')
        ax1.clabel(cnt_2000, inline=1, fontsize=11, fmt=levfmt)
        cnplot2 = cnt_2000
    else:
        cnplot2 = ax1.contourf(lat2d, density2d, var, levels=levels, cmap=cmap, extend='both')

    # # -- Draw mean contours
    # if (action == 'total_mme'  and var_mean != None) or action == 'total':
    #     cpplot21 = ax1.contour(lat2d, density2d, var_mean, clevsm, colors='black', linewidths=0.5)
    #     cpplot22 = ax1.contour(lat2d, density2d, var_mean, clevsm_bold, colors='black', linewidths=2)
    #     ax1.clabel(cpplot22, inline=1, fontsize=12, fontweight='bold', fmt=levfmt)

    # -- Draw areas where signal is not significant for D&W
    if action == 'total' :
        error_plot = ax1.contourf(lat2d, density2d, not_signif_change, levels=[0.25, 0.5, 1.5], colors='None',
                                  hatches=['', '....'], edgecolor='0.6', linewidth=0.0)

    if action != 'var_2000_hr' and action != 'var_2000_sig_hr' and action!='total' and action!='total_mme' and action != 'ToE':
        ax1.plot(lat, bowl, color='black')

    if action == 'ToE' or action == 'total_mme':
        ax1.plot(lat, bowl2, linestyle = '--', linewidth=2, color='black',label=label2)
        ax1.plot(lat, bowl1, linewidth=2, color='black',label=label1)


    ax1.set_ylim([domrho[1], domrho[2]])
    ax1.set_xlim([domlat[0], domlat[1]])
    ax1.invert_yaxis()
    ax1.tick_params(
        axis='x',  # changes apply to the x axis
        which='both',  # both major and minor ticks are affected
        top='off')  # ticks along the bottom edge are off

    # # For selecting specific boxes, make a fine grid
    # xminorLocator = AutoMinorLocator(4)
    # yminorLocator = AutoMinorLocator(5)
    # ax1.xaxis.set_minor_locator(xminorLocator)
    # ax1.yaxis.set_minor_locator(yminorLocator)
    # ax1.grid(True, which='minor')
    # ax1.grid(True, which='major', ls='-')

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


    return cnplot2



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
        # -- Plot areas where data is not significant
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
    f2 = open_ncfile('/home/ysilvy/Density_bining/Yona_analysis/data/140807_WOD13_masks.nc', 'r')
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
#          Average in lat/rho domain
# -----------------------------------------------

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


# -----------------------------------------------
#          Remap to Depth coordinates
# -----------------------------------------------

def remapToZ(fieldr,depthr,volumr, targetz, bowlz, v, bathy):
    '''
    The remaToZ() function remaps a density bined zonal field back to z
    It starts from the surface and computes the mean field for each z level, using the zonal volume of isopycnals for weighting

    Author:    Eric Guilyardi : Eric.Guilyardi@locean-ipsl.upmc.fr

    Created on Thu Aug 18 09:50:29 CEST 2016

    Inputs:
    ------
    - fieldr     - field (T or S)       - 3D basin,density,latitude
    - depthr     - depth of isopycnals  - 3D basin,density,latitude
    - volumr     - volume of isopycnals - 3D basin,density,latitude
    - targetz    - target z grid        - 1D
    - bowlz      - depth of bowl        - basin, latitude array

    Output:
    - fieldz    - 3D depth,latitude

    Usage:
    ------

    Notes:
    -----
    - EG 18 Aug 2016   - Initial function write
    - March 2017 Yona Silvy : updating script and adding an interpolation on each depth column
   '''

    latN = fieldr.shape[2]
    rhoN = fieldr.shape[1]
    basN = fieldr.shape[0]

    levN = len(targetz)

    # init fieldx array
    fieldz  = np.ma.masked_all([basN, levN, latN], dtype='float32')
    volumez = np.ma.masked_all([basN, levN, latN], dtype='float32')

    if v != 'V': # If not volume
        for ibasin in range(1,4):
            for j in range(latN):
                # Initialize local variables for interpolation to save levels that are not missing data
                iz_notempty = 0
                z_notempty = np.array([])
                fieldz_notempty = np.array([])
                for k in range(levN-1):
                    field_int = 0.
                    volum_int = 0.
                    for r in range(rhoN):
                        if volumr[ibasin,r,j] != 0:
                            if depthr[ibasin,r,j] >= targetz[k] and depthr[ibasin,r,j] < targetz[k+1]:
                                field_int = field_int + fieldr[ibasin,r,j]*volumr[ibasin,r,j]
                                volum_int = volum_int + volumr[ibasin,r,j]
                                depthr[ibasin,r,j] = -100. # to speed up search for next depths
                    if volum_int != 0.:
                        fieldz[ibasin,k,j] = field_int / volum_int
                        volumez[ibasin,k,j] = volum_int
                        # Save which levels are not missing data for extrapolating
                        z_notempty = np.append(z_notempty, targetz[k])
                        fieldz_notempty = np.append(fieldz_notempty, fieldz[ibasin,k,j])
                        iz_notempty = iz_notempty + 1
                    # Search bowl index for masking data later
                    if bowlz[ibasin,j] >= targetz[k] and bowlz[ibasin,j] < targetz[k+1] :
                        kbowl = k+1
                print('lat index', j)
                # print(iz_notempty-1, z_notempty.shape)
                if np.ma.is_masked(bowlz[ibasin,j]) == False :
                    # Interpolate the data on the depth column
                    if iz_notempty > 3:
                        # print 'Interpolate'
                        spl = InterpolatedUnivariateSpline(z_notempty, fieldz_notempty)
                        fieldz_new = spl(targetz)
                        fieldz[ibasin,:,j] = fieldz_new
                    # Mask field above the bowl
                    fieldz[ibasin,0:kbowl,j] = np.ma.masked
                # Mask bottom
                if np.ma.is_mask(bathy[ibasin,j]) == False and bathy[ibasin,j] < targetz[-1]:
                    bathy_mask = np.ma.nonzero(targetz>=bathy[ibasin,j])[0]
                    # print(targetz[bathy_mask[0]:])
                    fieldz[ibasin,bathy_mask[0]:,j] = np.ma.masked


    else :
        for ibasin in range(1,4):
            for j in range(latN):
                iz_notempty = 0
                z_notempty = np.array([])
                fieldz_notempty = np.array([])
                for k in range(levN-1):
                    field_int = 0.
                    for r in range(rhoN):
                        if depthr[ibasin,r,j] >= targetz[k] and depthr[ibasin,r,j] < targetz[k+1]:
                            field_int = field_int + fieldr[ibasin,r,j]
                            depthr[ibasin,r,j] = -100. # to speed up search for next depths
                    if field_int != 0.:
                        fieldz[ibasin,k,j] = field_int
                        # Save which levels are not missing data for extrapolating
                        z_notempty = np.append(z_notempty, targetz[k])
                        fieldz_notempty = np.append(fieldz_notempty,fieldz[ibasin,k,j])
                        iz_notempty = iz_notempty + 1
                    # Search bowl index for masking data later
                    if bowlz[ibasin,j] >= targetz[k] and bowlz[ibasin,j] < targetz[k+1] :
                        kbowl = k
                print('lat index', j)
                print(iz_notempty-1, z_notempty.shape)
                if np.ma.is_masked(bowlz[ibasin,j]) == False :
                    # Interpolate the data on the depth column
                    if iz_notempty > 3:
                        # print 'Interpolate'
                        spl = InterpolatedUnivariateSpline(z_notempty, fieldz_notempty)
                        fieldz_new = spl(targetz)
                        fieldz[ibasin,:,j] = fieldz_new
                    # Mask field below the bowl
                    fieldz[ibasin,0:kbowl,j] = np.ma.masked
                # Mask bottom
                if np.ma.is_mask(bathy[ibasin,j]) == False and bathy[ibasin,j] < targetz[-1]:
                    bathy_mask = np.ma.nonzero(targetz>=bathy[ibasin,j])[0]
                    print(targetz[bathy_mask[0]:])
                    fieldz[ibasin,bathy_mask[0]:,j] = np.ma.masked

    return fieldz



# ----------------------------------------------------
#   Build zonal latitude/depth plot
# ----------------------------------------------------

def zon_2Dz(plt, ax0, ax1, ticks, lat, lev, varBasin, clevsm, cmap, levels, domzed):

    # -- variables
    var = varBasin['var_change']
    bowl1 = varBasin['bowl1']
    bowl2 = varBasin['bowl2']

    # -- title and bowl labels
    label1 = varBasin['labBowl'][0]
    label2 = varBasin['labBowl'][1]

    # -- contour levels
    zedmin = domzed[0]
    zedmid = domzed[1]
    zedmax = domzed[2]

    # -- Latmin/max
    domlat = [-70, 70]

    # Create meshgrid
    lat2d, lev2d = np.meshgrid(lat, lev)

    #
    # ====  Upper panel  ===================================================
    #
    ax0.axis([domlat[0], domlat[1], zedmin, zedmid])
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

    # -- Format for contour labels
    levfmt='%.0f'
    if abs(clevsm[1]-clevsm[0]) < 1:
        levfmt='%.1f'
    if abs(clevsm[1]-clevsm[0]) < 0.1:
        levfmt='%.2f'

    # -- draw filled contours of period diff
    cnplot = ax0.contourf(lat2d, lev2d, var, cmap=cmap, levels=levels, extend='both')

    # # -- draw mean contours --> TODO compute mean field first
    # cpplot = ax0.contour(lat2d, lev2d, var, clevsm, colors='black', linewidths=0.5)
    # ax0.clabel(cpplot, inline=1, fontsize=10, fmt=levfmt)

    # -- draw bowl
    ax0.plot(lat, bowl1, color='black', linewidth=2, label=label1)
    ax0.plot(lat, bowl2, color='black', linewidth=2, linestyle='--', label=label2)

    if varBasin['name'] == 'Indian':
        ax0.legend(loc='upper right', title='Bowl', fontsize=12)

    #
    # ====  Lower panel   ===================================================
    #
    ax1.axis([domlat[0], domlat[1], zedmid, zedmax])
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

    # -- Re-label x-axis
    xlabels = ['', '60S', '40S', '20S', '0', '20N', '40N', '60N']
    ax1.set_xticklabels(xlabels)
    # -- Set y ticks
    ax1.set_yticks([500,1000,1500,2000])
    yminorLocator = AutoMinorLocator(4)
    ax1.yaxis.set_minor_locator(yminorLocator)

    # -- draw filled contours
    cnplot = ax1.contourf(lat2d, lev2d, var, cmap=cmap, levels=levels, extend='both')

    # # -- draw mean contours
    # cpplot = ax1.contour(lat2d, lev2d, var, clevsm, colors='black', linewidths=0.5)
    # ax1.clabel(cpplot, inline=1, fontsize=10, fmt=levfmt)

    # -- draw bowl
    ax1.plot(lat, bowl1, color='black', linewidth=2, label=label1)
    ax1.plot(lat, bowl2, color='black', linewidth=2, linestyle='--', label=label2)


    # Remove intersecting tick at zedmid
    yticks = ax1.yaxis.get_major_ticks()
    if ticks == 'left':
        yticks[0].label1.set_visible(False)
    if ticks == 'right':
        yticks[0].label2.set_visible(False)


    # -- add plot title
    ax1.text(domlat[0] + 10, zedmax-100, varBasin['name'], fontsize=16, fontweight='bold')

    return cnplot
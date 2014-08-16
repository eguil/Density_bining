#!/usr/local/uvcdat/latest/bin/cdat
##!/Users/ericg/Projets/CMIP/Metrics/WGNE/bin/python
#
# 
# Program to compute density bins and replace vertical z coordinate by neutral density
# Reads in netCDF T(x,y,z,t) and S(x,y,z,t) files and writes 
#  - T or S(x,y,sigma,t)
#  - D(x,y,sigma,t) (depth of isopycnal)
#  - Z(x,y,sigma,t) (thickness of isopycnal)
#
#
# Uses McDougall and Jackett 2005 EOS (IDL routine provided by G. Madec)
# Inspired from IDL density bin routines by G. Roullet and G. Madec (1999)
#
# --------------------------------------------------------------
#  E. Guilyardi - while at LBNL/LLNL -  March 2014
# 
# git add densit_bin.py
# git commit -m "with Paul"
# git push
# git checkout master
# git checkout <branch>

import cdms2 as cdm
import MV2 as mv
import os, sys
#import socket
import argparse
import string
import numpy as npy
import numpy.ma as ma
import cdutil as cdu
from genutil import statistics
import support_density as sd
import time as timc
import timeit
import resource
import ESMP
from cdms2 import CdmsRegrid
#import ZonalMeans
#from regrid import Regridder
#import matplotlib.pyplot as plt
#
# Keep track of time (CPU and elapsed)
ti0 = timc.clock()
te0 = timeit.default_timer()
#
# netCDF compression (use 0 for netCDF3)
comp = 0
cdm.setNetcdfShuffleFlag(comp)
cdm.setNetcdfDeflateFlag(comp)
cdm.setNetcdfDeflateLevelFlag(comp)
cdm.setAutoBounds('on')
#
# == Arguments
#
# 
# == Inits
#
npy.set_printoptions(precision = 2)
home   = os.getcwd()
outdir = os.getcwd()
hist_file_dir=home
#
# == Arguments
#
# == get command line options
parser = argparse.ArgumentParser(description = 'Script to perform density bining analysis')
parser.add_argument('-d', help = 'toggle debug mode', action = 'count', default = 0)
#parser.add_argument('-r','--sigma_range', help='neutral sigma range', required=True)
#parser.add_argument('-s','--sigma_increment', help='neutral sigma increment', required=True)
#parser.add_argument('-i','--input', help='input directory', default="./")
#parser.add_argument('-o','--output',help='output directory', default="./")
parser.add_argument('-t','--timeint', help='specify time domain in bining <init_idx>,<ncount>', default="all")
parser.add_argument('-n','--nomthoutput', help = 'no monthly output', action = 'count', default = 0)
#parser.add_argument('string', metavar='root for T and S files', type=str, help='netCDF input files root')
args = parser.parse_args()
#
# Write command line in history file
filer=hist_file_dir+'/z_density_hist.txt'
with open(filer, 'a') as f:
    f.write('\n\r'+str(sys.argv).translate(None, "',[]"))
# 
# read values
debug        = str(args.d)
#indir        = args.input
#outdir       = args.output
#sigma_range  = args.sigma_range 
#delta_sigma  = args.sigma_increment
timeint      = args.timeint
mthout       = args.nomthoutput
#file_root    = args.string
#
# Define T and S file names (local mac...)
#file_T = indir+'/'+file_root+'_thetao.nc'
#file_S = indir+'/'+file_root+'_so.nc'
#file_fx = indir+'/areacello_fx_IPSL-CM5A-LR_piControl_r0i0p0.nc'
#
#if socket.gethostname() == 'crunchy.llnl.gov':
#    file_T = indir+'/thetao/cmip5.'+file_root+'.thetao.ver-v20111010.latestX.xml'
#    file_S = indir+'/so/cmip5.'+file_root+'.so.ver-v20111010.latestX.xml'

#
# IPSL-CM5A-LR
#
file_fx = '/work/cmip5/fx/fx/areacello/cmip5.IPSL-CM5A-LR.piControl.r0i0p0.fx.ocn.fx.areacello.ver-v20120430.latestX.xml'
file_T = '/work/cmip5/historical/ocn/mo/thetao/cmip5.IPSL-CM5A-LR.historical.r1i1p1.mo.ocn.Omon.thetao.ver-v20111119.latestX.xml'
file_S = '/work/cmip5/historical/ocn/mo/so/cmip5.IPSL-CM5A-LR.historical.r1i1p1.mo.ocn.Omon.so.ver-v20111119.latestX.xml'
modeln = 'IPSL-CM5A-LR'
#
# GFDL-CM2p1
#
#file_fx = '/work/cmip5/fx/fx/areacello/cmip5.GFDL-CM2p1.historical.r0i0p0.fx.ocn.fx.areacello.ver-v20110601.latestX.xml'
#file_T  = '/work/cmip5/historical/ocn/mo/thetao/cmip5.GFDL-CM2p1.historical.r1i1p1.mo.ocn.Omon.thetao.ver-v20110601.latestX.xml'
#file_S  = '/work/cmip5/historical/ocn/mo/so/cmip5.GFDL-CM2p1.historical.r1i1p1.mo.ocn.Omon.so.ver-v20110601.latestX.xml'
#modeln = 'GFDL-CM2p1'
#
# CCSM4
#
#file_fx = '/work/cmip5/fx/fx/areacello/cmip5.CCSM4.historical.r0i0p0.fx.ocn.fx.areacello.ver-v20130312.latestX.xml'
#file_T = '/work/cmip5/historical/ocn/mo/thetao/cmip5.CCSM4.historical.r1i1p1.mo.ocn.Omon.thetao.ver-v20121128.latestX.xml'
#file_S = '/work/cmip5/historical/ocn/mo/so/cmip5.CCSM4.historical.r1i1p1.mo.ocn.Omon.so.ver-v20121128.latestX.xml'
#modeln = 'CCSM4'

#file_fx = '/Users/ericg/Desktop/Data/CMIP5/piControl/test_3d_ocn/areacello_fx_IPSL-CM5A-LR_piControl_r0i0p0.nc'
#file_T = '/Users/ericg/Desktop/Data/CMIP5/piControl/test_3d_ocn/IPSL-CM5A-LR_piControl_r1i1p1_180001-180012_Omon_thetao.nc'
#file_S = '/Users/ericg/Desktop/Data/CMIP5/piControl/test_3d_ocn/IPSL-CM5A-LR_piControl_r1i1p1_180001-180012_Omon_so.nc'
#
if debug >= '1':
    print ' Debug - File names:'
    print '    ', file_T
    print '    ', file_S
    debugp = True
else:
    debugp = False
#
# Open files
ft  = cdm.open(file_T)
fs  = cdm.open(file_S)
timeax = ft.getAxis('time')
#
# Dates to read
if timeint == 'all':
    tmin = 0
    tmax = timeax.shape[0]
else:
    tmin = int(timeint.split(',')[0]) - 1
    tmax = tmin + int(timeint.split(',')[1])

if debugp:
    #print; print ' Debug - Read only first month...'
    print; print ' Debug mode'
    #tmin = 0
    #tmax = 1
#
print '  time interval: ', tmin, tmax - 1
#
# Define temperature and salinity arrays
temp = ft('thetao', time = slice(0,1))
so   = fs('so', time = slice(0,1))
#
# Read file attributes
list_file=ft.attributes.keys()
file_dic={}
for i in range(0,len(list_file)):
    file_dic[i]=list_file[i],ft.attributes[list_file[i] ]
#
# Read masking value
valmask = so._FillValue
#
# Read time and grid
lon  = temp.getLongitude()
lat  = temp.getLatitude()
depth = temp.getLevel()
bounds = ft('lev_bnds')
ingrid = temp.getGrid()
#
# Read cell area
ff = cdm.open(file_fx)
area = ff('areacello')
ff.close()
#
# Define dimensions
N_i = int(lon.shape[1])
N_j = int(lon.shape[0])
N_z = int(depth.shape[0])
#N_t = int(time.shape[0])
#
# Define sigma grid
rho_min = 19
rho_max = 28
del_s = 0.2
s_s = npy.arange(rho_min, rho_max, del_s)
N_s = len(s_s)
sigma_bnds = mv.asarray([[s_s[:]],[s_s[:]+del_s]]) # make bounds for zonal mean computation
s_s = npy.tile(s_s, N_i*N_j).reshape(N_i*N_j,N_s).transpose() # make 3D for matrix computation
#
# Define zonal grid
delta_lat = 1.

#w=sys.stdin.readline() # stop the code here. [Ret] to keep going
#
# start density bining (loop on t and horizontal grid)
imin = 0
jmin = 0
imax = temp.shape[3]
jmax = temp.shape[2]
#
# test point  
itest = 80 
jtest = 60
ijtest = jtest*N_i+itest
#imin = itest
#jmin = jtest
#imax = itest+1
#jmax = jtest+1
#
# Define time read interval (as function of 3D array size)
grdsize = N_i * N_j * N_z
# define number of months in each chunk
print 'Grid size:', grdsize
if grdsize <= 1.e6:
    tcdel = min(120, tmax)
elif grdsize <= 1.e7:
    tcdel = min(24, tmax)
nyrtc = tcdel/12
tcmax = (tmax-tmin)/tcdel ; # number of time chunks
print ' ==> tcdel, tcmax:', tcdel, tcmax
#
# inits
# z profiles:
z_zt = depth[:]
z_zw = bounds.data[:,0]
#bowl_s = 0.
#
print
toc = timc.clock()
toc2 = timeit.default_timer()
#
# File output inits
#
s_sd = npy.arange(rho_min, rho_max+del_s, del_s, dtype = npy.float32)
s_axis = cdm.createAxis(s_sd, id = 'rhon')
s_axis.long_name = 'Neutral density'
s_axis.units = ''
s_axis.designateLevel()
grd = temp.getGrid()
# Monthly mean of T,S, thickness and depth on neutral density bins on source grid
file_out = outdir+'/'+modeln+'_out_density.nc'
if os.path.exists(file_out):
    os.remove(file_out)
if mthout == 0:
    g = cdm.open(file_out,'w+')
# Annual zonal mean of T,S, thick, depth and volume per basin on WOA grid
filez_out = outdir+'/'+modeln+'_outz_density.nc'
if os.path.exists(filez_out):
    os.remove(filez_out)
gz = cdm.open(filez_out,'w+')
# Annual mean zonal mean of persistence on WOA grid + volume of persistent domain
filep_out = outdir+'/'+modeln+'_out_persist.nc'
if os.path.exists(filep_out):
    os.remove(filep_out)
gp = cdm.open(filep_out,'w+')
#
# output arrays for each chunk
depth_bin = npy.ma.ones([tcdel, N_s+1, N_j*N_i], dtype='float32')*valmask 
depth_bin = mv.masked_where(mv.equal(depth_bin,valmask), depth_bin)
thick_bin = npy.ma.ones([tcdel, N_s+1, N_j*N_i], dtype='float32')*valmask 
thick_bin = mv.masked_where(mv.equal(thick_bin,valmask), thick_bin)
x1_bin = npy.ma.ones([tcdel, N_s+1, N_j*N_i], dtype='float32')*valmask 
x1_bin = mv.masked_where(mv.equal(x1_bin,valmask), x1_bin)
x2_bin = npy.ma.ones([tcdel, N_s+1, N_j*N_i], dtype='float32')*valmask 
x2_bin = mv.masked_where(mv.equal(x2_bin,valmask), x2_bin)

# target horizonal grid for interp 
fileg = '/work/guilyardi/Density_bining/WOD13_masks.nc'
#fileg = '/Users/ericg/Projets/Density_bining/WOD13_masks.nc'
gt = cdm.open(fileg)
maskg = gt('basinmask3')
outgrid = maskg.getGrid()
# global mask
maski = maskg.mask
# regional masks
maskAtl = maski*1 ; maskAtl[...] = True
idxa = npy.argwhere(maskg == 1).transpose()
maskAtl[idxa[0],idxa[1]] = False
maskPac = maski*1 ; maskPac[...] = True
idxp = npy.argwhere(maskg == 2).transpose()
maskPac[idxp[0],idxp[1]] = False
maskInd = maski*1 ; maskInd[...] = True
idxi = npy.argwhere(maskg == 3).transpose()
maskInd[idxi[0],idxi[1]] = False
#
loni = maskg.getLongitude()
lati = maskg.getLatitude()
Nii = int(loni.shape[0])
Nji = int(lati.shape[0])
# Compute area of target grid and zonal sums
areai = sd.compute_area(loni[:], lati[:])
#areai   = gt('basinmask3_area').data*1.e6
gt.close()
#
areazt  = cdu.averager(areai*maski  , axis=1, action='sum')
areazta = cdu.averager(areai*maskAtl, axis=1, action='sum')
areaztp = cdu.averager(areai*maskPac, axis=1, action='sum')
areazti = cdu.averager(areai*maskInd, axis=1, action='sum')
#
# Interpolation init
ESMP.ESMP_Initialize()
regridObj = CdmsRegrid(ingrid, outgrid, depth_bin.dtype, missing = valmask, regridMethod = 'linear', regridTool = 'esmf')
#
# Global arrays init
depthBini = npy.ma.ones([nyrtc, N_s+1, Nji, Nii], dtype='float32')*valmask 
thickBini = npy.ma.ones([nyrtc, N_s+1, Nji, Nii], dtype='float32')*valmask 
x1Bini    = npy.ma.ones([nyrtc, N_s+1, Nji, Nii], dtype='float32')*valmask 
x2Bini    = npy.ma.ones([nyrtc, N_s+1, Nji, Nii], dtype='float32')*valmask 
# Atl
# TODO: this is a lot of arrays - maybe there is a better way of doing this
depthBinia = npy.ma.ones([nyrtc, N_s+1, Nji, Nii], dtype='float32')*valmask 
thickBinia = npy.ma.ones([nyrtc, N_s+1, Nji, Nii], dtype='float32')*valmask 
x1Binia    = npy.ma.ones([nyrtc, N_s+1, Nji, Nii], dtype='float32')*valmask 
x2Binia    = npy.ma.ones([nyrtc, N_s+1, Nji, Nii], dtype='float32')*valmask 
# Pac
depthBinip = npy.ma.ones([nyrtc, N_s+1, Nji, Nii], dtype='float32')*valmask 
thickBinip = npy.ma.ones([nyrtc, N_s+1, Nji, Nii], dtype='float32')*valmask 
x1Binip    = npy.ma.ones([nyrtc, N_s+1, Nji, Nii], dtype='float32')*valmask 
x2Binip    = npy.ma.ones([nyrtc, N_s+1, Nji, Nii], dtype='float32')*valmask 
# Ind
depthBinii = npy.ma.ones([nyrtc, N_s+1, Nji, Nii], dtype='float32')*valmask 
thickBinii = npy.ma.ones([nyrtc, N_s+1, Nji, Nii], dtype='float32')*valmask 
x1Binii    = npy.ma.ones([nyrtc, N_s+1, Nji, Nii], dtype='float32')*valmask 
x2Binii    = npy.ma.ones([nyrtc, N_s+1, Nji, Nii], dtype='float32')*valmask 
#
# Persistence arrays
persist   = npy.ma.ones([nyrtc, N_s+1, N_j, N_i], dtype='float32')*valmask
persisti  = npy.ma.ones([nyrtc, N_s+1, Nji, Nii], dtype='float32')*valmask 
#
# loop on time chunks
for tc in range(tcmax):
    tuc = timc.clock()
    # read tcdel month by tcdel month to optimise memory
    trmin = tmin + tc*tcdel ; # define as function of tc and tcdel
    trmax = tmin + (tc+1)*tcdel ; # define as function of tc and tcdel
    print ' --> time chunk (bounds) = ',tc, ' (',trmin,trmax-1,')'
    temp = ft('thetao', time = slice(trmin,trmax))
    so   = fs('so'    , time = slice(trmin,trmax))
    time  = temp.getTime()
    # Kelvin or celsius ?
    tempmin = min(temp.data[0,:,N_j/2,N_i/2])
    print ' tempmin = ',tempmin
    if tempmin > valmask/10.:
        tempmin = min(temp.data[0,:,N_j/4,N_i/2])
        print ' tempmin = ',tempmin
    if tempmin > 273.:
        temp = temp -273.15
        print ' Change unit to celsius'
    #
    # Compute neutral density 
    rhon = sd.eos_neutral(temp,so)-1000.
    #
    # reorganise i,j dims in single dimension data
    temp = npy.reshape(temp, (tcdel, N_z, N_i*N_j))
    so   = npy.reshape(so  , (tcdel, N_z, N_i*N_j))
    rhon = npy.reshape(rhon, (tcdel, N_z, N_i*N_j))
    #
    # init output arrays
    depth_bin[...] = valmask
    thick_bin[...] = valmask
    x1_bin[...] = valmask
    x2_bin[...] = valmask
    #
    # Loop on time within chunk tc
    for t in range(trmax-trmin): 
        tac0 = timc.clock()
        if (t/20*20 == t): 
            print '      t = ',t
        # x1 contents on vertical (not yet implemented - may be done to ensure conservation)
        x1_content = temp.data[t] 
        x2_content = so.data[t] 
        vmask_3D = mv.masked_values(temp.data[t], valmask).mask 
        # find non-masked points
        nomask = npy.equal(vmask_3D[0],0)
        # init arrays
        z_s   = npy.ones((N_s+1, N_i*N_j))*valmask
        c1_s  = npy.ones((N_s+1, N_i*N_j))*valmask
        c2_s  = npy.ones((N_s+1, N_i*N_j))*valmask
        szmin  = npy.ones((N_i*N_j))*valmask
        szmax  = npy.ones((N_i*N_j))*valmask
        i_min = npy.ones((N_i*N_j))*0
        i_max = npy.ones((N_i*N_j))*0
        delta_rho = npy.ones((N_i*N_j))*valmask
        # find bottom level
        i_bottom = vmask_3D.argmax(axis=0)-1
        z_s [N_s, nomask] = z_zw[i_bottom[nomask]+1] ; # Cell depth limit
        c1_s[N_s, nomask] = x1_content[N_z-1,nomask] ; # Cell bottom temperature/salinity
        c2_s[N_s, nomask] = x2_content[N_z-1,nomask] ; # Cell bottom tempi_profilerature/salinity
        # init arrays = f(z)
        s_z = rhon.data[t]
        c1_z = x1_content
        c2_z = x2_content
        #
        # Extract a strictly increasing sub-profile
        i_min[nomask] = s_z.argmin(axis=0)[nomask]
        i_max[nomask] = s_z.argmax(axis=0)[nomask]-1
        i_min[i_min > i_max] = i_max[i_min > i_max]
        # Test on bottom - surface stratification
        delta_rho[nomask] = s_z[i_bottom[nomask],nomask] - s_z[0,nomask]
        i_min[delta_rho < del_s] = 0
        i_max[delta_rho < del_s] = i_bottom[delta_rho < del_s]
        #
        # General case
        # find min/max of density for each z profile
        # 
        tac00 = timc.clock()
        for i in range(N_i*N_j):
            if nomask[i]:
                szmin[i] = s_z[i_min[i],i]
                szmax[i] = s_z[i_max[i],i]
            else:
                szmin[i] = 0.
                szmax[i] = rho_max+10.
        tac3 = timc.clock()
        # Find indices between min and max of density (costing ~ 12% CPU)
        #
        # Construct arrays of szm/c1m/c2m = s_z[i_min[i]:i_max[i],i] and 'NaN' otherwise
        # same for zztm from z_zt (30% CPU)
        szm = s_z*1. ; szm[...] = 'NaN'
        zzm = s_z*1. ; zzm[...] = 'NaN'
        c1m = c1_z*1. ; c1m[...] = 'NaN'
        c2m = c2_z*1. ; c2m[...] = 'NaN'
        
        tac4 = timc.clock()        
        for k in range(N_z):
            k_ind = i_min*1.; k_ind[:] = valmask
            k_ind = npy.argwhere( (k >= i_min) & (k <= i_max))
            szm[k,k_ind] = s_z [k,k_ind]
            c1m[k,k_ind] = c1_z[k,k_ind]
            c2m[k,k_ind] = c2_z[k,k_ind]
            zzm[k,k_ind] = z_zt[k]
            
        tac5 = timc.clock()
        # interpolate depth(z) (z_zt) to depth(s) at s_s densities (z_s) using density(z) s_z
        # TODO: no loop (35% CPU)
        for i in range(N_i*N_j):
            if nomask[i]:
                z_s [0:N_s,i] = npy.interp(s_s[:,i], szm[:,i], zzm[:,i]) ; # consider spline           
                c1_s[0:N_s,i] = npy.interp(z_s[0:N_s,i], zzm[:,i], c1m[:,i]) 
                c2_s[0:N_s,i] = npy.interp(z_s[0:N_s,i], zzm[:,i], c2m[:,i]) 
        tac6 = timc.clock()
        if debugp and (t==0):
                print '    CPU stage 0 :', tac00-tac0
                print '    CPU stage 1 :', tac3-tac00
                print '    CPU stage 2 :', tac4-tac3
                print '    CPU stage 3 :', tac5-tac4
                print '    CPU stage 4 :', tac6-tac5
        #
        # if level in s_s has lower density than surface, isopycnal is put at surface (z_s=0)
        inds = npy.argwhere(s_s < szmin).transpose()
        z_s [inds[0],inds[1]] = 0.
        c1_s[inds[0],inds[1]] = valmask
        c2_s[inds[0],inds[1]] = valmask
        # if level of s_s has higher density than bottom density, 
        # isopycnal is set to bottom (z_s=z_zw[i_bottom])
        inds = npy.argwhere(s_s > szmax).transpose()
        z_s [inds[0],inds[1]] = z_s[N_s,inds[1]]
        c1_s[inds[0],inds[1]] = valmask
        c2_s[inds[0],inds[1]] = valmask
 
 
        depth_bin [t,:,:]     = z_s
        thick_bin [t,0,:]     = z_s[0,:]
        thick_bin [t,1:N_s,:] = z_s[1:N_s,:]-z_s[0:N_s-1,:]
        x1_bin    [t,:,:]     = c1_s
        x2_bin    [t,:,:]     = c2_s
    #
    # end of loop on t <===      
    #        
    # Reshape i*j back to i,j
    depth_bino = npy.reshape(depth_bin, (tcdel, N_s+1, N_j, N_i))
    thick_bino = npy.reshape(thick_bin, (tcdel, N_s+1, N_j, N_i))
    x1_bino    = npy.reshape(x1_bin,    (tcdel, N_s+1, N_j, N_i))
    x2_bino    = npy.reshape(x2_bin,    (tcdel, N_s+1, N_j, N_i))
    #
    # Wash mask over variables
    maskb = mv.masked_values(x1_bino, valmask).mask
    depth_bino.mask = maskb
    thick_bino.mask = maskb
    x1_bino.mask = maskb
    x2_bino.mask = maskb
    #
    tucf = timc.clock()
    #
    if debugp and (tc == 0):
        # test write
        i = itest
        j = jtest
        #print 'ind = ',ind
        print "test point",i,j, area[j,i]
        print "lon,lat",lon[j,i],lat[j,i]
        print 'depth_bin', depth_bino[0,:,j,i]
        print 'thick_bin', thick_bino[0,:,j,i]
        print 'x1_bin', x1_bino[0,:,j,i]
        print 'x2_bin', x2_bino[0,:,j,i]
    tic = timc.clock()
    #
    # Output files as netCDF
    # Def variables 
    # QQ??: only do for tc==0 ? depth_bin update enought for tc >= 1 ?
    depthBin = cdm.createVariable(depth_bino, axes = [time, s_axis, grd], id = 'isondepth')
    thickBin = cdm.createVariable(thick_bino, axes = [time, s_axis, grd], id = 'isonthick')
    x1Bin    = cdm.createVariable(x1_bino   , axes = [time, s_axis, grd], id = 'thetao')
    x2Bin    = cdm.createVariable(x2_bino   , axes = [time, s_axis, grd], id = 'so')
    if mthout == 0:
        if tc == 0:
            depthBin.long_name = 'Depth of isopycnal'
            depthBin.units = 'm'
            #
            thickBin.long_name = 'Thickness of isopycnal'
            thickBin.units = 'm'
            x1Bin.long_name = temp.long_name
            x1Bin.units = 'C'
            x2Bin.long_name = so.long_name
            x2Bin.units = so.units
            #
            g.write(area) ; # Added area so isonvol can be computed
            # write global attributes (inherited from thetao file)
            for i in range(0,len(file_dic)):
                dm=file_dic[i]
                setattr(g,dm[0],dm[1])
                setattr(g,'Post_processing_history','Density bining via densit_bin.py using delta_sigma = '+str(del_s))
                setattr(gz,'Post_processing_history','Zonal mean annual Density bining via densit_bin.py using monthly means and delta_sigma = '+str(del_s))
    #
    # Compute annual mean, persistence, make zonal mean and write
    # 
    ticz = timc.clock()
    if tcdel >= 12:
        # Annual mean
        # Note: HUGE COST: 30-60 sec for 12 months !!! and 120 sec for 24 !!!
        dy  = cdu.YEAR(depthBin)
        ty  = cdu.YEAR(thickBin)
        x1y = cdu.YEAR(x1Bin)
        x2y = cdu.YEAR(x2Bin)
        # this is 5 times cheaper but no grid is passed and ZonalMeans fails
        #dy  = cdu.averager(npy.reshape (depthBin, (nyrtc, 12, N_s+1, N_j, N_i)), axis=1)
        #ty  = cdu.averager(npy.reshape (thickBin, (nyrtc, 12, N_s+1, N_j, N_i)), axis=1)
        #x1y = cdu.averager(npy.reshape (x1Bin,    (nyrtc, 12, N_s+1, N_j, N_i)), axis=1)
        #xy2 = cdu.averager(npy.reshape (x2Bin,    (nyrtc, 12, N_s+1, N_j, N_i)), axis=1)
        toz = timc.clock()
        #if debugp:
        #
        # Compute annual persistence of isopycnal bins (from their thickness)
        #  = percentage of time bin is occupied during each year (annual bowl if % < 100)
        for t in range(nyrtc):
            # inits
            idxvm = npy.ma.ones([12, N_s+1, N_j, N_i], dtype='float32')*valmask 
            inim = (nyrtc-1)*12
            finm = (nyrtc-1)*12 + 11
            idxvm = 1-mv.masked_values(thick_bino[inim:finm,:,:,:], valmask).mask 
            persist[t,:,:,:] = cdu.averager(idxvm, axis=0)*100.
            # mask where value is zero
            persist._FillValue = valmask
            persist = mv.masked_where(persist <= 1.e-6, persist)
            persist[mv.masked_values(persist, valmask).mask]=valmask
            persbin = cdm.createVariable(persist, axes = [dy.getAxis(0), s_axis, grd], id = 'isonpers')           
            # regrid
            for ks in range(N_s+1):
                persisti[t,ks,:,:] = persbin [t,ks,:,:].regrid(outgrid, regridTool='ESMF', regridMethod='linear')
                #persisti[t,ks,:,:] = regridObj(persbin [t,ks,:,:])
                persisti[t,ks,:,:].mask = maski
            # Compute zonal mean
            # Global
            persistiz = cdu.averager(persisti, axis=3)
            # TO DO:
            #     - make zonal mean per basins (2D)
            #     - compute volume/temp/salinity of persistent ocean (global, per basin) (1D)
            #     - write output in file 
        #
        # Write persistence variables
        dbpz  = cdm.createVariable(persistiz, axes = [dy.getAxis(0), s_axis, lati], id = 'isonpers')
        if tc == 0:
            # Global attributes
            persbin.long_name = 'persistence of isopycnal bins'
            persbin.units = '% of time'
            dbpz.long_name = 'zonal persistence of isopycnal bins'
            dbpz.units = '% of time'
        # Write & append
        #gp.write(persbin , extend = 1, index = (trmin-tmin)/12)
        gp.write(dbpz    , extend = 1, index = (trmin-tmin)/12)
        #
        tozp = timc.clock()
            
        # Interpolate onto common grid (~90% of total CPU !!)
        for t in range(nyrtc):
            for ks in range(N_s+1):
                # Global
                #depthBini[t,ks,:,:] = dy [t,ks,:,:].regrid(outgrid, regridTool='ESMF', regridMethod='linear', diag = diag)
                #thickBini[t,ks,:,:] = ty [t,ks,:,:].regrid(outgrid, regridTool='ESMF', regridMethod='linear', diag = diag)
                #x1Bini   [t,ks,:,:] = x1y[t,ks,:,:].regrid(outgrid, regridTool='ESMF', regridMethod='linear', diag = diag)
                #x2Bini   [t,ks,:,:] = x2y[t,ks,:,:].regrid(outgrid, regridTool='ESMF', regridMethod='linear', diag = diag)
                depthBini[t,ks,:,:] = regridObj(dy [t,ks,:,:])
                thickBini[t,ks,:,:] = regridObj(ty [t,ks,:,:])
                x1Bini   [t,ks,:,:] = regridObj(x1y[t,ks,:,:])
                x2Bini   [t,ks,:,:] = regridObj(x2y[t,ks,:,:])
                #
                depthBini[t,ks,:,:].mask = maski
                thickBini[t,ks,:,:].mask = maski
                x1Bini   [t,ks,:,:].mask = maski
                x2Bini   [t,ks,:,:].mask = maski
                # Atl
                # TODO: many arrays - there is maybe a way to optimize this
                depthBinia[t,ks,:,:] = depthBini[t,ks,:,:]*1.
                thickBinia[t,ks,:,:] = thickBini[t,ks,:,:]*1.
                x1Binia[t,ks,:,:]    = x1Bini[t,ks,:,:]*1.
                x2Binia[t,ks,:,:]    = x2Bini[t,ks,:,:]*1.
                depthBinia[t,ks,:,:].mask = maskAtl
                thickBinia[t,ks,:,:].mask = maskAtl
                x1Binia   [t,ks,:,:].mask = maskAtl
                x2Binia   [t,ks,:,:].mask = maskAtl
                # Pac
                depthBinip[t,ks,:,:] = depthBini[t,ks,:,:]*1.
                thickBinip[t,ks,:,:] = thickBini[t,ks,:,:]*1.
                x1Binip   [t,ks,:,:] = x1Bini[t,ks,:,:]*1.
                x2Binip   [t,ks,:,:] = x2Bini[t,ks,:,:]*1.
                depthBinip[t,ks,:,:].mask = maskPac
                thickBinip[t,ks,:,:].mask = maskPac
                x1Binip   [t,ks,:,:].mask = maskPac
                x2Binip   [t,ks,:,:].mask = maskPac
                # Ind
                depthBinii[t,ks,:,:] = depthBini[t,ks,:,:]*1.
                thickBinii[t,ks,:,:] = thickBini[t,ks,:,:]*1.
                x1Binii   [t,ks,:,:] = x1Bini[t,ks,:,:]*1.
                x2Binii   [t,ks,:,:] = x2Bini[t,ks,:,:]*1.
                depthBinii[t,ks,:,:].mask = maskInd
                thickBinii[t,ks,:,:].mask = maskInd
                x1Binii   [t,ks,:,:].mask = maskInd
                x2Binii   [t,ks,:,:].mask = maskInd
        # Global
        depthBini._FillValue = valmask
        depthBini = mv.masked_where(depthBini > 1.e6, depthBini)
        thickBini._FillValue = valmask
        thickBini = mv.masked_where(thickBini > 1.e6, thickBini)
        x1Bini._FillValue = valmask
        x1Bini = mv.masked_where(depthBini > 1.e6, x1Bini)
        x2Bini._FillValue = valmask
        x2Bini = mv.masked_where(depthBini > 1.e6, x2Bini)
        # Atl
        depthBinia._FillValue = valmask
        depthBinia = mv.masked_where(depthBinia > 1.e6, depthBinia)
        thickBinia._FillValue = valmask
        thickBinia = mv.masked_where(thickBinia > 1.e6, thickBinia)
        x1Binia._FillValue = valmask
        x1Binia = mv.masked_where(depthBinia > 1.e6, x1Binia)
        x2Binia._FillValue = valmask
        x2Binia = mv.masked_where(depthBinia > 1.e6, x2Binia)
        # Pac
        depthBinip._FillValue = valmask
        depthBinip = mv.masked_where(depthBinip > 1.e6, depthBinip)
        thickBinip._FillValue = valmask
        thickBinip = mv.masked_where(thickBinip > 1.e6, thickBinip)
        x1Binip._FillValue = valmask
        x1Binip = mv.masked_where(depthBinip > 1.e6, x1Binip)
        x2Binip._FillValue = valmask
        x2Binip = mv.masked_where(depthBinip > 1.e6, x2Binip)
        # Ind
        depthBinii._FillValue = valmask
        depthBinii = mv.masked_where(depthBinii > 1.e6, depthBinii)
        thickBinii._FillValue = valmask
        thickBinii = mv.masked_where(thickBinii > 1.e6, thickBinii)
        x1Binii._FillValue = valmask
        x1Binii = mv.masked_where(depthBinii > 1.e6, x1Binii)
        x2Binii._FillValue = valmask
        x2Binii = mv.masked_where(depthBinii > 1.e6, x2Binii)

        tozi = timc.clock()
        # 
        #
        # Compute zonal mean
        # Global
        depthBinz = cdu.averager(depthBini, axis=3)
        thickBinz = cdu.averager(thickBini, axis=3)
        x1Binz    = cdu.averager(x1Bini,    axis=3)
        x2Binz    = cdu.averager(x2Bini,    axis=3)
        # Atl
        depthBinza = cdu.averager(depthBinia, axis=3)
        thickBinza = cdu.averager(thickBinia, axis=3)
        x1Binza    = cdu.averager(x1Binia,    axis=3)
        x2Binza    = cdu.averager(x2Binia,    axis=3)
        # Pac
        depthBinzp = cdu.averager(depthBinip, axis=3)
        thickBinzp = cdu.averager(thickBinip, axis=3)
        x1Binzp    = cdu.averager(x1Binip,    axis=3)
        x2Binzp    = cdu.averager(x2Binip,    axis=3)
        # Ind
        depthBinzi = cdu.averager(depthBinii, axis=3)
        thickBinzi = cdu.averager(thickBinii, axis=3)
        x1Binzi    = cdu.averager(x1Binii,    axis=3)
        x2Binzi    = cdu.averager(x2Binii,    axis=3)

        toziz = timc.clock()
        # 
        #areaz , depthBinz, inv = ZonalMeans.compute(dy , area=area, delta_band=delta_lat)
        #areazt, thickBinz, inv = ZonalMeans.compute(ty , area=area, delta_band=delta_lat)
        #areaz , x1Binz   , inv = ZonalMeans.compute(x1y, area=area, delta_band=delta_lat)
        #areaz , x2Binz   , inv = ZonalMeans.compute(x2y, area=area, delta_band=delta_lat)
        #
        # Compute volume of isopycnals
        volBinz  = thickBinz  * areazt
        volBinza = thickBinza * areazta
        volBinzp = thickBinzp * areaztp
        volBinzi = thickBinzi * areazti
        #
        #
        # Init zonal mean output variables
        # Global
        dbz  = cdm.createVariable(depthBinz, axes = [dy.getAxis(0), s_axis, lati], id = 'isondepth')
        tbz  = cdm.createVariable(thickBinz, axes = [dy.getAxis(0), s_axis, lati], id = 'isonthick')
        vbz  = cdm.createVariable(volBinz*1.e-12, axes = [dy.getAxis(0), s_axis, lati], id = 'isonvol')
        x1bz = cdm.createVariable(x1Binz   , axes = [dy.getAxis(0), s_axis, lati], id = 'thetao')
        x2bz = cdm.createVariable(x2Binz   , axes = [dy.getAxis(0), s_axis, lati], id = 'so')
        # Atl
        dbza  = cdm.createVariable(depthBinza, axes = [dy.getAxis(0), s_axis, lati], id = 'isondeptha')
        tbza  = cdm.createVariable(thickBinza, axes = [dy.getAxis(0), s_axis, lati], id = 'isonthicka')
        vbza  = cdm.createVariable(volBinza*1.e-12, axes = [dy.getAxis(0), s_axis, lati], id = 'isonvola')
        x1bza = cdm.createVariable(x1Binza   , axes = [dy.getAxis(0), s_axis, lati], id = 'thetaoa')
        x2bza = cdm.createVariable(x2Binza   , axes = [dy.getAxis(0), s_axis, lati], id = 'soa')
        # Pac
        dbzp  = cdm.createVariable(depthBinzp, axes = [dy.getAxis(0), s_axis, lati], id = 'isondepthp')
        tbzp  = cdm.createVariable(thickBinzp, axes = [dy.getAxis(0), s_axis, lati], id = 'isonthickp')
        vbzp  = cdm.createVariable(volBinzp*1.e-12, axes = [dy.getAxis(0), s_axis, lati], id = 'isonvolp')
        x1bzp = cdm.createVariable(x1Binzp   , axes = [dy.getAxis(0), s_axis, lati], id = 'thetaop')
        x2bzp = cdm.createVariable(x2Binzp   , axes = [dy.getAxis(0), s_axis, lati], id = 'sop')
        # Ind
        dbzi  = cdm.createVariable(depthBinzi, axes = [dy.getAxis(0), s_axis, lati], id = 'isondepthi')
        tbzi  = cdm.createVariable(thickBinzi,  axes = [dy.getAxis(0), s_axis, lati], id = 'isonthicki')
        vbzi  = cdm.createVariable(volBinzi*1.e-12, axes = [dy.getAxis(0), s_axis, lati], id = 'isonvoli')
        x1bzi = cdm.createVariable(x1Binzi   , axes = [dy.getAxis(0), s_axis, lati], id = 'thetaoi')
        x2bzi = cdm.createVariable(x2Binzi   , axes = [dy.getAxis(0), s_axis, lati], id = 'soi')
        if tc == 0:
            # Global attributes
            dbz.long_name = 'Global zonal depth of isopycnal'
            dbz.units = 'm'
            tbz.long_name = 'Global zonal thickness of isopycnal'
            tbz.units = 'm'
            vbz.long_name = 'Volume of isopycnal'
            vbz.units = '10.e12 m^3'
            x1bz.long_name = temp.long_name
            x1bz.units = 'C'
            x2bz.long_name = so.long_name
            x2bz.units = so.units
            # Atl
            dbza.long_name = 'Atl. zonal depth of isopycnal'
            dbza.units = dbz.units
            tbza.long_name = 'Atl. zonal thickness of isopycnal'
            tbza.units = 'm'
            vbza.long_name = 'Atl. volume of isopycnal'
            vbza.units = '10.e12 m^3'
            x1bza.long_name = temp.long_name
            x1bza.units = 'C'
            x2bza.long_name = so.long_name
            x2bza.units = so.units
            # Pac
            dbzp.long_name = 'Pac. zonal depth of isopycnal'
            dbzp.units = dbz.units
            tbzp.long_name = 'Pac. zonal thickness of isopycnal'
            tbzp.units = 'm'
            vbzp.long_name = 'Pac. volume of isopycnal'
            vbzp.units = '10.e12 m^3'
            x1bzp.long_name = temp.long_name
            x1bzp.units = 'C'
            x2bzp.long_name = so.long_name
            x2bzp.units = so.units
            # Ind
            dbzi.long_name = 'Ind. zonal depth of isopycnal'
            dbzi.units = dbz.units
            tbzi.long_name = 'Ind. zonal thickness of isopycnal'
            tbzi.units = 'm'
            vbzi.long_name = 'Ind. volume of isopycnal'
            vbzi.units = '10.e12 m^3'
            x1bzi.long_name = temp.long_name
            x1bzi.units = 'C'
            x2bzi.long_name = so.long_name
            x2bzi.units = so.units
        # Write & append
        gz.write(dbz , extend = 1, index = (trmin-tmin)/12)
        gz.write(tbz , extend = 1, index = (trmin-tmin)/12)
        gz.write(vbz , extend = 1, index = (trmin-tmin)/12)
        gz.write(x1bz, extend = 1, index = (trmin-tmin)/12)
        gz.write(x2bz, extend = 1, index = (trmin-tmin)/12)
        # Atl
        gz.write(dbza , extend = 1, index = (trmin-tmin)/12)
        gz.write(tbza , extend = 1, index = (trmin-tmin)/12)
        gz.write(vbza , extend = 1, index = (trmin-tmin)/12)
        gz.write(x1bza, extend = 1, index = (trmin-tmin)/12)
        gz.write(x2bza, extend = 1, index = (trmin-tmin)/12)
        # Pac
        gz.write(dbzp , extend = 1, index = (trmin-tmin)/12)
        gz.write(tbzp , extend = 1, index = (trmin-tmin)/12)
        gz.write(vbzp , extend = 1, index = (trmin-tmin)/12)
        gz.write(x1bzp, extend = 1, index = (trmin-tmin)/12)
        gz.write(x2bzp, extend = 1, index = (trmin-tmin)/12)
        # Ind
        gz.write(dbzi , extend = 1, index = (trmin-tmin)/12)
        gz.write(tbzi , extend = 1, index = (trmin-tmin)/12)
        gz.write(vbzi , extend = 1, index = (trmin-tmin)/12)
        gz.write(x1bzi, extend = 1, index = (trmin-tmin)/12)
        gz.write(x2bzi, extend = 1, index = (trmin-tmin)/12)

    #    
    # Write/append to file
    if mthout == 0:
        g.write(depthBin, extend = 1, index = trmin-tmin)
        g.write(thickBin, extend = 1, index = trmin-tmin)
        g.write(x1Bin,    extend = 1, index = trmin-tmin)
        g.write(x2Bin,    extend = 1, index = trmin-tmin)
    #
    print '   CPU of density bining =', ticz-tuc
    if tcdel >= 12:
        print '   CPU of annual mean compute =', toz-ticz
        print '   CPU of persistence compute =', tozp-toz
        print '   CPU of interpolation =', tozi-tozp
        print '   CPU of zonal mean =', toziz-tozi

#
# end loop on tc <===
print
print ' Max memory use',resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1.e6,'GB'
ratio =  12.*float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)/float(grdsize*tmax)
print ' Ratio to grid*nyears',ratio,'kB/unit(size*nyears)'
#
ft.close()
fs.close()
if mthout == 0:
    g.close()
    print ' Wrote file: ', file_out
gz.close()
gp.close()
if tcdel >= 12:
    print ' Wrote file: ', filez_out
    print ' Wrote file: ', filep_out

print ' CPU use, elapsed', timc.clock() - ti0, timeit.default_timer() - te0
ratio = 1.e6*(timc.clock() - ti0)/float(grdsize*tmax)
print ' Ratio to grid*nyears',ratio,'1.e-6 sec/unit(size*nyears)'
print
# -----------------------------------------------------------------------------

# multiprocs:
# http://stackoverflow.com/questions/20820367/more-complex-parallel-for-loop-in-python

# ----------------------------------------
# some useful commands....

# array or tuple to list:
#array.tolist()
# find index: eg: where(s_s LT s_z[i_min])
#   ind = next(x[0] for x in enumerate(s_s) if x[1] < s_z[i_min])
#   see support procs: ind = sd.whereLT(s_s, s_z[i_min])


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
#  TO DO list:
#    - add bowl interpolation on density
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
import ZonalMeans
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
parser = argparse.ArgumentParser(description='Script to perform density bining analysis')
parser.add_argument('-d', help='toggle debug mode', action='count', default=0)
#parser.add_argument('-r','--sigma_range', help='neutral sigma range', required=True)
#parser.add_argument('-s','--sigma_increment', help='neutral sigma increment', required=True)
#parser.add_argument('-i','--input', help='input directory', default="./")
#parser.add_argument('-o','--output',help='output directory', default="./")
parser.add_argument('-t','--timeint', help='specify time domain in bining <init_idx>,<ncount>', default="all")
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
file_fx = '/work/cmip5/fx/fx/areacello/cmip5.IPSL-CM5A-LR.piControl.r0i0p0.fx.ocn.fx.areacello.ver-v20120430.latestX.xml'
file_T = '/work/cmip5/historical/ocn/mo/thetao/cmip5.IPSL-CM5A-LR.historical.r1i1p1.mo.ocn.Omon.thetao.ver-v20111119.latestX.xml'
file_S = '/work/cmip5/historical/ocn/mo/so/cmip5.IPSL-CM5A-LR.historical.r1i1p1.mo.ocn.Omon.so.ver-v20111119.latestX.xml'
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
temp = ft('thetao', time = slice(0,1))-273.15
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
if grdsize <= 1.e6:
    tcdel = min(120, tmax)
elif grdsize <= 1.e7:
    tcdel = min(12, tmax)
nyrtc = tcdel/12
tcmax = tmax/tcdel ; # number of time chunks
print '==> tcdel, tcmax:', tcdel, tcmax
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
file_out = outdir+'/out_density.nc'
if os.path.exists(file_out):
    os.remove(file_out)
g = cdm.open(file_out,'w+')
filez_out = outdir+'/outz_density.nc'
if os.path.exists(filez_out):
    os.remove(filez_out)
gz = cdm.open(filez_out,'w+')
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
#fileg='/work/guilyardi/database/ORAS4/ORAS4_1mm_01_12_1958-2009_grid1_so.nc'
#gt = cdm.open(fileg)
#sog = gt('so', time=slice(0,0))
#outgrid=sog.getGrid()
#gt.close()
#regridF = cdm.mvCdmsRegrid.regrid2.Regridder(ingrid,outgrid)
#outgrdz = cdm.createZonalGrid(ingrid)
#regridF = Regridder(ingrid,outgrdz)
#   
# loop on time chunks
for tc in range(tcmax):
    tuc = timc.clock()
    tuc2 = timeit.default_timer()
    # read tcdel month by tcdel month to optimise memory
    trmin = tc*tcdel ; # define as function of tc and tcdel
    trmax = (tc+1)*tcdel ; # define as function of tc and tcdel
    print ' --> time chunk (bounds) = ',tc, ' (',trmin,trmax-1,')'
    temp = ft('thetao', time = slice(trmin,trmax))-273.15
    so   = fs('so', time = slice(trmin,trmax))
    time  = temp.getTime()
    tur = timc.clock()
    print '     read  CPU:',tur-tuc
    # Compute neutral density (TODO optimize: 22 % CPU)
    rhon = sd.eos_neutral(temp,so)-1000.
    turn = timc.clock()
    print '     rhon compute CPU:',turn-tur
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
    #bowl_bin  = npy.ma.zeros([N_j, N_i]) # dim: i,j 
    tac = timc.clock()
    tac2 = timeit.default_timer()
    print '     read, rhon compute and array init CPU, elapsed:',tac-tuc, tac2-tuc2
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
        if debugp:
        #if t >= 0:
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
        # accumulate for annual mean
        #iyr = (tc*tcdel/12) + t/12*12
        #print iyr
        #dy[iyr,...] = dy[iyr,...] + depth_bin [t,...]
        #
        # debug
        if t == -1:
            ir=range(int(i_min[ijtest]),int(i_max[ijtest])+1)
            print 'test point',ijtest
            print ' i_bottom',i_bottom[ijtest]
            print ' i_min,i_max',i_min[ijtest],i_max[ijtest]
            #print ' ind',ind[0][npy.where(ind[1] == ijtest)]
            print ' i_profil',ir
            print ' s_z[i_profil] ', szm[ir,ijtest]
            print ' s_s[ind] ', s_s[ind[0][npy.where(ind[1]==ijtest)],ijtest]
            print ' z_zt[i_profil] ', zzm[ir,ijtest]
            print ' z_s[ind] ', z_s[ind[0][npy.where(ind[1] == ijtest)],ijtest]
            print ' c1_s[ind] ', c1_s[ind[0][npy.where(ind[1] == ijtest)],ijtest]
            print ' c2_s[ind] ', c2_s[ind[0][npy.where(ind[1] == ijtest)],ijtest]
    #
    # end of loop on t <===      
    #        
    # Reshape i*j back to i,j
    depth_bin = npy.reshape(depth_bin, (tcdel, N_s+1, N_j, N_i))
    thick_bin = npy.reshape(thick_bin, (tcdel, N_s+1, N_j, N_i))
    x1_bin    = npy.reshape(x1_bin,    (tcdel, N_s+1, N_j, N_i))
    x2_bin    = npy.reshape(x2_bin,    (tcdel, N_s+1, N_j, N_i))
    #
    # Wash mask over variables
    maskb = mv.masked_values(x1_bin, valmask).mask
    depth_bin.mask = maskb
    thick_bin.mask = maskb
    x1_bin.mask = maskb
    x2_bin.mask = maskb
    #
        
    #
    if tc == -1:
        # test write
        i = itest
        j = jtest
        print 'ind = ',ind
        print "test point",i,j, area[j,i]
        print "lon,lat",lon[j,i],lat[j,i]
        print 'depth_bin', depth_bin[0,:,j,i]
        print 'thick_bin', thick_bin[0,:,j,i]
        print 'x1_bin', x1_bin[0,:,j,i]
        print 'x2_bin', x2_bin[0,:,j,i]
    tic = timc.clock()
    tic2 = timeit.default_timer()
    print '   Loop on t done - CPU & elapsed total (per month) = ',tic-tac, tic2-tac2, '(',(tic-tac)/float(tcdel),(tic2-tac2)/float(tcdel),')'
    print '   Rhon computation vs. rest and % ',turn-tur,tic-tac,100.*(turn-tur)/(tic-tuc)
    #
    # Output files as netCDF
    # Def variables 
    # QQ??: only do for tc==0 ? depth_bin update enought for tc >= 1 ?
    depthBin = cdm.createVariable(depth_bin, axes = [time, s_axis, grd], id = 'isondepth')
    thickBin = cdm.createVariable(thick_bin, axes = [time, s_axis, grd], id = 'isonthick')
    x1Bin    = cdm.createVariable(x1_bin   , axes = [time, s_axis, grd], id = 'thetao')
    x2Bin    = cdm.createVariable(x2_bin   , axes = [time, s_axis, grd], id = 'so')
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
    # Compute annual mean, make zonal mean and write
    # TODO: optimize as VERY expensive (x2-4 preceeding loop !)
    ticz = timc.clock()
    if tcdel >= 12:
        # TODO: HUGE COST: 30-60 sec for 12 months !!! and 120 sec for 24 !!!
        dy = cdu.YEAR(depthBin)
        ty = cdu.YEAR(thickBin)
        x1y = cdu.YEAR(x1Bin)
        x2y = cdu.YEAR(x2Bin)
        # this is 5 times cheaper but no grid is passed and ZonalMeans fails
        #dy  = cdu.averager(npy.reshape (depthBin, (nyrtc, 12, N_s+1, N_j, N_i)), axis=1)
        #ty  = cdu.averager(npy.reshape (thickBin, (nyrtc, 12, N_s+1, N_j, N_i)), axis=1)
        #x1y = cdu.averager(npy.reshape (x1Bin,    (nyrtc, 12, N_s+1, N_j, N_i)), axis=1)
        #xy2 = cdu.averager(npy.reshape (x2Bin,    (nyrtc, 12, N_s+1, N_j, N_i)), axis=1)

        toz = timc.clock()
       # 10 sec for 12 months
        if debugp:
            print '   CPU of annual mean compute =', toz-ticz
            print '   test '
            print dy[0,:,80,60]
        areaz , depthBinz, inv = ZonalMeans.compute(dy , area=area, delta_band=delta_lat)
        areazt, thickBinz, inv = ZonalMeans.compute(ty , area=area, delta_band=delta_lat)
        areaz , x1Binz   , inv = ZonalMeans.compute(x1y, area=area, delta_band=delta_lat)
        areaz , x2Binz   , inv = ZonalMeans.compute(x2y, area=area, delta_band=delta_lat)
        # try to add degenerated x dimension for IDL read 
        #depthBinz.reshape = (nyrtc,180./delta_lat,1)
        #thickBinz.reshape = (nyrtc,180./delta_lat,1)
        #x1Binz.reshape    = (nyrtc,180./delta_lat,1)
        #x2Binz.reshape    = (nyrtc,180./delta_lat,1)
        if tc == 0:
            volBinz =  thickBinz*areazt
            volBinz.id='isonvol'
            volBinz.long_name = 'Volume of isopycnal'
            volBinz.units = 'm^3'
        gz.write(depthBinz, extend = 1, index = trmin/12)
        gz.write(thickBinz, extend = 1, index = trmin/12)
        gz.write(volBinz  , extend = 1, index = trmin/12)
        gz.write(x1Binz   , extend = 1, index = trmin/12)
        gz.write(x2Binz   , extend = 1, index = trmin/12)

    ticza = timc.clock()
    print '   CPU of zonal mean compute and write =', ticza-toz
    #    
    # Write/append to file
    g.write(depthBin, extend = 1, index = trmin)
    g.write(thickBin, extend = 1, index = trmin)
    g.write(x1Bin,    extend = 1, index = trmin)
    g.write(x2Bin,    extend = 1, index = trmin)
    #

#
# end loop on tc <===
print
print ' Max memory use',resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1.e6,'GB'
ratio =  12.*float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)/float(grdsize*tmax)
print ' Ratio to grid*nyears',ratio,'kB/unit(size*nyears)'
#
ft.close()
fs.close()
g.close()
gz.close()

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


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
#

import cdms2 as cdm
import MV2 as mv
import os, sys  
import socket, argparse
import string
import numpy as npy
import numpy.ma as ma
import cdutil as cdu
from genutil import statistics
import support_density as sd
import time as timc
import timeit
#import matplotlib.pyplot as plt

# netCDF compression (use 0 for netCDF3)
comp = 0
cdm.setNetcdfShuffleFlag(comp)
cdm.setNetcdfDeflateFlag(comp)
cdm.setNetcdfDeflateLevelFlag(comp)

#
cdm.setAutoBounds('on')
#
# == Arguments
#
# 
# == Inits
#

home='/Users/ericg/Projets/Density_bining'

if socket.gethostname() == 'crunchy.llnl.gov':
    home='/work/guilyardi/Density_bining'

hist_file_dir=home


# == get command line options
 
parser = argparse.ArgumentParser(description='Script to perform density bining analysis')
parser.add_argument('-d', help='toggle debug mode', action='count', default=0)
#parser.add_argument('-r','--sigma_range', help='neutral sigma range', required=True)
#parser.add_argument('-s','--sigma_increment', help='neutral sigma increment', required=True)
parser.add_argument('-i','--input', help='input directory', default="./")
parser.add_argument('-o','--output',help='output directory', default="./")
parser.add_argument('-t','--timeint', help='specify time domain in bining <init_idx>,<ncount>', default="all")
parser.add_argument('string', metavar='root for T and S files', type=str, help='netCDF input files root')
args = parser.parse_args()

# Write command line in history file

filer=hist_file_dir+'/z_density_hist.txt'

with open(filer, 'a') as f:
    f.write('\n\r'+str(sys.argv).translate(None, "',[]"))
 
## read values
debug        = str(args.d)
indir        = args.input
outdir       = args.output
#sigma_range  = args.sigma_range 
#delta_sigma  = args.sigma_increment
timeint      = args.timeint
file_root    = args.string

if debug >= '1': 
    print; print ' Debug - Args =', args

tic = timc.clock()

# Define T and S file names (local mac...)
file_T = indir+'/'+file_root+'_thetao.nc'
file_S = indir+'/'+file_root+'_so.nc'
file_fx = indir+'/areacello_fx_IPSL-CM5A-LR_piControl_r0i0p0.nc'

if socket.gethostname() == 'crunchy.llnl.gov':
    file_T = indir+'/thetao/cmip5.'+file_root+'.thetao.ver-v20111010.latestX.xml'
    file_S = indir+'/so/cmip5.'+file_root+'.so.ver-v20111010.latestX.xml'
    file_fx = '/work/cmip5/fx/fx/areacello/cmip5.IPSL-CM5A-LR.piControl.r0i0p0.fx.ocn.fx.areacello.ver-v20120430.latestX.xml'

if debug >= '1':
    print ' Debug - File names:'
    print '    ', file_T
    print '    ', file_S

# Open files

ft  = cdm.open(file_T)
fs  = cdm.open(file_S)

timeax = ft.getAxis('time')

# Dates to read

if timeint == 'all':
    tmin = 0
    tmax = timeax.shape[0]
else:
    tmin = int(timeint.split(',')[0]) - 1
    tmax = tmin + int(timeint.split(',')[1])

if debug >= '1':
    print; print ' Debug - Read only first month...'
    tmin = 0
    tmax = 1

print '  time interval: ', tmin, tmax - 1

# Define temperature and salinity arrays
# TODO: read month by month to optimise memory ?

temp = ft('thetao', time = slice(tmin,tmax))-273.15
so   = fs('so', time = slice(tmin,tmax))

# Read file attributes
list_file=ft.attributes.keys()
file_dic={}
for i in range(0,len(list_file)):
    file_dic[i]=list_file[i],ft.attributes[list_file[i] ]

# Read masking value
valmask = so._FillValue

# Read time and grid
time  = temp.getTime()
lon  = temp.getLongitude()
lat  = temp.getLatitude()

depth = temp.getLevel()

bounds = ft('lev_bnds')

# Read cell area

ff = cdm.open(file_fx)
area = ff('areacello')

toc = timc.clock()
print '   ... read CPU: ', toc-tic

ft.close()
fs.close()
ff.close()

# Define dimensions

N_i = int(lon.shape[1])
N_j = int(lon.shape[0])
N_z = int(depth.shape[0])
N_t = int(time.shape[0])

# Compute neutral density

rhon = sd.eos_neutral(temp,so)-1000.
tic = timc.clock()

# decide which field to bin
# TO DO: bring to args (or bin both T and S ?)

var='temp'
if var == 'temp':
    x1 = temp
    x1_name = 'Temperature'
    x1_units = temp.units

print; print 'Rhon computed (time = ',tic-toc,')'
print '  rho min/max ', npy.min(rhon), npy.max(rhon)

# Define sigma grid

rho_min = 19
rho_max = 28
del_s = 0.2
s_s = npy.arange(rho_min, rho_max, del_s).tolist()
N_s = len(s_s)

#w=sys.stdin.readline() # stop the code here. [Ret] to keep going

# start density bining (loop on t and horizontal grid)

imin = 0
jmin = 0

imax = temp.shape[3]
jmax = temp.shape[2]

# test point  
itest = 80 
jtest = 60  
           
#imin = itest
#jmin = jtest
#imax = itest+1
#jmax = jtest+1

# inits
# z profiles:
c1_z = [float('NaN')]*N_z
s_z  = [float('NaN')]*N_z
z_zt = depth[:]
z_zw = bounds.data[:,0]
#bowl_s = 0.

# density profiles:
z_s  = [float('NaN')]*(N_s+1)
c1_s = [float('NaN')]*(N_s+1)
z1_s = [float('NaN')]*(N_s+1)

# output arrays
depth_bin = npy.ma.ones([N_t, N_s+1, N_j, N_i], dtype='float32')*1.e+20 
depth_bin = mv.masked_where(depth_bin==1.e+20, depth_bin)
thick_bin = depth_bin.copy() 
vol_bin   = depth_bin.copy()
x1_bin    = depth_bin.copy() 
#bowl_bin  = npy.ma.zeros([N_j, N_i]) # dim: i,j 

print
toc = timc.clock()
toc2 = timeit.default_timer()


# loop on time
for t in range(tmin,tmax):
    print ' --> t = ',t
# TODO: read month by month to optimise memory ?
# x1 contents on vertical (not yet implemented - may be done to ensure conservation)
    x1_content = x1.data[t,:,:,:] # dims: i,j,k
    vmask_3D = mv.masked_values(temp[t, ...], 0) 

    # Loop on horizontal grid (TO DO: to be optimized !)
    # (TODO: reorganize arrays to collapse i j dims ? on keep only indices of ocean points ?)
    # (TODO: order of loops ok ?)
    for j in range(jmin,jmax):
        for i in range(imin,imax):
            # test on masked points
#            vmask = npy.nonzero(temp[t,:,j,i])
            vmask = vmask_3D[:,j,i].mask

            z_s = npy.asarray([float(0.0)]*(N_s+1)) 
            c1_s = npy.asarray([float('NaN')]*(N_s+1))
  
            #bowl_s = float('NaN') 

            if not vmask[0]: # check point is not masked
                # find bottom level
                i_bottom = npy.where(vmask)[0][0] - 1
                z_s[N_s] = z_zw[i_bottom+1]
                c1_s[N_s] = x1_content[N_z-1,j,i]

                s_z = rhon[t,:,j,i].data
                c1_z = x1_content[:,j,i]

                # extract a strictly increasing sub-profile
                # first test on bottom - surface stratification
                delta_rho = s_z[i_bottom] - s_z[0]
                if delta_rho < del_s:
                    i_min = 0
                    i_max = i_bottom
                else:
                    irange = range(i_bottom+1)
                    mini = min(s_z[irange])
                    maxi = max(s_z[irange])
                    i_min = (s_z[irange]).tolist().index(mini)
                    i_max = (s_z[irange]).tolist().index(maxi)

                if i_min > i_max:
                    print '*** i_min > i_max ', i_min,i_max
                    exit(1)

                # if level in s_s has lower density than surface, isopycnal is put at surface (z_s=0)
                ind = sd.whereLT(s_s, s_z[i_min])
                z_s[ind] = 0.
                c1_s[ind] = valmask

                # if level of s_s has higher density than bottom density, isopycnal is set to bottom (z_s=z_zw[i_bottom])
                ind = sd.whereGT(s_s, s_z[i_max])
                z_s[ind] = z_s[N_s]
                c1_s[ind] = valmask

                # General case
                ind = sd.where_between(s_s, s_z[i_min], s_z[i_max])
                if len(ind) >= 1:
                    i_profil = irange[i_min:i_max+1]

                # interpolate depth(z) (z_zt) to depth(s) at s_s densities (z_s) using density(z) s_z
                
                    z_s[ind] = npy.interp(npy.asarray(s_s)[ind], s_z[i_profil], z_zt[i_profil])
                        
                    c1_s[ind] = npy.interp(z_s[ind], z_zt[i_profil], c1_z[i_profil]) 

                    idt = sd.whereLT ( (z_s[1:N_s]-z_s[0:N_s-1]), -0.1 )
                    if len(idt) >= 1:
                        print 'ind = ',ind
                        print 'i_min,i_max  ', i_min,i_max
                        print 'i_profil ', i_profil
                        print "non increasing depth profile",i,j
                        print "lon,lat",lon[j,i],lat[j,i]
                        print " s_z[i_profil] ", s_z[i_profil]
                        print " z_s[ind] ", z_s[ind]
                        print " s_s[ind] ", npy.asarray(s_s)[ind]
                        print " z_zt[i_profil] ", z_zt[i_profil]

                
   #             print 'z_s = ',z_s
   #             print 'c1_s = ',c1_s

                # TO DO: bowl depth bining
                # IF sig_bowl EQ 1 THEN BEGIN
                #   bowl_s = interpol(s_z[i_profil], z_zt[i_profil], sobwlmax[i, j])
                # ENDIF ELSE BEGIN
                #   bowl_s = 0
                # ENDELSE
            # end if masked point

            depth_bin [t,:,j,i]     = z_s
            thick_bin [t,0,j,i]     = z_s[0]
            thick_bin [t,1:N_s,j,i] = z_s[1:N_s]-z_s[0:N_s-1]
            vol_bin   [t,:,j,i]     = thick_bin [t,:,j,i] * area[j,i]
#* e1t * e2t
            x1_bin    [t,:,j,i]     = c1_s
            #bowl_bin  [j, i]        = bowl_s

    # end loop on i,j

# end loop on t
#   
# test write
i = itest
j = jtest
print 'ind = ',ind
print "test point",i,j
print "lon,lat",lon[j,i],lat[j,i]
print 'thick_bin', thick_bin[0,:,j,i]
print 'vol_bin', vol_bin[0,:,j,i]
print 'x1_bin', x1_bin[0,:,j,i]

tic = timc.clock()
tic2 = timeit.default_timer()
print 'Loop on t,i,j done (CPU & elapsed = ',tic-toc, tic2-toc2, ')'

# Output files as netCDF

s_sd = npy.arange(rho_min, rho_max+del_s, del_s, dtype = npy.float32)
s_axis = cdm.createAxis(s_sd, id = 'rhon')
s_axis.long_name = 'Neutral density'
s_axis.units = ''
s_axis.designateLevel()

grd = temp.getGrid()

# Def variables

depthBin = cdm.createVariable(depth_bin, axes = [time, s_axis, grd], id = 'isondepth')
thickBin = cdm.createVariable(thick_bin, axes = [time, s_axis, grd], id = 'isonthick')
volBin   = cdm.createVariable(vol_bin, axes = [time, s_axis, grd], id = 'isonvol')
x1Bin    = cdm.createVariable(x1_bin, axes = [time, s_axis, grd], id = 'thetao')

depthBin.long_name = 'Depth of isopycnal'
depthBin.units = 'm'

thickBin.long_name = 'Thickness of isopycnal'
thickBin.units = 'm'

volBin.long_name = 'Volume of isopycnal'
volBin.units = 'm3'

x1Bin.long_name = 'Bined '+x1_name
x1Bin.units = x1_units

file_out = outdir+'/out_density.nc'
g = cdm.open(file_out,'w+')
g.write(depthBin)
g.write(thickBin)
g.write(volBin)
g.write(x1Bin)

# write global attributes (inherited from thetao file)
for i in range(0,len(file_dic)):
    dm=file_dic[i]
    setattr(g,dm[0],dm[1])
setattr(g,'Post-processing history','Density bining via densit_bin.py using delta_sigma = '+str(del_s))
g.close()

toc = timc.clock()
print 'Wrote file', file_out
print ' time write = ',toc-tic

# -----------------------------------------------------------------------------
# plt.plot(x, y, '.-')
# plt.plot(xs, ys)
# plt.show()
#next(x[0] for x in enumerate(s_s) if x[1] < s_z[i_min])
# use numpy first !                


# multiprocs:
# http://stackoverflow.com/questions/20820367/more-complex-parallel-for-loop-in-python

# ----------------------------------------
# some useful commands....

# d.info
# t=d.getTime() or t1=d.getAxis(0)
# print all: levs[:]
# import scipy
# import genutil (local dev)
# dir(genutil) doc for completion

#print so.attributes.keys()
#valmask = so._FillValue

# array or tuple to list:
#array.tolist()
# find index: eg: where(s_s LT s_z[i_min])
#   ind = next(x[0] for x in enumerate(s_s) if x[1] < s_z[i_min])
#   see support procs: ind = sd.whereLT(s_s, s_z[i_min])

#
# == detect time dimension and length
#
#time=ft[temp].getTime()

#if time is None:  
#    print "*** no time dimension in file ",file_T
#    count = raw_input("Enter number of time steps: ")
#else:
#    timename = time.id
#    count    = time.shape[0]



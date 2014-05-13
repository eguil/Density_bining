#!/usr/local/uvcdat/latest/bin/cdat
##!/Users/ericg/Projets/CMIP/Metrics/WGNE/bin/python
# 
# Program to compute density bins and replace vertical z coordinate by neutral density
# Reads in netCDF T(x,y,z,t) and S(x,y,z,t) files and writes 
#  - T or S(x,y,sigma,t)
#  - D(x,y,sigma,t) (depth of isopycnal)
#  - V(x,y,sigma,t) (volume of isopycnal)
#
#  TO DO list:
#    - add bowl interpolation on density
#
# Uses McDougall and Jackett 2005 (IDL routine provided by Gurvan Madec)
# Inspired from IDL density bin routines (by G. Roullet and G. Madec 1999)
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

if debug == '1': 
    print; print ' Debug - Args =', args

tic = timc.clock()

# Define T and S file names (local mac...)
file_T = indir+'/'+file_root+'_thetao.nc'
file_S = indir+'/'+file_root+'_so.nc'

if socket.gethostname() == 'crunchy.llnl.gov':
    file_T = indir+'/thetao/cmip5.'+file_root+'.thetao.ver-v20111010.latestX.xml'
    file_S = indir+'/so/cmip5.'+file_root+'.so.ver-v20111010.latestX.xml'

if debug == '1':
    print ' Debug - File names:',file_T, file_S

# Open files

ft  = cdm.open(file_T)
fs  = cdm.open(file_S)

# Define temperature and salinity arrays

if debug == '1':
    nread = 1
    print; print ' Debug - Read only first ',nread,' month(s)...'
    temp = ft('thetao',slice(0,nread-1))-273.15
    so   = fs('so', slice(0,nread-1))
else:
    temp = ft('thetao')-273.15
    so   = fs('so')

valmask = so._FillValue[0]

time  = temp.getTime()
lon  = temp.getLongitude()
lat  = temp.getLatitude()

depth = temp.getLevel()
bounds = ft('lev_bnds')

toc = timc.clock()
print ' ...time read = ', toc-tic

# Define dimensions

N_i = int(lon.shape[1])
N_j = int(lon.shape[0])
N_z = int(depth.shape[0])
N_t = int(time.shape[0])

# Compute neutral density

rhon = sd.eos_neutral(temp,so)-1000.
tic = timc.clock()

# decide which field to bin
# TO DO: bring to args

var='temp'
if var == 'temp':
    x1 = temp
    x1_name = 'Temperature'
    x1_units = temp.units

print
print 'Rhon computed (time = ',tic-toc,')'
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

imax = temp.shape[3] - 1
jmax = temp.shape[2] - 1

# test point                
#imin = 80
#jmin = 60
#imax = 80+1
#jmax = 60+1

if timeint == 'all':
    tmin = 0
    tmax = temp.shape[0] - 1
else:
    tmin = int(timeint.split(',')[0]) - 1
    tmax = tmin + int(timeint.split(',')[1])

print '  time interval: ', tmin, tmax - 1

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
x1_bin    = depth_bin.copy() 
#bowl_bin  = npy.ma.zeros([N_j, N_i]) # dim: i,j 

print
toc = timc.clock()
toc2 = timeit.default_timer()
# loop on time
for t in range(tmin,tmax):
    print ' --> t=',t
# x1 contents on vertical (not yet implemented - may be done to ensure conservation)
    x1_content = x1.data[t,:,:,:] # dims: i,j,k
    
    # Loop on horizontal grid (TO DO: to be optimized !)
    # (Paul: reorganize arrays to collapse i j dims ? on keep only indices of ocean points ?)
    # (Paul: order of loops ok ?)
    for j in range(jmin,jmax):
        for i in range(imin,imax):
            # loop on vertical axis to define in which density bins the vertical levels are
            # for k in range(temp.shape[1]):
            # test on masked points
            vmask = npy.nonzero(temp[t,:,j,i])
            z_s = npy.asarray([float(0.0)]*(N_s+1)) # simpler way to define an array of floats of dim N_s+1 ?
            c1_s = npy.asarray([float('NaN')]*(N_s+1))
            #bowl_s = float('NaN') 
            if len(vmask[0]) > 0: # i.e. point is not masked
                #
                i_bottom = vmask[0][len(vmask[0])-1]
                z_s[N_s] = z_zw[i_bottom]
                c1_s[N_s] = x1_content[N_z-1,j,i]

                s_z = rhon[t,:,j,i].data
                c1_z = x1_content[:,j,i]

                # extract a strictly increasing sub-profile
                # first test on bottom - surface stratification
                delta_rho = s_z[vmask[0]][i_bottom-1] - s_z[vmask[0]][0]
                if delta_rho < del_s:
                    i_min = 0
                    i_max = i_bottom
                else:
                    mini = min(s_z[vmask[0]])
                    maxi = max(s_z[vmask[0]])
                    i_min = (s_z[vmask[0]]).tolist().index(mini)
                    i_max = (s_z[vmask[0]]).tolist().index(maxi)

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
                c1_s[ind] = c1_s[N_s]
                c1_s[ind] = valmask

                # General case
                ind = sd.where_between(s_s, s_z[i_min], s_z[i_max])
                if len(ind) >= 1:
                    i_profil = vmask[0][i_min:i_max]

                # interpolate depth(z) (z_zt) to depth(s) at s_s densities (z_s) using density(z) s_z
#                    if delta_rho < del_s:
#                        print 'ind = ',ind
#                        print "test point",i,j
#                        print "lon,lat",lon[j,i],lat[j,i]
#                        print 'i_min,i_max ', i_min,i_max
#                        print 'density profile s_z', s_z
#                        print 'delta_rho ', delta_rho
#                        print ' '
                #                    print 'field profile c1_z', c1_z
                
                    z_s[ind] = npy.interp(npy.asarray(s_s)[ind], s_z[i_profil], z_zt[i_profil])
                        
                    c1_s[ind] = npy.interp(z_s[ind], z_zt[i_profil], c1_z[i_profil]) 
                
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
            x1_bin    [t,:,j,i]     = c1_s
            #bowl_bin  [j, i]        = bowl_s

    # end loop on i,j

# end loop on t

tic = timc.clock()
tic2 = timeit.default_timer()
print 'Loop on t,i,j done (CPU & elapsed = ',tic-toc, tic2-toc2, ')'

# Output files as netCDF

s_sd = npy.arange(rho_min, rho_max+del_s, del_s, dtype=npy.float32)
s_axis = cdm.createAxis(s_sd)
s_axis.id = 'Neutral_density'
s_axis.units = ''
s_axis.designateLevel()

# Def variables

depthBin = cdm.createVariable(depth_bin, axes=[time, s_axis, lat, lon])
thickBin = cdm.createVariable(thick_bin, axes=[time, s_axis, lat, lon])
x1Bin    = cdm.createVariable(x1_bin, axes=[time, s_axis, lat, lon])

depthBin.id = 'isodepth'
depthBin.long_name = 'Depth of isopycnal'
depthBin.units = 'm'

thickBin.id = 'isothick'
thickBin.long_name = 'Thickness of isopycnal'
thickBin.units = 'm'

x1Bin.id = 'thetao'
x1Bin.long_name = 'Bined '+x1_name
x1Bin.units = x1_units

file_out = outdir+'/density_out.nc'
g = cdm.open(file_out,'w+')
g.write(depthBin)
g.write(thickBin)
g.write(x1Bin)
#g.description = 'Density bining via densit_bin.py using delta_sigam = ', del_s
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



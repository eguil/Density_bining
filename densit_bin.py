#!/Users/ericg/Projets/CMIP/Metrics/WGNE/bin/python
# 
# Program to compute density bins and replace vertical z coordinate by neutral density
# Reads in netCDF T(x,y,z,t) and S(x,y,z,t) files and writes 
#  - T or S(x,y,sigma,t)
#  - D(x,y,sigma,t) (depth of isopycnal)
#  - V(x,y,sigma,t) (volume of isopycnal)
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

#
# == Arguments
#
# 
# == Inits
#

home='/Users/ericg/Projets/Density_bining'
hist_file_dir=home

if socket.gethostname() == 'crunchy.llnl.gov':
    home='/work/guilyardi'

hist_file_dir=home
toolpath=home+"/STL_analysis"

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

with open(filer, "a") as f:
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
    print args


# Define T and S file names
file_T=file_root+'_thetao.nc'
file_S=file_root+'_so.nc'


if debug == "1":
    print file_T, file_S
  
ft  = cdm.open(indir+"/"+file_T)
fs  = cdm.open(indir+"/"+file_S)

# Define temperature and salinity arrays

temp = ft("thetao")-273.15
so   = fs("so")
valmask = so._FillValue[0]

time  = temp.getTime()
lon  = temp.getLongitude()
lat  = temp.getLatitude()

depth = temp.getLevel()
bounds = ft("lev_bnds")

N_z = int(depth.shape[0])

# Compute neutral density

rhon = sd.eos_neutral(temp,so)-1000.

# decide which field to bin

x1 = temp

print 'Rhon computed'

# Define sigma grid

s_s = npy.arange(19,28,.2).tolist()
N_s = len(s_s)

#w=sys.stdin.readline() # stop the code here. [Ret] to keep going

# start density bining (loop on t and horizontal grid)

imin = 1
jmin = 1
tmin = 1
imax = temp.shape[2]
jmax = temp.shape[3]
kmax = temp.shape[1]
tmax = temp.shape[0]

# test point                
imin = 60
jmin = 80
imax = 60+1
jmax = 80+1
tmax = 1+1

# inits
# z profiles:
c1_z = [float('NaN')]*N_z
s_z  = [float('NaN')]*N_z
z_zt = depth[:]
z_zw = bounds.data[:,0]
bowl_s = 0.

# density profiles:
z_s  = [float('NaN')]*(N_s+1)
c1_s = [float('NaN')]*(N_s+1)
z1_s = [float('NaN')]*(N_s+1)

# output arrays
depth_bin = [] # dim: i,j,Ns_+1
thick_bin = [] # dim: i,j,Ns_+1
x1_bin    = [] # dim: i,j,Ns_+1
bowl_bin  = [] # dim: i,j 

# loop on time
for t in range(tmin,tmax):

# x1 contents on vertical
    x1_content = x1[t,:,:,:] # dims: i,j,k
    
    # Loop on horizontal grid (to be optimized !)
    # (Paul: reorganize arrays to collapse i j dims ? on keep only indices of ocean points ?)
    # (Paul: dims are always this order ?)
    for j in range(jmin,jmax):
        for i in range(imin,imax):
            # loop on vertical axis to define in which density bins the vertical levels are
            # for k in range(temp.shape[1]):
            # test on masked points
            vmask = npy.nonzero(temp[t,:,j,i])
            z_s = npy.asarray([float(0.0)]*(N_s+1))
            c1_s = npy.asarray([float('NaN')]*(N_s+1))
            bowl_s = float("NaN") 
            if len(vmask[0]) > 0: # i.e. point is not masked
                print "test point",i,j
                print "lon,lat",lon[j,i],lat[j,i]
                print "depths ",depth[:]
                #
                i_bottom = vmask[0][len(vmask[0])-1]
                z_s[N_s] = z_zw[i_bottom]
                c1_s[N_s] = x1_content[N_z-1,j,i]

                s_z = rhon[t,:,j,i].data
                c1_z = x1_content[:,j,i]

                print 'density profile s_z', s_z
                print 'field profile c1_z', c1_z

                # extract a strictly increasing sub-profile
                mini = min(s_z[vmask[0]])
                maxi = max(s_z[vmask[0]])
                i_min = (s_z[vmask[0]]).tolist().index(mini)
                i_max = (s_z[vmask[0]]).tolist().index(maxi)

                # largest of min indices and smaller of max indices
                # to do

                # if level in s_s has lower density than surface, isopycnal is put at surface (z_s=0)
                ind = sd.whereLT(s_s, s_z[i_min])
                z_s[ind] = 0.
                c1_s[ind] = 0.

                # if level of s_s has higher density than bottom density, isopycnal is set to bottom (z_s=z_zw[i_bottom])
                ind = sd.whereGT(s_s, s_z[i_max])
                z_s[ind] = z_s[N_s]
                c1_s[ind] = c1_s[N_s]
                c1_s[ind] = float("NaN")

                # General case
                ind = sd.where_between(s_s, s_z[i_min], s_z[i_max])

                print 'ind = ',ind


#next(x[0] for x in enumerate(s_s) if x[1] < s_z[i_min])
                

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



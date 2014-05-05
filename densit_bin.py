#!/Users/ericg/Projets/CMIP/Metrics/WGNE/bin/python
# 
# Program to compute density bins and replace vertical z coordinate by neutral density
# Reads in netCDF T(x,y,z,t) and S(x,y,z,t) files and writes 
#  - T(x,y,sigma,t)
#  - S(x,y,sigma,t)
#
# Uses McDougall and Jackett 2005 (IDL routine provided by Gurvan Madec)
# Inspired from IDL density bin routines
#
# --------------------------------------------------------------
#  E. Guilyardi - while at LBNL/LLNL -  March 2014
# 
# git add densit_bin.py
# git commit -m "with Paul"
# git push
#

import cdms2 as cdm
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

temp = ft("thetao")
so = fs("so")

print so.attributes.keys()
valmask = so._FillValue

w=sys.stdin.readline() # stop the code here. [Ret] to keep going

maskVar=ma.masked_greater_equal(temp, valmask).mask
mask=1.-maskVar.astype(float)
rhon=sd.eos_neutral(temp,so,mask)

# Define sigma grid

sig = np.arange(19,28,.2)

#
# == detect time dimension and length
#
time=f[temp].getTime()

if time is None:  
    print "*** no time dimension in file ",file_T
    count = raw_input("Enter number of time steps: ")
else:
    timename = time.id
    count    = time.shape[0]


# target grid
fileg='/Users/ericg/Desktop/Data/ORAS4/ORAS4_1mm_01_12_1958-2009_grid1_so.nc'
g = cdm.open(fileg)
so = g('so', time=slice(0,0))
grd=so.getGrid()
#w=sys.stdin.readline() # stop the code here. [Ret] to keep going

for l in range(len(lst)):
  i=lst[l]
  mod = string.split(i,'.')[1]
  f = cdms.open(pathin + var + '/' + i)
  d = f(var)
  d=d[0,0,...] # needs to loop on vertical levels

# d=f(var, latitude=(-5.,5.),longitude=(-140.,-110.),level=(0,40))
#  print d.info
# interpolate - regrid
#  if rgrid_meth == 'regrid2':
#    do = do.regrid(obsg,regridTool='regrid2')
#  if rgrid_meth in ['linear','conservative']:
  diag = {}
  d = d.regrid(grd,regridTool='ESMF',regridMethod='linear')
  aavg = cdu.averager(d, axis = '123')
  stdev=float(statistics.std(aavg))

  print mod,'  ', aavg.shape, ' ',stdev 

  f.close()

  aavg.id=var
  g=cdms.open(mod+"_out.nc","w+")
  g.write(aavg)
  g.close()



# d.info
# t=d.getTime() or t1=d.getAxis(0)
# print all: levs[:]
# import scipy
# import genutil (local dev)
# dir(genutil) doc for completion




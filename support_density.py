#!/Users/ericg/Projets/CMIP/Metrics/WGNE/bin/python
# 
#  ---------------------------------------
#  Support python procs for density bining
#  ---------------------------------------
#  EG May 2014
#
# -------------------------------------------------------------------------------
#
# NAME: EOS
#
# PURPOSE:
#       computes rho  (in situ volumic mass) 
#
# CATEGORY:
#       calculation
#
# CALLING SEQUENCE:
#       tableau=eos(t,s,tmask)
#
# INPUTS:
#       t : temperature
#       s : salinity
#       mask
#
# REFERENCE:
#	Compute the neutral volumic mass (Kg/m3) from known potential 
#       temperature and salinity fields using McDougall and Jackett 2005
#       equation of state.
#              potential temperature         t        deg celsius
#              salinity                      s        psu
#              nutral density                rho      kg/m**3
#
#         Check value: rho(35,20) = 1024.59416751197 kg/m**3 
#          t = 20 deg celcius, s=35 psu
#
#       McDougall and Jackett, J. Mar Res., 2005
#
#   MODIFICATION HISTORY:
#-      Gurvan Madec (04/14/2005) 
#-      Eric Guilyardi (02/05/2014) python version
#------------------------------------------------------------
#------------------------------------------------------------
#------------------------------------------------------------
import numpy as npy
import MV2 as mv

def eos_neutral(t, s):
#
# mask t and s fields

    zt=t
    zs=s
#
# neutral density
#
#   square root salinity
     	
    zsr= npy.ma.sqrt(zs)

# Numerator
# T-Polynome:                    T^3                       T^2                         T                      cst
    zr1= ( ( -4.3159255086706703e-4*zt+8.1157118782170051e-2 )*zt+2.2280832068441331e-1 )*zt+1002.3063688892480
# S-T Polynome:                  S^2                        T S                     S
    zr2= ( -1.7052298331414675e-7*zs-3.1710675488863952e-3*zt-1.0304537539692924e-4 )*zs

#Denominator
# T-Polynome:                     T^4                      T^3                       T^2                         T                      cst
    zr3= ( ( (-2.3850178558212048e-9*zt -1.6212552470310961e-7 )*zt+7.8717799560577725e-5 )*zt+4.3907692647825900e-5 )*zt+     1.0
# S T-Polynome:                     T^3 S                  T S                       S
    zr4= ( ( -2.2744455733317707e-9*zt*zt+6.0399864718597388e-6)*zt-5.1268124398160734e-4 )*zs
# S T-Polynome:                     T^2 S^3/2               S^3/2
    zr5= ( -1.3409379420216683e-9*zt*zt-3.6138532339703262e-5)*zs*zsr
#
#   ... masked neutral density
    zrho= ( zr1 + zr2 ) / ( zr3 + zr4 + zr5 )
    return zrho 

def whereEQ(elements, reference):
    # find index(indices)
    index_list=[]
    for index, elem in enumerate(elements):
        if elem == reference:
            index_list.append(index)
    return npy.asarray(index_list).astype(int)
    raise ValueError("No index found whereEQ")


def whereLT(elements, reference):
    # find index(indices)
    index_list=[]
    for index, elem in enumerate(elements):
        if elem < reference:
            index_list.append(index)
    return npy.asarray(index_list).astype(int)
    raise ValueError("No index found whereLT")

def whereGT(elements, reference):
    # find index(indices)
    index_list=[]
    for index, elem in enumerate(elements):
        if elem > reference:
            index_list.append(index)
    return npy.asarray(index_list).astype(int)
    raise ValueError("No index found whereGT")

def where_between(elements, reference_low, reference_high):
    # find index(indices)
    index_list=[]
    for index, elem in enumerate(elements):
        if elem > reference_low and elem < reference_high:
            index_list.append(index)
    return npy.asarray(index_list).astype(int)
    raise ValueError("No index found where_between")

def interp_mask(ifield, ofield, regrido, maski):
    ofield = regrido(ifield)
    ofield.mask = maski
    return ofield

def mask_val(field, valmask):
    field [npy.isnan(field.data)] = valmask
    field._FillValue = valmask
    field = mv.masked_where(field > valmask/10, field)
    return field

def compute_area(lon, lat):
    # compute area of grid cells on earth
    # use mid points and formula:
    # area = R^2(lon2-lon1)*(sin(lat2) - sin(lat1))
    rade = 6371000.
    radconv = npy.pi/180.
    N_i = int(lon.shape[0])
    N_j = int(lat.shape[0])
    area = npy.ma.ones([N_j, N_i], dtype='float32')*0.
    lonr = lon[:] * radconv
    latr = lat[:] * radconv
    #loop
    for i in range(1,N_i-1):
        lonm1 = (lonr[i-1] + lonr[i]  )*0.5
        lonp1 = (lonr[i]   + lonr[i+1])*0.5
        for j in range(1,N_j-1):
            latm1 = (latr[j-1] + latr[j]  )*0.5
            latp1 = (latr[j]   + latr[j+1])*0.5
            area[j,i] = npy.float(rade**2 * (lonp1 - lonm1) * (npy.sin(latp1) - npy.sin(latm1)))
        # North and south bounds
        latm1 = ((-90.*radconv) + latr[0] )*0.5
        latp1 = (latr[0]        + latr[1] )*0.5
        area[0,i] = npy.float(rade**2 * (lonp1 - lonm1) * (npy.sin(latp1) - npy.sin(latm1)))
        latm1 = (latr[N_j-2] + latr[N_j-1])*0.5
        latp1 = (latr[N_j-1] + (90.*radconv)  )*0.5
        area[N_j-1,i] = npy.float(rade**2 * (lonp1 - lonm1) * (npy.sin(latp1) - npy.sin(latm1)))
    # East and west bounds
    area[:,0]     = area[:,1]
    area[:,N_i-1] = area[:,N_i-2] 

    return area

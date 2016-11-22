
def defModels():
#
# List available models and properties
#
#  name, props=[Nb of hist members, nb of HistNat members, idx of common interval [1861-2005]]
#      picontrol=[length of run], correctFile [idx_i,idx_i1,jmax] for longitude correction
#

    models = [
#        {'name':'ACCESS1-0'     ,'props':[2,0,11,156], 'picontrol':[500],'correctFile':[0,0,0]}, # 0
#        {'name':'ACCESS1-3'     ,'props':[3,0,11,156], 'picontrol':[500],'correctFile':[0,0,0]}, # 1
#        {'name':'bcc-csm1-1-m'  ,'props':[3,0,11,156], 'picontrol':[0],'correctFile':[0,0,0]}, # 2
         {'name':'bcc-csm1-1'    ,'props':[3,1,11,156], 'picontrol':[0],'correctFile':[0,0,0]}, # 3
#        {'name':'bcc-csm1-1'    ,'props':[3,1,11,16], 'picontrol':[0],'correctFile':[0,0,0]}, # 3
#        {'name':'BNU-ESM'       ,'props':[1,0,11,156], 'picontrol':[559],'correctFile':[0,0,0]}, # 4
         {'name':'CanESM2'       ,'props':[5,5,11,156], 'picontrol':[996],'correctFile':[179,180,180]}, # 5
         {'name':'CCSM4'         ,'props':[6,4,11,156], 'picontrol':[1051],'correctFile':[139,140,145]}, # 6
#        {'name':'CESM1-BGC'     ,'props':[1,0,11,156], 'picontrol':[500],'correctFile':[0,0,0]}, # 7
         {'name':'CESM1-CAM5'    ,'props':[3,2,11,156], 'picontrol':[319],'correctFile':[139,140,145]}, # 8
#        {'name':'CESM1-FASTCHEM','props':[3,0,11,156], 'picontrol':[175],'correctFile':[0,0,0]}, # 9
#        {'name':'CESM1-WACCM'   ,'props':[1,0,11,156], 'picontrol':[200],'correctFile':[0,0,0]}, # 10
#        {'name':'CMCC-CESM'     ,'props':[1,0,11,156], 'picontrol':[277],'correctFile':[0,0,0]}, # 11
#        {'name':'CMCC-CM'       ,'props':[1,0,11,156], 'picontrol':[330],'correctFile':[0,0,0]}, # 12
#        {'name':'CMCC-CMS'      ,'props':[1,0,11,156], 'picontrol':[500],'correctFile':[0,0,0]}, # 13
#        {'name':'CNRM-CM5-2'    ,'props':[1,0,11,156], 'picontrol':[410],'correctFile':[0,0,0]}, # 14
         {'name':'CNRM-CM5'      ,'props':[9,6,11,156], 'picontrol':[850],'correctFile':[0,0,0]}, # 15
         {'name':'CSIRO-Mk3-6-0' ,'props':[9,5,11,156], 'picontrol':[500],'correctFile':[178,179,180]}, # 16
#        {'name':'CSIRO-Mk3L-1-2','props':[2,0,10,155], 'picontrol':[0],'correctFile':[0,0,0]}, # 17
#        {'name':'EC-EARTH'      ,'props':[6,0,11,156], 'picontrol':[452],'correctFile':[0,0,0]}, # 18
         {'name':'FGOALS-g2'     ,'props':[4,3,11,156], 'picontrol':[700],'correctFile':[0,0,0]}, # 19
#        {'name':'FGOALS-s2'     ,'props':[3,0,11,156], 'picontrol':[501],'correctFile':[0,0,0]}, # 20
#        {'name':'GFDL-CM2p1'    ,'props':[9,0,11,156], 'picontrol':[0],'correctFile':[0,0,0]}, # 21
         {'name':'GFDL-CM3'      ,'props':[4,3, 1,146], 'picontrol':[0],'correctFile':[0,0,0]}, # 22
#        {'name':'GFDL-ESM2G'    ,'props':[1,0, 1,146], 'picontrol':[500],'correctFile':[0,0,0]}, # 23
         {'name':'GFDL-ESM2M'    ,'props':[1,1, 0,145], 'picontrol':[500],'correctFile':[0,0,0]}, # 24
#         {'name':'GISS-E2-H'     ,'props':[0,11,11,156],'picontrol':[780],'correctFile':[0,0,0]},# 25
#        {'name':'GISS-E2-H-CC'  ,'props':[1,0,11,156], 'picontrol':[251],'correctFile':[0,0,0]}, # 26
         {'name':'GISS-E2-R'     ,'props':[16,11,11,156],'picontrol':[846],'correctFile':[179,180,180]},# 27
#        {'name':'GISS-E2-R-CC'  ,'props':[1,0,11,156], 'picontrol':[251],'correctFile':[0,0,0]}, # 28
#        {'name':'HadCM3'        ,'props':[9,0, 1,146], 'picontrol':[0],'correctFile':[0,0,0]}, # 29
#        {'name':'HadGEM2-CC'    ,'props':[1,0, 1,146], 'picontrol':[240],'correctFile':[0,0,0]}, # 30
         {'name':'HadGEM2-ES'    ,'props':[3,3, 1,146], 'picontrol':[576],'correctFile':[179,179,180]}, # 31
         {'name':'IPSL-CM5A-LR'  ,'props':[6,3,11,156], 'picontrol':[1000],'correctFile':[0,0,0]}, # 32
         {'name':'IPSL-CM5A-MR'  ,'props':[3,3,11,156], 'picontrol':[300],'correctFile':[0,0,0]}, # 33
#        {'name':'IPSL-CM5B-LR'  ,'props':[1,0,11,156], 'picontrol':[300],'correctFile':[0,0,0]}, # 34
         {'name':'MIROC-ESM-CHEM','props':[1,1,11,156], 'picontrol':[255],'correctFile':[179,180,180]}, # 35
         {'name':'MIROC-ESM'     ,'props':[3,3,11,156], 'picontrol':[100],'correctFile':[179,180,180]}, # 36
#        {'name':'MPI-ESM-LR'    ,'props':[3,0,11,156], 'picontrol':[0],'correctFile':[0,0,0]}, # 37
#        {'name':'MPI-ESM-MR'    ,'props':[3,0,11,156], 'picontrol':[0],'correctFile':[0,0,0]}, # 38
#        {'name':'MPI-ESM-P'     ,'props':[2,0,11,156], 'picontrol':[0],'correctFile':[0,0,0]}, # 39
#        {'name':'NorESM1-ME'    ,'props':[1,0,11,156], 'picontrol':[0],'correctFile':[0,0,0]}, # 40
#        {'name':'NorESM1-M'     ,'props':[2,0,11,156], 'picontrol':[0],'correctFile':[0,0,0]}, # 41
        ]

    return models

# models to check
"""
    models = [
        ]
"""

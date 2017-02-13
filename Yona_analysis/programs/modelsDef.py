
def defModels():
#
# List available models and properties
#
#  name, props=[Nb of hist members, nb of HistNat members, idx of common interval [1861-2005]]
#      picontrol=[length of run], correctFile [idx_i,idx_i1,jmax] for longitude correction
#

    models = [
        {'name':'bcc-csm1-1'    ,'props':[3,1,11,156], 'picontrol':[0],'correctFile':[0,0,0], 'file_end':'v20130329'}, #0
         {'name':'CanESM2'       ,'props':[5,5,11,156], 'picontrol':[996],'correctFile':[179,180,180], 'file_end':'1'}, #1
         {'name':'CCSM4'         ,'props':[6,4,11,156], 'picontrol':[1051],'correctFile':[139,140,145], 'file_end':'v20121128'}, #2
         {'name':'CESM1-CAM5'    ,'props':[3,2,11,156], 'picontrol':[319],'correctFile':[139,140,145], 'file_end':'v20140822'}, #3
         {'name':'CNRM-CM5'      ,'props':[9,6,11,156], 'picontrol':[850],'correctFile':[0,0,0], 'file_end':'v20130101'}, #4
         {'name':'CSIRO-Mk3-6-0' ,'props':[9,5,11,156], 'picontrol':[500],'correctFile':[178,179,180], 'file_end':'1'}, #5
         {'name':'FGOALS-g2'     ,'props':[4,3,11,156], 'picontrol':[700],'correctFile':[0,0,0], 'file_end':'v1'}, #6
         {'name':'GFDL-CM3'      ,'props':[4,3, 1,146], 'picontrol':[0],'correctFile':[0,0,0], 'file_end':'v20110601'}, #7
         {'name':'GFDL-ESM2M'    ,'props':[1,1, 0,145], 'picontrol':[500],'correctFile':[0,0,0], 'file_end':'v20130226'}, #8
         {'name':'GISS-E2-R'     ,'props':[16,11,11,156],'picontrol':[846],'correctFile':[179,180,180], 'file_end':'v20121015'},#9
         {'name':'HadGEM2-ES'    ,'props':[3,3, 1,146], 'picontrol':[576],'correctFile':[179,179,180], 'file_end':'v20110916'}, #10
         {'name':'IPSL-CM5A-LR'  ,'props':[6,3,11,156], 'picontrol':[1000],'correctFile':[0,0,0], 'file_end':'v20111119'}, #11
         {'name':'IPSL-CM5A-MR'  ,'props':[3,3,11,156], 'picontrol':[300],'correctFile':[0,0,0], 'file_end':'v20120804'}, #12
         {'name':'MIROC-ESM-CHEM','props':[1,1,11,156], 'picontrol':[255],'correctFile':[179,180,180], 'file_end':'1'}, #13
         {'name':'MIROC-ESM'     ,'props':[3,3,11,156], 'picontrol':[100],'correctFile':[179,180,180], 'file_end':'1'}, #14
        ]

    return models


def defModelsCO2piC():
    models = [
        {'name':'ACCESS1-0'     , 'picontrol':[500],'correctFile':[0,0,0], 'file_end_CO2':'1',
         'file_end_piC':'1'}, # 0
        {'name':'ACCESS1-3'     , 'picontrol':[500],'correctFile':[0,0,0], 'file_end_CO2':'1',
         'file_end_piC':'1'}, # 1
        {'name':'BNU-ESM'       , 'picontrol':[559],'correctFile':[0,0,0], 'file_end_CO2':'1',
         'file_end_piC':'1'}, # 2
        {'name':'CCSM4'         , 'picontrol':[1051],'correctFile':[139,140,145], 'file_end_CO2':'v20121128',
         'file_end_piC':'v20130513'}, # 3
        {'name':'CESM1-BGC'     , 'picontrol':[500],'correctFile':[0,0,0], 'file_end_CO2':'v20140822',
         'file_end_piC':'v20140822'}, # 4
        {'name':'CESM1-CAM5'    , 'picontrol':[319],'correctFile':[139,140,145], 'file_end_CO2':'v20121129',
         'file_end_piC':'v20140822'}, #5
        {'name':'CNRM-CM5'      , 'picontrol':[850],'correctFile':[0,0,0], 'file_end_CO2':'v20130101',
         'file_end_piC':'v20121001'}, # 6
        {'name':'CNRM-CM5-2'    , 'picontrol':[410],'correctFile':[0,0,0], 'file_end_CO2':'v20130401',
         'file_end_piC':'v20130402'}, # 7
        {'name':'CSIRO-Mk3-6-0' , 'picontrol':[500],'correctFile':[178,179,180], 'file_end_CO2':'v20111221',
         'file_end_piC':'1'}, # 8
        {'name':'FGOALS-g2'     , 'picontrol':[700],'correctFile':[0,0,0], 'file_end_CO2':'1',
         'file_end_piC':'v1'}, # 9
        {'name':'GFDL-ESM2G'    , 'picontrol':[500],'correctFile':[0,0,0], 'file_end_CO2':'v20120820',
         'file_end_piC':'v20110601'}, # 10
        {'name':'GFDL-ESM2M'    , 'picontrol':[500],'correctFile':[0,0,0], 'file_end_CO2':'v20130226',
         'file_end_piC':'v20130226'}, # 11
        {'name':'GISS-E2-H'     ,'picontrol':[780],'correctFile':[0,0,0], 'file_end_CO2':'v20130925',
         'file_end_piC':'1'}, # 12
        {'name':'HadGEM2-ES'    , 'picontrol':[576],'correctFile':[179,179,180], 'file_end_CO2':'v20111017',
         'file_end_piC':'v20110928'}, # 13
        {'name':'IPSL-CM5A-LR'  , 'picontrol':[1000],'correctFile':[0,0,0], 'file_end_CO2':'v20120114.latestX.WARN2.xml',
         'file_end_piC':'v20111010'}, # 14
        {'name':'IPSL-CM5A-MR'  , 'picontrol':[300],'correctFile':[0,0,0], 'file_end_CO2':'v20111119',
         'file_end_piC':'v20111119'}, # 15
        {'name':'IPSL-CM5B-LR'  , 'picontrol':[300],'correctFile':[0,0,0], 'file_end_CO2':'v20120430',
         'file_end_piC':'v20120114'}, # 16
    ]

    return models

# models to check
"""
    models = [
        ]
"""

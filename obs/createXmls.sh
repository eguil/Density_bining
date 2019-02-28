#!/bin/bash
# -*- coding: utf-8 -*-
#"""
#Created on Sun Apr 24 11:31:35 2016
#
# PJD 24 Apr 2016   - Added host check
# PJD 24 Apr 2016   - Using 'include' var argument
# PJD 24 Apr 2016   - Moved xml cleanup into if blocks
# PJD 24 Apr 2016   - Further tweaked arg to deal with ORAS4 garbled CF-1 inputs
# PJD 24 Apr 2016   - Added chdir
# PJD 24 Apr 2016   - Updated SODA224 to use --time-linear
# PJD 13 Feb 2019   - Update to latest data - oceanonly
#                   - TODO:
#                   - Ishii v7.2 data has problems with grib2
#                   - Fix issue with ORAS4 time dimension, not CF/CDMS cannot read
#                   - Add in time dimension for SODA224
#
#@author: durack1
#"""

#/work/cmip5/historical/ocn/mo/thetao/cmip5.ACCESS1-3.historical.r1i1p1.mo.ocn.Omon.thetao.ver-1.latestX.xml
if [ "$HOST" == "oceanonly.llnl.gov" ]; then
    echo ocean
    cd /work/durack1/Shared/190213_data_density
    export PATH=/export/durack1/anaconda2/envs/cdat80py2/bin:$PATH
    # Obs
    # EN4
    \rm -f obs.EN4.historical.r0i0p0.mo.ocn.Omon.*.xml ; # Cleanup existing *.xmls
    cdscan --include so -x obs.EN4.historical.r0i0p0.mo.ocn.Omon.so.ver-1.latestX.xml /work/durack1/Shared/obs_data/EN4/180723*.nc
    cdscan --include thetao -x obs.EN4.historical.r0i0p0.mo.ocn.Omon.thetao.ver-1.latestX.xml /work/durack1/Shared/obs_data/EN4/180723*.nc
    # IPRC
    \rm -f obs.IPRC.historical.r0i0p0.mo.ocn.Omon.*.xml ; # Cleanup existing *.xmls
    cdscan --include so -a alias_IPRC.txt -x obs.IPRC.historical.r0i0p0.mo.ocn.Omon.so.ver-1.latestX.xml /work/durack1/Shared/obs_data/Argo/IPRC/180719/*.nc
    cdscan --include thetao -a alias_IPRC.txt -x obs.IPRC.historical.r0i0p0.mo.ocn.Omon.thetao.ver-1.latestX.xml /work/durack1/Shared/obs_data/Argo/IPRC/180719/*.nc
    # Ishii
#    \rm -f obs.Ishii.historical.r0i0p0.mo.ocn.Omon.*.xml ; # Cleanup existing *.xmls
#    cdscan --include so -a alias_Ishii.txt -x obs.Ishii.historical.r0i0p0.mo.ocn.Omon.so.ver-1.latestX.xml /work/durack1/Shared/obs_data/Ishii/130702*.nc
#    cdscan --include thetao -a alias_Ishii.txt -x obs.Ishii.historical.r0i0p0.mo.ocn.Omon.thetao.ver-1.latestX.xml /work/durack1/Shared/obs_data/Ishii/130702*.nc
    # JAMSTEC
    \rm -f obs.JAMSTEC.historical.r0i0p0.mo.ocn.Omon.*.xml ; # Cleanup existing *.xmls
    cdscan --include so -a alias_JAMSTEC.txt -x obs.JAMSTEC.historical.r0i0p0.mo.ocn.Omon.so.ver-1.latestX.xml /work/durack1/Shared/obs_data/Argo/JAMSTEC/180910*.nc
    cdscan --include thetao -a alias_JAMSTEC.txt -x obs.JAMSTEC.historical.r0i0p0.mo.ocn.Omon.thetao.ver-1.latestX.xml /work/durack1/Shared/obs_data/Argo/JAMSTEC/180910*.nc
    # SmithAndMurphy2007
    \rm -f obs.SmithAndMurphy2007.historical.r0i0p0.mo.ocn.Omon.*.xml ; # Cleanup existing *.xmls
    cdscan --include so -x obs.SmithAndMurphy2007.historical.r0i0p0.mo.ocn.Omon.so.ver-1.latestX.xml /work/durack1/Shared/obs_data/SmithAndMurphy2007/131011*.nc
    cdscan --include thetao -x obs.SmithAndMurphy2007.historical.r0i0p0.mo.ocn.Omon.thetao.ver-1.latestX.xml /work/durack1/Shared/obs_data/SmithAndMurphy2007/131011*.nc
    # UCSD
    \rm -f obs.UCSD.historical.r0i0p0.mo.ocn.Omon.*.xml ; # Cleanup existing *.xmls
    cdscan --include so -a alias_UCSD.txt -x obs.UCSD.historical.r0i0p0.mo.ocn.Omon.so.ver-1.latestX.xml /work/durack1/Shared/obs_data/Argo/UCSD/180719*.nc
    cdscan --include thetao -a alias_UCSD.txt -x obs.UCSD.historical.r0i0p0.mo.ocn.Omon.thetao.ver-1.latestX.xml /work/durack1/Shared/obs_data/Argo/UCSD/180719*.nc
    # Reanalysis
elif [ "$HOST" == 'crunchy.llnl.gov' ] ; then
    echo crunch
    cd /work/durack1/Shared/160421_data_density
    export PATH=/usr/local/uvcdat/2.4.1/bin:$PATH
    # ORAS4
    # '-t "t"' or '-e 't.id = "time"' options do not work
    \rm -f ocerean.ORAS4.historical.r0i0p0.mo.ocn.Omon.*.xml ; # Cleanup existing *.xmls
    #cdscan  -a alias_ORAS4.txt -e 'deptht.positive = "down"' --include so --time-linear '0,1,months since 1958,gregorian' -x ocerean.ORAS4.historical.r0i0p0.mo.ocn.Omon.so.ver-1.latestX.xml /work/guilyardi/database/ORAS4/so*.nc
    #cdscan -a alias_ORAS4.txt -e 'deptht.positive = "down"' --include so -x ocerean.ORAS4.historical.r0i0p0.mo.ocn.Omon.so.ver-1.latestX.xml /work/guilyardi/database/ORAS4/so*.nc
    #cdscan -a alias_ORAS4.txt -e 'deptht.positive = "down"' --include thetao -x ocerean.ORAS4.historical.r0i0p0.mo.ocn.Omon.thetao.ver-1.latestX.xml /work/guilyardi/database/ORAS4/thetao*.nc
    # SODA224
    \rm -f ocerean.SODA224.historical.r0i0p0.mo.ocn.Omon.*.xml ; # Cleanup existing *.xmls
    cdscan -a alias_SODA224.txt --include so --time-linear '0,1,months since 1870,gregorian' -x ocerean.SODA224.historical.r0i0p0.mo.ocn.Omon.so.ver-1.latestX.xml /work/guilyardi/database/SODA224/*.cdf
    cdscan -a alias_SODA224.txt --include thetao --time-linear '0,1,months since 1870,gregorian' -x ocerean.SODA224.historical.r0i0p0.mo.ocn.Omon.thetao.ver-1.latestX.xml /work/guilyardi/database/SODA224/*.cdf
else
    echo "Platform unrecognized"
fi
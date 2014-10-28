#!/bin/env python
# -*- coding: utf-8 -*-
#
# <------------------------------------------------>
#       Test driver for surface transformation
# <------------------------------------------------>
from surface_transf import surfTransf
#
#
# MPI
#
file_fx = '/work/cmip5/fx/fx/areacello/cmip5.MPI-ESM-LR.historical.r0i0p0.fx.ocn.fx.areacello.ver-v20111006.latestX.xml'
file_tos = '/work/cmip5/historical/ocn/mo/tos/cmip5.MPI-ESM-LR.historical.r1i1p1.mo.ocn.Omon.tos.ver-1.latestX.xml'
file_sos = '/work/cmip5/historical/ocn/mo/sos/cmip5.MPI-ESM-LR.historical.r1i1p1.mo.ocn.Omon.sos.ver-1.latestX.xml'
file_hef = '/work/cmip5/historical/ocn/mo/hfds/cmip5.MPI-ESM-LR.historical.r1i1p1.mo.ocn.Omon.hfds.ver-1.latestX.xml'
file_wfo = '/work/cmip5/historical/ocn/mo/wfo/cmip5.MPI-ESM-LR.historical.r1i1p1.mo.ocn.Omon.wfo.ver-1.latestX.xml'
outfileSurfDen = 'test/cmip5.MPI-ESM-LR.historical.r1i1p1.mo.ocn.Omon.surfden.ver-1.latest.nc'

#
# -----------------------------------------
#  Compute density flux and transformation
# -----------------------------------------
#
surfTransf(file_fx, file_tos, file_sos, file_hef, file_wfo, outfileSurfDen, debug=True,timeint='all')
#surfTransf(file_fx, file_tos, file_sos, file_hef, file_wfo, outfileSurfDen, debug=True,timeint='1,3')


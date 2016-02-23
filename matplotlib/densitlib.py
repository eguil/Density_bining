"""
Densitlib for matplotlib for density plots
(c) Eric Guilyardi
Feb 2016

"""

import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
#
# Build zonal mean with zoom
#

def zon_2dom(plt,ax0,ax1,lat,lev,varBasin,varSigma,unit,minmax,clevsm,cmap,domrho,title,agreelev,noax):

#-- variables
    var   = varBasin['diffBowl']
    varm  = varBasin['meanBowl']
    varag = varBasin['agree']

#-- title
    title = title+' '+varBasin['name']
#-- contour levels
    rhomin = domrho[0]
    rhomid = domrho[1]
    rhomax = domrho[2]

#-- Latmin/max
    latmin = -80.
    latmax = 80.
    deltalat = 20
    labels = ['','60S','40S','20S','Eq','20N','40N','60N','']
#
# ====   Window 1  ===================================================
#
    ax0.axis([latmin,latmax,rhomin,rhomid])
    ax0.invert_yaxis()
    ax0.xaxis.set_ticks(np.arange(latmin, latmax, deltalat))
    ax0.set_xticklabels(labels)
    if noax == 'T':
        ax0.set_yticklabels([])
        ax1.set_yticklabels([])
    if noax == 'R': 
        ax0.yaxis.tick_right()
        ax1.yaxis.tick_right()

#-- levels for diff plot
    levels = MaxNLocator(nbins=minmax[2]).tick_values(minmax[0], minmax[1])

#-- draw filled contours of period diff
    cnplot = ax0.contourf(lat,lev,var,cmap=cmap, levels=levels)

#-- draw agreement contour > agreement level (agreelev)
    cmapbl = LinearSegmentedColormap('cmapbl',bluecol())
    chplot = ax0.contourf(lat,lev,varag,levels=[-agreelev,agreelev],hatches=['..'],colors='none')
    cpplot = ax0.contour(lat,lev,varag,[agreelev-.0001,agreelev+0.00001],cmap=cmapbl, linewidths=2)
    cpplot = ax0.contour(lat,lev,varag,[-agreelev-.0001,-agreelev+0.00001],cmap=cmapbl, linewidths=2)

#-- draw mean contours
    cmapb = LinearSegmentedColormap('cmapb',blkcol())
    cpplot = ax0.contour(lat,lev,varm,clevsm,cmap=cmapb)
    ax0.clabel(cpplot, inline=1, fontsize=10, fmt='%.1f')

#-- draw ptopsigma for 2 periods (yr1 = ref, yr2 = end of serie)
    lnplot1 = ax0.plot(lat,varSigma['yr1'],linestyle='--', color='black', linewidth=2)
    lnplot2 = ax0.plot(lat,varSigma['yr2'],linestyle='-', color='black', linewidth=2)

# 
# ====   Window 2  ===================================================
#
    ax1.axis([latmin,latmax,rhomid,rhomax])
    ax1.invert_yaxis()
    ax1.xaxis.set_ticks(np.arange(latmin, latmax, deltalat))
    ax1.set_xticklabels(labels)

#-- draw filled contours
    cnplot = ax1.contourf(lat,lev,var,cmap=cmap, levels=levels)

#-- draw agreement contour > agreement level (agreelev)
    cmapbl = LinearSegmentedColormap('cmapbl',bluecol())
    chplot = ax1.contourf(lat,lev,varag,levels=[-agreelev,agreelev],hatches=['..'],colors='none')
    cpplot = ax1.contour(lat,lev,varag,[agreelev-.0001,agreelev+0.00001],cmap=cmapbl,linewidths=2)
    cpplot = ax1.contour(lat,lev,varag,[-agreelev-.0001,-agreelev+0.00001],cmap=cmapbl,linewidths=2)

#-- draw mean contours
    cmapb = LinearSegmentedColormap('cmapb',blkcol())
    cpplot = ax1.contour(lat,lev,varm,clevsm,cmap=cmapb)
    ax1.clabel(cpplot, inline=1, fontsize=10, fmt='%.1f')

#-- draw ptopsigma for 2 periods (yr1 = ref, yr2 = end of serie)
    lnplot1b = ax1.plot(lat,varSigma['yr1'],linestyle='--', color='black', linewidth=2, label='<1950')
    lnplot2b = ax1.plot(lat,varSigma['yr2'],linestyle='-', color='black', linewidth=2, label='2000')

#-- Add legend for bowl position
    plt.legend(loc='upper right', title='Bowl')

#-- add plot title
    ax0.set_title(title)

    return [cnplot,lnplot1,lnplot2]

#
# Convert GMT color palette form post-it to pylab
#
# http://scipy.github.io/old-wiki/pages/Cookbook/Matplotlib/Loading_a_colormap_dynamically.html
#
####################
def gmtColormap(fileName,GMTPath = None):
      """
      gmtColormap(fileName,GMTPath='~/python/cmaps/')
      
      Returns a dict for use w/ LinearSegmentedColormap
      cdict = gmtcolormapPylab.gmtcolormapPylab('spectrum-light')
      colormap = pylab.cm.colors.LinearSegmentedColormap('spectrum-light',cdict)
      """
      import colorsys
      import numpy as N

      if type(GMTPath) == type(None):
          filePath = "/usr/local/cmaps/"+ fileName+".cpt"
      else:
          filePath = GMTPath+"/"+ fileName +".cpt"
      try:
          f = open(filePath)
      except:
          print "file ",filePath, "not found"
          return None
 
      lines = f.readlines()
      f.close()

      x = N.array([])
      r = N.array([])
      g = N.array([])
      b = N.array([])
      colorModel = "RGB"
      for l in lines:
          ls = l.split()
          if l[0] == "#":
             if ls[-1] == "HSV":
                 colorModel = "HSV"
                 continue
             else:
                 continue
          if ls[0] == "B" or ls[0] == "F" or ls[0] == "N":
             pass
          else:
              x=N.append(x,float(ls[0]))
              r=N.append(r,float(ls[1]))
              g=N.append(g,float(ls[2]))
              b=N.append(b,float(ls[3]))
              xtemp = float(ls[4])
              rtemp = float(ls[5])
              gtemp = float(ls[6])
              btemp = float(ls[7])
              

      x=N.append(x,xtemp)
      r=N.append(r,rtemp)
      g=N.append(g,gtemp)
      b=N.append(b,btemp)

      nTable = len(r)
      if colorModel == "HSV":
         for i in range(r.shape[0]):
             rr,gg,bb = colorsys.hsv_to_rgb(r[i]/360.,g[i],b[i])
             r[i] = rr ; g[i] = gg ; b[i] = bb
      if colorModel == "HSV":
         for i in range(r.shape[0]):
             rr,gg,bb = colorsys.hsv_to_rgb(r[i]/360.,g[i],b[i])
             r[i] = rr ; g[i] = gg ; b[i] = bb
      if colorModel == "RGB":
          r = r/255.
          g = g/255.
          b = b/255.
      #print shape(x)

      xNorm = (x - x[0])/(x[-1] - x[0])

      red = []
      blue = []
      green = []
      for i in range(len(x)):
          red.append([xNorm[i],r[i],r[i]])
          green.append([xNorm[i],g[i],g[i]])
          blue.append([xNorm[i],b[i],b[i]])
      colorDict = {"red":red, "green":green, "blue":blue}
      return (colorDict)    
   
#- def black color bar
def blkcol(): 

      cdict = {'red':   ((0.0,  0.0, 0.0),
                         (1.0,  1.0, 1.0)),
               
               'green': ((0.0,  0.0, 0.0),
                         (1.0,  1.0, 1.0)),
               
               'blue':  ((0.0,  0.0, 0.0),
                         (1.0,  1.0, 1.0))}

      return (cdict)
#- def blues color bar
def bluecol(): 

      cdict = {'red':   ((0.0,  0.0, 0.0),
                         (1.0,  0.0, 0.0)),
               
               'green': ((0.0,  0.0, 0.0),
                         (1.0,  0.0, 0.0)),
               
               'blue':  ((0.0,  0.0, 1.0),
                         (1.0,  1.0, 1.0))}


      return (cdict)


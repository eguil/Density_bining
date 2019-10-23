'''
Compute, in zonal means, how much of the ocean (per basin) has emerged during from 1860 to 2100
'''

import os
import glob
from netCDF4 import Dataset as open_ncfile
import matplotlib.pyplot as plt
import numpy as np
import datetime
import pickle

# -- Read result
emerge = pickle.load( open( "/home/ysilvy/Density_bining/Yona_analysis/data/percentage_emergence_medians_meanhistNat.pkl", "rb" ) )

# -- Median and range
median_emerge = np.ma.median(emerge,axis=1)
pc25_emerge = np.percentile(emerge,25,axis=1)
pc75_emerge = np.percentile(emerge,75,axis=1)
time = np.arange(1860,2100)

# -- Plot

fig, axes = plt.subplots(1,3,sharex=True,sharey=True,figsize=(16, 5))
axes[0].plot(time,median_emerge[:,1],color='k')
axes[0].fill_between(time,pc25_emerge[:,1],pc75_emerge[:,1],color='0.8') #alpha=0.3
axes[0].set_xlim([1920,2080])
axes[0].set_ylim([0,83])
axes[0].axvline(x=2005,ls='--',color='k',lw=0.5)
axes[0].grid(axis='y')
#axes[0].set_title('Atlantic',fontweight='bold')
axes[0].text(2050,5,'Atlantic',fontweight='bold',fontsize=15,va='center',ha='center')
axes[1].plot(time,median_emerge[:,2],color='k')
axes[1].fill_between(time,pc25_emerge[:,2],pc75_emerge[:,2],color='0.8')
axes[1].axvline(x=2005,ls='--',color='k',lw=0.5)
#axes[1].set_title('Pacific',fontweight='bold')
axes[1].text(2050,5,'Pacific',fontweight='bold',fontsize=15,va='center',ha='center')
axes[1].grid(axis='y')
axes[2].plot(time,median_emerge[:,3],color='k')
axes[2].fill_between(time,pc25_emerge[:,3],pc75_emerge[:,3],color='0.8')
axes[2].axvline(x=2005,ls='--',color='k',lw=0.5)
#axes[2].set_title('Indian',fontweight='bold')
axes[2].text(2050,5,'Indian',fontweight='bold',fontsize=15,va='center',ha='center')
axes[2].grid(axis='y')
axes[0].set_ylabel('% of basin zonal mean',fontweight='bold',fontsize=14)
axes[0].set_xticks(np.arange(1920,2081,20))
axes[1].tick_params(axis='y', labelleft='on')
axes[2].tick_params(axis='y', labelleft='on')
plt.subplots_adjust(wspace=0.1,top=0.85,left=0.04, right=0.92)
plt.suptitle('Percentage of basin emergence in zonal means under the bowl',fontweight='bold', fontsize=13)
plt.figtext(.006,.95,'b',fontweight='bold',fontsize=16)

for i in range(3):
    plt.setp(axes[i].get_xticklabels(), fontweight='bold')
    plt.setp(axes[i].get_yticklabels(), fontweight='bold')
    axes[i].xaxis.set_tick_params(which='major',width=2)
    axes[i].yaxis.set_tick_params(which='major',width=2)

# Date
now = datetime.datetime.now()
date = now.strftime("%Y-%m-%d")

plotName = 'percentage_basin_emergence_medians_rcp85_meanhistNat_paper'
figureDir = 'models/ToE/'
plt.savefig('/home/ysilvy/figures/'+figureDir+plotName+'.png') #,bbox_inches='tight')
#plt.show()

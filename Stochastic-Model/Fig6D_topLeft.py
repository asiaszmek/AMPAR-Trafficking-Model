# -*- coding: utf-8 -*-
"""
Created on Tue Nov  3 10:13:49 2020

@author: Moritz
"""

"""Fig 6D

This script reproduces the plots seen in Fig 6D of "The biophysical basis underlying the maintenance of early phase long-term potentiation".

This script requires that `read_roi`,`numpy`,`matplotlib` and `seaborn` are installed within the Python environment you are running this script in.
"""

import sys
sys.path.append('../')

import numpy as np
import matplotlib.pyplot as plt
from ampartrafficking.frap import FRAP
import seaborn as sns
sns.set_style("ticks")
from matplotlib.font_manager import FontProperties
fontLgd = FontProperties()
fontLgd.set_size('small')
from read_roi import read_roi_zip

#%%


roi = read_roi_zip('Data\\Tatavarty2013_FRAP_Fig3D_RoiSet.zip')

xAxis=80/(roi['xAxis']['x2']-roi['xAxis']['x1'])
yAxis=100/(roi['yAxis']['y2']-roi['yAxis']['y1'])
dataPoints_x=((np.array(roi['dataPoints']['x'])-roi['yAxis']['x1'])*xAxis)-1.5
dataPoints_y=((np.array(roi['dataPoints']['y'])-roi['yAxis']['y1'])*yAxis)


roi_M = read_roi_zip('Data\\Makino2009_FRAP-GluR1_Fig1D_RoiSet.zip')

xAxis=30/(roi_M['xAxis']['x2']-roi_M['xAxis']['x1'])
yAxis=120/(roi_M['yAxis']['y2']-roi_M['yAxis']['y1'])
dataPoints_x_M=((np.array(roi_M['dataPoints']['x'])-roi_M['yAxis']['x1'])*xAxis)+0.1
dataPoints_y_M=((np.array(roi_M['dataPoints']['y'])-roi_M['yAxis']['y1'])*yAxis)-3.5  
        
#%%
    
#def main():

SaveFig=0

#Cooperative binding; FRAP dependence on PSD size P:

duration=15000#Duration in s

Nr_Trials=100
N_List=[5,9,12]#[6,10,14]#
UFP_List=[10]

beta=1
alpha=16
kUB=0.0005
kBU=0.1

A_spine_basal=0.898
    
B_N, U_N, B_notBleached_N, U_notBleached_N, PSD, Time=FRAP(N_List, UFP_List, beta, alpha, kUB, kBU, duration, Nr_Trials)


#%%

col=sns.color_palette('colorblind')[3::]

fig=plt.figure(figsize=(3.5,2), dpi=150)

for i,N in enumerate(N_List):

    A_spine=N**2/70*A_spine_basal
    
    for j,UFP_0 in enumerate(UFP_List):
        
        BMean=np.mean(np.mean((B_N[i][j]+B_notBleached_N[i][j]), axis=0)[int(len(Time)/2)::])
        UMean=np.mean(np.mean((U_N[i][j]+U_notBleached_N[i][j]), axis=0)[int(len(Time)/2)::])
        print(BMean)
        print(UMean)
        
        Y=(np.mean(B_notBleached_N[i][j], axis=0)+np.mean(U_notBleached_N[i][j], axis=0))/(BMean+UMean)*100
        
        plt.plot(Time-duration/2/60,Y, color=col[i], linewidth=2, label=r'P={0:1.0f}, $\langle B \rangle$={1:1.0f}'.format(N**2, BMean))
        
plt.text(0.3,0.55,r'$\langle U/A_{spine} \rangle$='+'{0:1.1f} $\#/ \mu m^2$'.format(UMean/A_spine), transform=fig.axes[0].transAxes)
        
#sns.lineplot(dataPoints_x,dataPoints_y, color='k', marker='s', label='Tatavarty 2013')
sns.lineplot(dataPoints_x_M,dataPoints_y_M, color='k', alpha=0.6, linewidth=2, marker='d', label='Makino 2009')

plt.axhline(100, color='k', zorder=0)
plt.xlim(0,35)
plt.ylim(0,120)
plt.xlabel('Time (min)')
plt.ylabel('AMPAR \n Recovery (%)')
plt.legend(ncol=2, prop=fontLgd)
sns.despine()
fig.tight_layout()
if SaveFig==1:
    print('Save')
    fig.savefig('Figures\\Fig6D_topLeft.png', bbox_inches="tight", dpi=400)
    fig.savefig('Figures\\Fig6D_topLeft.svg', bbox_inches="tight", dpi=400)


#%%


fig=plt.figure(figsize=(3.5,2), dpi=150)

for i,N in enumerate(N_List):

    A_spine=N**2/70*A_spine_basal
    
    for j,UFP_0 in enumerate(UFP_List):
        
        BMean=np.mean(np.mean((B_N[i][j]+B_notBleached_N[i][j]), axis=0)[int(len(Time)/2)::])
        UMean=np.mean(np.mean((U_N[i][j]+U_notBleached_N[i][j]), axis=0)[int(len(Time)/2)::])
        print(BMean)
        print(UMean)
        
        Y=(np.mean(B_notBleached_N[i][j], axis=0)+np.mean(U_notBleached_N[i][j], axis=0))/(BMean+UMean)*100
        Y2=(np.mean(B_N[i][j], axis=0)+np.mean(U_N[i][j], axis=0))/(BMean+UMean)*100
        
        plt.plot(Time-duration/2/60,Y+Y2, color=col[i], linewidth=2, label=r'P={0:1.0f}, $\langle B \rangle$={1:1.0f}'.format(N**2, BMean))
        
plt.text(0.3,0.55,r'$\langle U/A_{spine} \rangle$='+'{0:1.1f} $\#/ \mu m^2$'.format(UMean/A_spine), transform=fig.axes[0].transAxes)
        

plt.axhline(100, color='k', zorder=0)
# plt.xlim(0,35)
plt.ylim(0,120)
plt.xlabel('Time (min)')
plt.ylabel('AMPAR \n bleached+not bleached (%)')
plt.legend(ncol=2, prop=fontLgd)
sns.despine()
fig.tight_layout()

#%%

# if __name__ == "__main__":
#     main()

# -*- coding: utf-8 -*-
"""
Created on Wed Nov  4 21:42:22 2020

@author: Moritz
"""

"""Fig 6D

This script reproduces the plots seen in Fig 6D of "The biophysical basis underlying the maintenance of early phase long-term potentiation".

This script requires that `numpy`,`matplotlib` and `seaborn` are installed within the Python environment you are running this script in.
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

#%%

#No Cooperative binding; FRAP dependence on mobile receptor concentartion U/Aspine:

SaveFig=0#1

duration=15000#Duration in s
    
Nr_Trials=100
N_List=[12]#[9]
UFP_List=[10,30,60]

beta=0
alpha=0
kUB=0.0036
kBU=0.1

A_spine_basal=0.898
    
B_N, U_N, B_notBleached_N, U_notBleached_N, PSD, Time=FRAP(N_List, UFP_List, beta, alpha, kUB, kBU, duration, Nr_Trials)


plt.figure()
plt.imshow(PSD)
plt.colorbar()

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
        
        if j==0 and i==0:
            plt.plot(Time-duration/2/60,Y, color=col[j], linewidth=2, label=r'$\langle B \rangle={:1.0f}$, '.format(BMean)+r'$\langle U/A_{spine} \rangle$='+r'{0:1.1f} $\#/ \mu m^2$'.format(UMean/A_spine))
        else:
            plt.plot(Time-duration/2/60,Y, color=col[j], linewidth=2, label=r'$={:1.0f}$, '.format(BMean)+r'$=$'+r'{0:1.1f} $\#/ \mu m^2$'.format(UMean/A_spine))
            
        
plt.text(0.3,0.55,r'P={0:1.0f}'.format(N**2), transform=fig.axes[0].transAxes)

plt.axhline(100, color='k', zorder=0)
plt.xlim(0,35)
plt.ylim(0,120)
plt.xlabel('Time (min)')
plt.ylabel('AMPAR \n Recovery (%)')
plt.legend(ncol=1, prop=fontLgd)
sns.despine()
fig.tight_layout()
if SaveFig==1:
    print('Save')
    fig.savefig('Figures\\Fig6D_bottomRight.png', bbox_inches="tight", dpi=400)
    fig.savefig('Figures\\Fig6D_bottomRight.svg', bbox_inches="tight", dpi=400)
    
    
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
        
        if j==0 and i==0:
            plt.plot(Time-duration/2/60,Y+Y2, color=col[j], linewidth=2, label=r'$\langle B \rangle={:1.0f}$, '.format(BMean)+r'$\langle U/A_{spine} \rangle$='+r'{0:1.1f} $\#/ \mu m^2$'.format(UMean/A_spine))
        else:
            plt.plot(Time-duration/2/60,Y+Y2, color=col[j], linewidth=2, label=r'$={:1.0f}$, '.format(BMean)+r'$=$'+r'{0:1.1f} $\#/ \mu m^2$'.format(UMean/A_spine))
        
plt.text(0.3,0.55,r'P={0:1.0f}'.format(N**2), transform=fig.axes[0].transAxes)
        

plt.axhline(100, color='k', zorder=0)
# plt.xlim(0,35)
plt.ylim(0,120)
plt.xlabel('Time (min)')
plt.ylabel('AMPAR \n bleached+not bleached (%)')
plt.legend(ncol=1, prop=fontLgd)
sns.despine()
fig.tight_layout()
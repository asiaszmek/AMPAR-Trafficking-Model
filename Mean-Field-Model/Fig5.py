# -*- coding: utf-8 -*-
"""
Created on Tue Nov  3 17:01:58 2020

@author: Moritz
"""

"""Fig 5

This script reproduces the plots seen in Fig 5 of "The biophysical basis underlying the maintenance of early phase long-term potentiation".

This script requires that `numpy`,`scipy.integrate`,`matplotlib` and `seaborn` are installed within the Python environment you are running this script in.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import mf as mf
import seaborn as sns
sns.set_style("ticks")
plt.rcParams['svg.fonttype'] = 'none'
from matplotlib.font_manager import FontProperties
fontLgd = FontProperties()
fontLgd.set_size('xx-small')

#%%

# def main():
    
SaveFig=0

#Import experimental Data for comparison with the model. Taken from Penn et al 2017 and Barco et al. 2002 (see Paper for reference details):
Barco2002_x=np.loadtxt("Data\\Barco2002_x.csv", delimiter=",")
Barco2002_y=np.loadtxt("Data\\Barco2002_y.csv", delimiter=",")

#%%

#Color palettes used in plots:
col1=sns.dark_palette(sns.color_palette("colorblind", 10)[0],4)[1:4]

#%%

#Set parameter values and initial conditions:
    
kin=0.2
kout=0.018
kexo0=0.0018
kendo=0.0021
kBU=0.1
kUB0=0.0005
kin_RE=0.1
kout_RE=0.00769
Vspine0=0.08

kUB=mf.Parameter(kUB0)
kUB.timecourse(mf.Stim_Resp,[1,5,5,60])
kexo=mf.Parameter(kexo0)
kexo.timecourse(mf.Stim_Resp,[1,5,25,60])
Vspine=mf.Parameter(Vspine0)
Vspine.timecourse(mf.DV,[Vspine0,True])
    
t=np.arange(0,120*60)

#%%

#Fig 5 B:
    
initList = [[10,8.7,13],[10,29.5,13],[10,56.9,13]]
pList=[40,100,200]

sLTP=0
Cooperativity=1
    
fig=plt.figure(figsize=(3.5,2), dpi=150)
sns.lineplot(Barco2002_x,Barco2002_y, color='k', marker="s", ms=5, legend=False, label='E-LTP', zorder=0)
    
for i,P,Init in zip(range(len(pList)),pList,initList):

    
    Model=mf.Model_system()
    
    solve,infodict = odeint(Model.odes,Init,t,args=(Vspine,P,kin,kout,kexo,kendo,kUB,kBU,Cooperativity,sLTP,kin_RE,kout_RE),full_output=True)
    
    
    plt.plot(t/60,solve.T[1]/Init[1]*100, color=col1[i], linewidth=2, label='$P={0:1.1f},\, B^*{1:1.1f}$'.format(P,Init[1]))
    plt.xlabel('Time (min)')
    plt.ylabel('Bound AMPARs $B$ (%)')
    lgd=plt.legend(bbox_to_anchor=(0.15,0.6,0,0), loc="lower left", mode="normal", borderaxespad=1.5, ncol=2, prop=fontLgd)
    plt.axhline(100, color='k', linestyle='--', linewidth=0.5)
    plt.xlim(-3,130)
    plt.ylim(80)
    sns.despine()
    plt.tight_layout()
    if SaveFig==1:
        print('Save')
        fig.savefig('Figures\\Fig5B.png', bbox_inches="tight", dpi=400)
        fig.savefig('Figures\\Fig5B.svg', bbox_inches="tight", dpi=400)
        
#%%

#Fig 5 C:

vList=[0.039,0.08,0.376]
initList = [[6.2,10,13],[10,20,13],[28.0,56.5,13]]

sLTP=1
Cooperativity=1
    
fig=plt.figure(figsize=(3.5,1.25), dpi=150)

for i,Init, Vspine0 in zip(range(len(vList)),initList,vList):

    Vspine=mf.Parameter(Vspine0)
    Vspine.timecourse(mf.DV,[Vspine0,True])
    
    P=4*np.pi*(3*Vspine.base_value/(4*np.pi))**(2/3)/0.89*70
    
    Model=mf.Model_system()
    
    solve,infodict = odeint(Model.odes,Init,t,args=(Vspine,P,kin,kout,kexo,kendo,kUB,kBU,Cooperativity,sLTP,kin_RE,kout_RE),full_output=True)

    
    plt.plot(t/60,solve.T[1]-Init[1], color=col1[i], label='$V_{spine}^{0}=$'+'{:1.3f} $\mu m^3$'.format(Vspine0))
    plt.xlabel('Time (min)')
    plt.ylabel('$\Delta B$ (#)')
    lgd=plt.legend(bbox_to_anchor=(0.57,-0.2,0,0), loc="lower left", mode="normal", borderaxespad=1.5, ncol=1, prop=fontLgd)
    plt.axhline(0, color='k', linestyle='--', linewidth=0.5)
    sns.despine()
    plt.tight_layout()
    if SaveFig==1:
        print('Save')
        fig.savefig('Figures\\Fig5C.png', bbox_inches="tight", dpi=400)
        fig.savefig('Figures\\Fig5C.svg', bbox_inches="tight", dpi=400)    

#%%

#Fig 5 D:

vList=[0.039,0.08,0.376]
initList = [[5.9,8.5,13],[10,20,13],[38.9,142.4,13]]

sLTP=1
Cooperativity=1
    
fig=plt.figure(figsize=(3.5,1.25), dpi=150)

for i,Init, Vspine0 in zip(range(len(vList)),initList,vList):
    
    Vspine=mf.Parameter(Vspine0)
    Vspine.timecourse(mf.DV,[Vspine0,True])
    
    #Calculate the exocytosis event size Sexo from the spine volume:
    Sexo_0=13*Vspine0/0.08
    kin_RE=0.1
    kout_RE=kin_RE/Sexo_0
    
    Init[2]=Sexo_0
    
    P=4*np.pi*(3*Vspine.base_value/(4*np.pi))**(2/3)/0.89*70
    
    Model=mf.Model_system()
    
    solve,infodict = odeint(Model.odes,Init,t,args=(Vspine,P,kin,kout,kexo,kendo,kUB,kBU,Cooperativity,sLTP,kin_RE,kout_RE),full_output=True)
   
    
    plt.plot(t/60,solve.T[1]-Init[1], color=col1[i], label='$V_{spine}^{0}=$'+'{:1.3f} $\mu m^3$'.format(Vspine0))
    plt.xlabel('Time (min)')
    plt.ylabel('$\Delta B$ (#)')
    lgd=plt.legend(bbox_to_anchor=(0.57,-0.2,0,0), loc="lower left", mode="normal", borderaxespad=1.5, ncol=1, prop=fontLgd)
    plt.axhline(0, color='k', linestyle='--', linewidth=0.5)
    sns.despine()
    plt.tight_layout()
    if SaveFig==1:
        print('Save')
        fig.savefig('Figures\\Fig5D.png', bbox_inches="tight", dpi=400)
        fig.savefig('Figures\\Fig5D.svg', bbox_inches="tight", dpi=400)


# if __name__ == "__main__":
#     main()
    
    
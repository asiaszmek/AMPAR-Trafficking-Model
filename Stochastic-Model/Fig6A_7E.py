# -*- coding: utf-8 -*-
"""
Created on Thu Nov 12 08:18:23 2020

@author: Moritz
"""

"""Fig 6A and 7E

This script reproduces the plots seen in Fig 6A and 7E of "The biophysical basis underlying the maintenance of early phase long-term potentiation".

This script requires that `numpy`,`scipy.integrate`,`matplotlib` and `seaborn` are installed within the Python environment you are running this script in.
"""

import sys
sys.path.append('../')

import numpy as np
import matplotlib.pyplot as plt
from ampartrafficking import stochastic_model as sm
from ampartrafficking import rate_model as rm
import seaborn as sns
sns.set_style("ticks")
from matplotlib.font_manager import FontProperties 
fontLgd = FontProperties()
fontLgd.set_size('small')
plt.rcParams['svg.fonttype'] = 'none'

#%%

SaveFig=0

A_spine=0.898
    
duration=100000#Duration in s
kout=0.01796
kin=0.02
kexo0=0.0018
kendo=0.002058
S_exo=13

dt=0.5#the size of each time step

Nr_Trials=1

beta=1
alpha=16
kUB0=0.0005048
kBU=0.1


UFP_List=np.array([5,7,8.5,10,12,15,20,25,30])
N_List=[3,5,8,10,12,15]


#%%

B_N=[]
U_N=[]

ID_basal=1

for N in N_List:
    
    B_U=[]
    U_U=[]

    for UFP_0 in UFP_List:
        
        dt=sm.calcTimeStep(UFP_0,A_spine,kUB0,alpha,kBU,kout+kendo, kin+kexo0*S_exo/UFP_0)
        
        B_Tr=[]
        U_Tr=[]

        
        for Trial in range(0,Nr_Trials):
            
            if Trial%10==0:
                print('Trial:',Trial)
                print('P:',N**2, 'U:', UFP_0)
                
            
            PSD=np.zeros((N,N))

            Time=[]
            B_t=[]
            U_t=[]

            U=UFP_0

            for t in np.arange(0,duration+dt,dt):
                
                
                NN=sm.nearestNeighbours(PSD)
                
                Mbu=sm.kBUcoop(kBU, NN, PSD, ID_basal)*dt
                Mub=sm.kUBcoop(kUB0*U/A_spine, NN, PSD)*dt

                PSD,dBoff,dBon=sm.probabilityEval(Mub,Mbu,PSD,ID_basal)
                
                pout=(kout*U/A_spine+kendo*U/A_spine)*dt
                pin=(kin*UFP_0+kexo0*S_exo)*dt
                U=sm.update_mobilePool(U,pin,pout,dBoff,dBon)
                
                if t%0.5==0:
                    Time.append(t)
                    B_t.append(np.sum(PSD==1))
                    U_t.append(U)
                
            B_Tr.append(B_t)
            U_Tr.append(U_t)

        B_U.append(B_Tr)
        U_U.append(U_Tr)

    B_N.append(B_U)
    U_N.append(U_U)


#%%

plt.figure()
plt.imshow(PSD)
plt.colorbar()

#%%

BMean_N=[]
BStd_N=[]
UMean_N=[]
UStd_N=[]
for i,N in enumerate(N_List):
    BMean_U=[]
    BStd_U=[]
    UMean_U=[]
    UStd_U=[]
    for j,UFP_0 in enumerate(UFP_List):
        
        BMean=0
        BStd=0
        UMean=0
        UStd=0
        for k in range(0,Nr_Trials):
            BMean+=np.mean(B_N[i][j][k][int(len(Time)/2)::])
            BStd+=np.std(B_N[i][j][k][int(len(Time)/2)::])
            UMean+=np.mean(U_N[i][j][k][int(len(Time)/2)::])
            UStd+=np.std(U_N[i][j][k][int(len(Time)/2)::])
            
        BMean/=Nr_Trials
        BStd/=Nr_Trials
        UMean/=Nr_Trials
        UStd/=Nr_Trials
        
        print(BMean)
        print(BStd)
        print(UMean)
        print(UStd)
        
        BStd_U.append(BStd)
        BMean_U.append(BMean)
        UStd_U.append(UStd)
        UMean_U.append(UMean)

    BStd_N.append(BStd_U)
    BMean_N.append(BMean_U)
    UStd_N.append(UStd_U)
    UMean_N.append(UMean_U)
    
BStd_N=np.array(BStd_N)
BMean_N=np.array(BMean_N)
UStd_N=np.array(UStd_N)
UMean_N=np.array(UMean_N)
    
#%%

Cooperativity=1

BFP_N=[]
UFP_N=[]

for N in N_List:
    
    P=N**2
    BFP_U=[]
    UFP_U=[]

    for UFP_0 in np.linspace(0,max(UFP_List),40):
    
        UFP=rm.UFP_(kexo0,S_exo,kendo,kin*UFP_0,kout,A_spine)
        BFP=rm.BFP_(UFP,kUB0,kBU,P,A_spine,Cooperativity)
        
        BFP_U.append(BFP)
        UFP_U.append(UFP)

    BFP_N.append(BFP_U)
    UFP_N.append(UFP_U)

BFP_N=np.array(BFP_N)
UFP_N=np.array(UFP_N)

#%%

col=sns.color_palette('colorblind')

fig=plt.figure(figsize=(3.5,2), dpi=150)
for i in np.arange(1,len(N_List),2):
    plt.errorbar(x=UMean_N[i,:]/A_spine, y=BMean_N[i,:], xerr=UStd_N[i,:]/2, yerr=BStd_N[i,:]/2, color=col[i], linewidth=2, label='P={:1.0f}'.format(N_List[i]**2),zorder=0)        
    # plt.errorbar(x=np.array(UFP_List)/A_spine, y=BMean_N[i,:], yerr=BStd_N[i,:]/2, color=col[i], linewidth=2, label='P={:1.0f}'.format(N_List[i]**2),zorder=0)  
    
    plt.plot(UFP_N[i]/A_spine,BFP_N[i], color=np.array(col[i])*0.75, linestyle='--', linewidth=2, zorder=1)

plt.legend(prop=fontLgd)
plt.xlim(5,25)
plt.xlabel('Mobile AMPAR conc. $U/A_{spine} \, (\#/\mu m^2)$')
plt.ylabel('Bound AMPARs $B^*$ (#)')
sns.despine()
fig.tight_layout()
if SaveFig==1:
    fig.savefig('Figures\\Fig6A.png', bbox_inches="tight", dpi=400)
    fig.savefig('Figures\\Fig6A.svg', bbox_inches="tight", dpi=400)

#%%

plt.figure()
plt.plot(np.arange(0,duration+0.5,0.5)/60,np.array(U_t)/A_spine, label='Trace from stoch. sim.')
plt.axhline(np.mean(np.array(U_t)/A_spine), color='r', linewidth=6, label='Mean Stochastic Model')
plt.axhline(UFP_U[-1]/A_spine, color='g', linewidth=4, label='Mean Rate Model')
plt.axhline((kexo0*S_exo+kin*30)/(kendo+kout), color='k', linewidth=2, label='Mean from fixed point eq.')
plt.legend()
#plt.ylim(0,60)
plt.xlabel('Time (min)')
plt.ylabel('$U/A_{spine}$')

#%%



Cooperativity=1

BFP_N=[]
UFP_N=[]

for N in N_List:
    
    P=N**2
    BFP_U=[]
    UFP_U=[]

    for UFP_0 in [5.0,10.0,15.0,25.0]:
    
        UFP=rm.UFP_(kexo0,S_exo,kendo,kin*UFP_0,kout,A_spine)
        BFP=rm.BFP_(UFP,kUB0,kBU,P,A_spine,Cooperativity)
        
        BFP_U.append(BFP)
        UFP_U.append(UFP)

    BFP_N.append(BFP_U)
    UFP_N.append(UFP_U)

BFP_N=np.array(BFP_N)
UFP_N=np.array(UFP_N)
col=sns.color_palette('colorblind')

#%%


fig=plt.figure(figsize=(3,2), dpi=150)

i=-1
for j in np.arange(0,len(UFP_List)):
    if UFP_List[j] in [5,10,15,25]:
        i+=1
        if i==0:
            plt.errorbar(x=np.array(N_List)**2, y=BMean_N[::,j], yerr=BStd_N[::,j]/2, color=col[i], label=r'$\langle U/A_{spine} \rangle =$'+'{:1.1f}'.format(UMean_N[0][j]/A_spine)+'$\, \#/\mu m^2$',zorder=0)
        else:
            plt.errorbar(x=np.array(N_List)**2, y=BMean_N[::,j], yerr=BStd_N[::,j]/2, color=col[i], label='$=$'+'{:1.1f}'.format(UMean_N[0][j]/A_spine)+'$\, \#/\mu m^2$',zorder=0)
    print(UFP_List[j])


for i in range(4): 
    plt.plot(np.array(N_List)**2,BFP_N.T[i], color=np.array(col[i])*0.75, linestyle='--', linewidth=2, zorder=1)
    
plt.xlim(0,230)
plt.ylim(-10,225)
plt.legend(ncol=1,prop=fontLgd)
plt.xlabel('Slots $P$')
plt.ylabel('Bound AMPARs \n $B^*$ (#)')
sns.despine()
fig.tight_layout()
if SaveFig==1:
    fig.savefig('Figures\\Fig7E.png', bbox_inches="tight", dpi=400)
    fig.savefig('Figures\\Fig7E.svg', bbox_inches="tight", dpi=400)
    
#%%


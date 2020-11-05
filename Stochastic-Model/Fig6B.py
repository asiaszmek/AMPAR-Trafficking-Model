# -*- coding: utf-8 -*-
"""
Created on Tue Nov  3 11:43:19 2020

@author: Moritz
"""

"""Fig 6B

This script reproduces the plots seen in Fig 6B of "The biophysical basis underlying the maintenance of early phase long-term potentiation", which shows the coefficient of variation of the cooperative model for differenr values of bound receptors B and PSD sizes P.

This script requires that `numpy`,`matplotlib` and `seaborn` are installed within the Python environment you are running this script in.
"""

import numpy as np
import matplotlib.pyplot as plt
import sm as sm
import seaborn as sns
sns.set_style("ticks")
from matplotlib.font_manager import FontProperties 
fontLgd = FontProperties()
fontLgd.set_size('small')

#%%

# def main():

SaveFig=0
    
A_spine_basal=0.898
    
duration=100000#Duration in s
kBU=0.1
kout=0.018
kin=0.02
#dt=0.5#the size of each time step

Nr_Trials=1

beta=1#0#
alpha=16#0#
kUB=0.0005#*7
kBU=0.1#/100


N_List=[6,8,10]#[5,9,12]#
UFP_List=[5,10,20]#[10]#

N_List=[5,6,8,10]
UFP_List=np.array([3,5,7,10,12,14,18,30])
        


#%%

B_N=[]
U_N=[]

ID_basal=1

for N in N_List:
    
    A_spine=N**2/70*A_spine_basal
    
    B_U=[]
    U_U=[]

    for UFP_0 in UFP_List:
        
        dt=sm.calcTimeStep(UFP_0,A_spine,kUB,alpha,kBU,kout,kin)
        
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
                Mub=sm.kUBcoop(kUB*U/A_spine, NN, PSD)*dt

                PSD,dBoff,dBon=sm.probabilityEval(Mub,Mbu,PSD,ID_basal)
                
                pout=kout*U/A_spine*dt
                pin=kin*UFP_0*dt
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
BCV_N=[]
for i,N in enumerate(N_List):
    BMean_U=[]
    BCV_U=[]
    for j,UFP_0 in enumerate(UFP_List):
        
        BMean=0
        BStd=0
        for k in range(0,Nr_Trials):
            BMean+=np.mean(B_N[i][j][k][int(len(Time)/5)::])
            BStd+=np.std(B_N[i][j][k][int(len(Time)/5)::])
            
        BMean/=Nr_Trials
        BStd/=Nr_Trials
        
        print(BMean)
        print(BStd)
        print(BStd/BMean*100)
        
        BCV_U.append(BStd/BMean*100)
        BMean_U.append(BMean)

    BCV_N.append(BCV_U)
    BMean_N.append(BMean_U)

#%%


col=sns.color_palette('colorblind')

fig=plt.figure(figsize=(3.5,2), dpi=150)
for i,x,y in (zip(range(len(N_List)),BMean_N,BCV_N)):
    if i>0:
        plt.plot(x,y, marker='o', markersize=2, linewidth=1.5, linestyle='-', color=col[i], label='P={:1.0f}'.format(N_List[i]**2))

plt.xlim(0,70)
plt.ylim(0,100)
plt.xlabel('Bound AMPARs $B$ (#)')
plt.ylabel('Coefficient of variation \n CV (%)')
plt.legend(prop=fontLgd)
sns.despine(top=True)
fig.tight_layout()

if SaveFig==1:
    fig.savefig('Figures\\CV.png', bbox_inches="tight", dpi=400)
    fig.savefig('Figures\\CV.svg', bbox_inches="tight", dpi=400)
    

#%%


# if __name__ == "__main__":
#     main()



#%%  

# B_N=np.array(B_N)
# U_N=np.array(U_N)

# col=sns.color_palette('colorblind')[3::]

# dn=1#25

# from scipy import stats

# for i,N in enumerate(N_List):
    
#     A_spine=N**2/70*A_spine_basal
    
#     for j,UFP_0 in enumerate(UFP_List):
        
#         BMean=0
#         UMean=0
#         for k in range(len(B_N[i][j])):
#             BMean+=np.mean(B_N[i][j][k][int(len(Time)/5)::])
#             UMean+=np.mean(U_N[i][j][k][int(len(Time)/5)::])
#         BMean/=len(B_N[i][j])
#         print(BMean)
#         UMean/=len(U_N[i][j])
#         print(UMean)
        
#         B_err=stats.sem((B_N[i][j][0:,0::dn]+U_N[i][j][0:,0::dn])/(BMean+UMean)*100)
#         Y=(np.mean(B_N[i][j][0:,0::dn], axis=0)+np.mean(U_N[i][j][0:,0::dn], axis=0))/(BMean+UMean)*100
#         #.errorbar(Time[0::dn]-duration/2/60,Y, yerr=[B_err/2,B_err/2], label=r'P={0:1.0f}, $\langle B_0 \rangle$={1:1.0f}'.format(N**2, BMean))
        
#         fig=plt.figure(figsize=(3.5,2), dpi=150)
#         plt.plot(np.array(Time[0::dn])/60,Y, color=col[i], linewidth=2, label=r'P={0:1.0f}, $\langle B \rangle$={1:1.0f}'.format(N**2, BMean))
        
#         plt.text(0.3,0.55,r'$\langle U/A_{spine} \rangle$='+'{0:1.1f} $\#/ \mu m^2$'.format(UMean/A_spine), transform=fig.axes[0].transAxes)
                
        
#         plt.axvline(int(len(Time)/5)*0.5/60)
#         plt.axhline(100, color='k')
#         # plt.xlim(0,35)
#         # plt.ylim(0,120)
#         plt.xlabel('Time (min)')
#         plt.ylabel('AMPAR \n Recovery (%)')
#         #plt.legend(ncol=2, prop=fontLgd)
#         sns.despine()
#         fig.tight_layout()
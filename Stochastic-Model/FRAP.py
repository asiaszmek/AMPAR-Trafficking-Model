# -*- coding: utf-8 -*-
"""
Created on Thu Nov  5 15:35:31 2020

@author: Moritz
"""

import numpy as np
import sm as sm

#%%



def FRAP(N_List, UFP_List, beta, alpha, kUB, kBU, duration, Nr_Trials):

    """Returns evolution over time of photobleached and not photobleached mobile and bound receptors (FRAP simulation).

    Parameters
    ----------
    N_List : array_like
        List of PSD sizes P
    UFP_List : array_like
        List of mobile receptor pool fixed points. Sets the influx of receptors into spine. 
    beta : float
        cooperativity factor for the unbinding. (Should be set to 1 or 0)
    alpha : float
        cooperativity factor for the binding
    kUB : float
        bidning rate
    kBU : float
        unbidning rate
    duration : float
        Duration of the simulation
    Nr_Trials : integer
        Number of trials.

    Returns
    -------
    B_N : array_like
        Time evolution of bound receptors for different conditions and number of Trials; shape(len(N_List), len(UFP_List), Trials, duration/0.5+1)
    U_N : array_like
        Time evolution of mobile receptors for different conditions and number of Trials; shape(len(N_List), len(UFP_List), Trials, duration/0.5+1)
    B_notBleached_N
        Time evolution of bound receptors not photobleached for different conditions and number of Trials; shape(len(N_List), len(UFP_List), Trials, duration/0.5+1)
    U_notBleached_N
        Time evolution of mobile receptors not photobleached for different conditions and number of Trials; shape(len(N_List), len(UFP_List), Trials, duration/0.5+1)
    PSD : array_like
        Matrix representing the PSD grid and its receptors at the end of the simulation (bleached, not bleached).
    Time : arrayl_like
        Time; shape(duration/0.5+1,)
    
    """
    
    A_spine_basal=0.898
    
    kBU=0.1
    kout=0.018
    kin=0.02
    #dt=0.5#the size of each time step
    
    t_bleaching=duration/2# 1 #
    
    
    #%%
    
    B_N=[]
    U_N=[]
    B_notBleached_N=[]
    U_notBleached_N=[]
    
    ID_basal=1
    ID_notBleached=2
    
    for N in N_List:
        
        A_spine=N**2/70*A_spine_basal
        
        B_U=[]
        U_U=[]
        B_notBleached_U=[]
        U_notBleached_U=[]
    
        for UFP_0 in UFP_List:
            
            dt=sm.calcTimeStep(UFP_0,A_spine,kUB,alpha,kBU,kout,kin)
            
            B_Tr=[]
            U_Tr=[]
            B_notBleached_Tr=[]
            U_notBleached_Tr=[]
            
            for Trial in range(0,Nr_Trials):
                
                if Trial%10==0:
                    print('Trial:',Trial)
                    print('P:',N**2, 'U:', UFP_0)
                    
                
                PSD=np.zeros((N,N))
    
                UFP=UFP_0
                UFP_notBleached=0
                
                Time=[]
                B_t=[]
                U_t=[]
                B_notBleached_t=[]
                U_notBleached_t=[]
    
                U=UFP_0
                U_notBleached=0
    
                for t in np.arange(0,duration+dt,dt):
                    
                    if t>t_bleaching:
                        UFP=0
                        UFP_notBleached=UFP_0
                    
                    NN=sm.nearestNeighbours(PSD)
                    
                    Mbu=sm.kBUcoop(kBU, NN, PSD, ID_basal, beta)*dt
                    Mub=sm.kUBcoop(kUB*U/A_spine, NN, PSD, alpha)*dt
                    Mbu_notBleached=sm.kBUcoop(kBU, NN, PSD, ID_notBleached, beta)*dt
                    Mub_notBleached=sm.kUBcoop(kUB*U_notBleached/A_spine, NN, PSD, alpha)*dt
    
                    PSD,dBoff,dBon,dBoff_notBleached,dBon_notBleached=sm.probabilityEval(Mub,Mbu,PSD,ID_basal,Mub_notBleached,Mbu_notBleached,ID_notBleached)
                    
                    pout=kout*U/A_spine*dt
                    pout_notBleached=kout*U_notBleached/A_spine*dt
                    pin=kin*UFP*dt
                    pin_notBleached=kin*UFP_notBleached*dt
                    U,U_notBleached=sm.update_mobilePool(U,pin,pout,dBoff,dBon, U_notBleached,pin_notBleached,pout_notBleached,dBoff_notBleached,dBon_notBleached)
                    
                    if t%0.5==0:
                        Time.append(t)
                        B_t.append(np.sum(PSD==ID_basal))
                        U_t.append(U)
                        B_notBleached_t.append(np.sum(PSD==ID_notBleached))
                        U_notBleached_t.append(U_notBleached)
                    
                B_Tr.append(B_t)
                U_Tr.append(U_t)
                B_notBleached_Tr.append(B_notBleached_t)
                U_notBleached_Tr.append(U_notBleached_t)
                
            B_U.append(B_Tr)
            U_U.append(U_Tr)
            B_notBleached_U.append(B_notBleached_Tr)
            U_notBleached_U.append(U_notBleached_Tr)
    
        B_N.append(B_U)
        U_N.append(U_U)
        B_notBleached_N.append(B_notBleached_U)
        U_notBleached_N.append(U_notBleached_U)
        
    return np.array(B_N), np.array(U_N), np.array(B_notBleached_N), np.array(U_notBleached_N), PSD, np.array(Time)/60
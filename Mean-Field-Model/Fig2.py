# -*- coding: utf-8 -*-
"""
Created on Tue Nov  3 16:37:38 2020

@author: Moritz
"""

"""Fig 2

This script reproduces the plots seen in Fig 2 of "The biophysical basis underlying the maintenance of early phase long-term potentiation".

This script requires that `numpy`,`scipy.integrate`,`matplotlib` and `seaborn` are installed within the Python environment you are running this script in.

This file can also be imported as a module and contains the following
functions:

    * Stim_Resp - returns a function describing the response of a parameter to a stimulus.
    * DV - returns a function describing the spine volume change during sLTP.
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
fontLgd.set_size('x-small')


#%%Define functions for parameter stimulus response and structural plasticity:


# def main():

SaveFig=0

#Color palettes used in plots:
    
col1=sns.dark_palette(sns.color_palette("colorblind", 10)[0],4)[3]
col2=sns.dark_palette(sns.color_palette("colorblind", 10)[1],4)[3]

#%%

#Import experimental Data for comparison with the model. Taken from Penn et al 2017 and Barco et al. 2002 (see Paper for reference details):

Penn2017_x=np.loadtxt("Data\\Penn2017_x.csv", delimiter=",")
Penn2017_y=np.loadtxt("Data\\Penn2017_y.csv", delimiter=",")
Barco2002_x=np.loadtxt("Data\\Barco2002_x.csv", delimiter=",")
Barco2002_y=np.loadtxt("Data\\Barco2002_y.csv", delimiter=",")

#%%

#Set parameter values and initial conditions:

Init = [10,20,13]

P=70
kin=0.2
kout=0.018
kexo0=0.0018
kendo=0.0021
kBU=0.1
kin_RE=0.1
kout_RE=0.000615
V0=0.08

t=np.arange(0,130*60)

#%%

#Basic AMPAR trafficking

sLTP=0
Cooperativity=0
kUB0=0.0036

#E-LTP with exocytosis:
kUB=mf.Parameter(kUB0)
kUB.timecourse(mf.Stim_Resp,[1,30,5,60])
kexo=mf.Parameter(kexo0)
kexo.timecourse(mf.Stim_Resp,[1,5,25,60])
Vspine=mf.Parameter(V0)
Vspine.timecourse(mf.DV,[V0,True])

Model=mf.Model_system()

solve,infodict = odeint(Model.odes,Init,t,args=(Vspine,P,kin,kout,kexo,kendo,kUB,kBU,Cooperativity,sLTP,kin_RE,kout_RE),full_output=True)

#E-LTP without exocytosis:
kexo=mf.Parameter(0)
kexo.timecourse(mf.Stim_Resp,[1,0,25,60])
Vspine=mf.Parameter(V0)
Vspine.timecourse(mf.DV,[V0,False])

Model=mf.Model_system()

solve2,infodict2 = odeint(Model.odes,Init,t,args=(Vspine,P,kin,kout,kexo,kendo,kUB,kBU,Cooperativity,sLTP,kin_RE,kout_RE),full_output=True)


fig=plt.figure(figsize=(3.5,2), dpi=150)
sns.lineplot(Barco2002_x,Barco2002_y, color='k', marker="s", ms=5, label='E-LTP', legend=False, zorder=0)
sns.lineplot(Penn2017_x,Penn2017_y, color='k', marker="v", ms=5, label='LTP, $k_{exo}=0$', legend=False, zorder=0)
plt.plot(t/60,solve.T[1]/20*100, color=col1, linewidth=2, label='Basic Model')
plt.plot(t/60,solve2.T[1]/20*100, color=col2, linewidth=2, label='Basic Model, $k_{exo}=0$')
lgd=plt.legend(bbox_to_anchor=(0.05,0.65,0,0), loc="lower left", mode="normal", borderaxespad=1.5, ncol=2, prop=fontLgd)
plt.axhline(100, color='k', linestyle='--', linewidth=0.5)
plt.xlabel('Time (min)')
plt.ylabel('Bound AMPARs $B$ (%)')
plt.xlim(-3,130)
plt.ylim(40)
sns.despine()
fig.tight_layout()
if SaveFig==1:
    print('Save')
    fig.savefig('Figures\\Fig2A.png', bbox_inches="tight", dpi=400)
    fig.savefig('Figures\\Fig2A.svg', bbox_inches="tight", dpi=400)


#%%

#sLTP

sLTP=1
Cooperativity=0
kUB0=0.0036

#E-LTP with exocytosis:
kUB=mf.Parameter(kUB0)
kUB.timecourse(mf.Stim_Resp,[1,30,5,60])
kexo=mf.Parameter(kexo0)
kexo.timecourse(mf.Stim_Resp,[1,5,25,60])
Vspine=mf.Parameter(V0)
Vspine.timecourse(mf.DV,[V0,True])

Model=mf.Model_system()

solve,infodict = odeint(Model.odes,Init,t,args=(Vspine,P,kin,kout,kexo,kendo,kUB,kBU,Cooperativity,sLTP,kin_RE,kout_RE),full_output=True)

#E-LTP without exocytosis:
kexo=mf.Parameter(0)
kexo.timecourse(mf.Stim_Resp,[1,0,25,60])
Vspine=mf.Parameter(V0)
Vspine.timecourse(mf.DV,[V0,False])

Model=mf.Model_system()

solve2,infodict2 = odeint(Model.odes,Init,t,args=(Vspine,P,kin,kout,kexo,kendo,kUB,kBU,Cooperativity,sLTP,kin_RE,kout_RE),full_output=True)


fig=plt.figure(figsize=(3.5,2), dpi=150)
sns.lineplot(Barco2002_x,Barco2002_y, color='k', marker="s", ms=5, legend=False, zorder=0)
sns.lineplot(Penn2017_x,Penn2017_y, color='k', marker="v", ms=5, legend=False, zorder=0)
plt.plot(t/60,solve.T[1]/20*100, color=col1, linewidth=2, label='sLTP Model')
plt.plot(t/60,solve2.T[1]/20*100, color=col2, linewidth=2, label='sLTP Model, $k_{exo}=0$')
lgd=plt.legend(bbox_to_anchor=(0.4,0.65,0,0), loc="lower left", mode="normal", borderaxespad=1.5, ncol=1, prop=fontLgd)
plt.axhline(100, color='k', linestyle='--', linewidth=0.5)
plt.xlabel('Time (min)')
plt.ylabel('Bound AMPARs $B$ (%)')
plt.xlim(-3,130)
plt.ylim(40)
sns.despine()
fig.tight_layout()
if SaveFig==1:
    print('Save')
    fig.savefig('Figures\\Fig2B.png', bbox_inches="tight", dpi=400)
    fig.savefig('Figures\\Fig2B.svg', bbox_inches="tight", dpi=400)

#%%

#Cooperative receptor binding

sLTP=0
Cooperativity=1
kUB0=0.0005#0.0036#

#E-LTP with exocytosis:
kUB=mf.Parameter(kUB0)
kUB.timecourse(mf.Stim_Resp,[1,5,5,60])
kexo=mf.Parameter(kexo0)
kexo.timecourse(mf.Stim_Resp,[1,5,25,60])
Vspine=mf.Parameter(V0)
Vspine.timecourse(mf.DV,[V0,True])

Model=mf.Model_system()

solve,infodict = odeint(Model.odes,Init,t,args=(Vspine,P,kin,kout,kexo,kendo,kUB,kBU,Cooperativity,sLTP,kin_RE,kout_RE),full_output=True)

#E-LTP without exocytosis:
kexo=mf.Parameter(0)
kexo.timecourse(mf.Stim_Resp,[1,0,25,60])
Vspine=mf.Parameter(V0)
Vspine.timecourse(mf.DV,[V0,False])

Model=mf.Model_system()

solve2,infodict2 = odeint(Model.odes,Init,t,args=(Vspine,P,kin,kout,kexo,kendo,kUB,kBU,Cooperativity,sLTP,kin_RE,kout_RE),full_output=True)


fig=plt.figure(figsize=(3.5,2), dpi=150)
sns.lineplot(Barco2002_x,Barco2002_y, color='k', marker="s", ms=5, legend=False, zorder=0)
sns.lineplot(Penn2017_x,Penn2017_y, color='k', marker="v", ms=5, legend=False, zorder=0)
plt.plot(t/60,solve.T[1]/20*100, color=col1, linewidth=2, label='Coop. Model')
plt.plot(t/60,solve2.T[1]/20*100, color=col2, linewidth=2, label='Coop. Model, $k_{exo}=0$')
lgd=plt.legend(bbox_to_anchor=(0.4,0.65,0,0), loc="lower left", mode="normal", borderaxespad=1.5, ncol=1, prop=fontLgd)
plt.axhline(100, color='k', linestyle='--', linewidth=0.5)
plt.xlabel('Time (min)')
plt.ylabel('Bound AMPARs $B$ (%)')
plt.xlim(-3,130)
plt.ylim(40)
sns.despine()
fig.tight_layout()
if SaveFig==1:
    print('Save')
    fig.savefig('Figures\\Fig2C.png', bbox_inches="tight", dpi=400)
    fig.savefig('Figures\\Fig2C.svg', bbox_inches="tight", dpi=400)


#%%

#sLTP + cooperative receptor binding

sLTP=1
Cooperativity=1
kUB0=0.0005#0.0036#

#E-LTP with exocytosis:
kUB=mf.Parameter(kUB0)
kUB.timecourse(mf.Stim_Resp,[1,5,5,60])
kexo=mf.Parameter(kexo0)
kexo.timecourse(mf.Stim_Resp,[1,5,25,60])
Vspine=mf.Parameter(V0)
Vspine.timecourse(mf.DV,[V0,True])

Model=mf.Model_system()

solve,infodict = odeint(Model.odes,Init,t,args=(Vspine,P,kin,kout,kexo,kendo,kUB,kBU,Cooperativity,sLTP,kin_RE,kout_RE),full_output=True)

#E-LTP without exocytosis:
kexo=mf.Parameter(0)
kexo.timecourse(mf.Stim_Resp,[1,0,25,60])
Vspine=mf.Parameter(V0)
Vspine.timecourse(mf.DV,[V0,False])

Model=mf.Model_system()

solve2,infodict2 = odeint(Model.odes,Init,t,args=(Vspine,P,kin,kout,kexo,kendo,kUB,kBU,Cooperativity,sLTP,kin_RE,kout_RE),full_output=True)


fig=plt.figure(figsize=(3.5,2), dpi=150)
sns.lineplot(Barco2002_x,Barco2002_y, color='k', marker="s", ms=5, legend=False, zorder=0)
sns.lineplot(Penn2017_x,Penn2017_y, color='k', marker="v", ms=5, legend=False, zorder=0)
plt.plot(t/60,solve.T[1]/20*100, color=col1, linewidth=2, label='sLTP+Coop. Model')
plt.plot(t/60,solve2.T[1]/20*100, color=col2, linewidth=2, label='sLTP+Coop. Model, $k_{exo}=0$')
lgd=plt.legend(bbox_to_anchor=(0.3,0.65,0,0), loc="lower left", mode="normal", borderaxespad=1.5, ncol=1, prop=fontLgd)
plt.axhline(100, color='k', linestyle='--', linewidth=0.5)
plt.xlabel('Time (min)')
plt.ylabel('Bound AMPARs $B$ (%)')
plt.xlim(-3,130)
plt.ylim(40)
sns.despine()
fig.tight_layout()
if SaveFig==1:
    print('Save')
    fig.savefig('Figures\\Fig2D.png', bbox_inches="tight", dpi=400)
    fig.savefig('Figures\\Fig2D.svg', bbox_inches="tight", dpi=400)


# if __name__ == "__main__":
#     main()

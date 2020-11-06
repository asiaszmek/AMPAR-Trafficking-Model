# -*- coding: utf-8 -*-
"""
Created on Wed Nov  4 11:43:38 2020

@author: Moritz
"""


import numpy as np


class Parameter():
    """This class defines some basic properties and functionality for parameters that change during LTP induction.
    
    Attributes
    ----------
    base_value : float
        Base value of the correspconding model parameter, i.e. before LTP-induction.

    Methods
    -------
    timecourse(function,pars)
        Sets up time course of the parameter during LTP-induction. Requires a function and the parameters of that function to be passed.
    update(t)
        Updates the parameter to the current value according to the function passed to timecourse.
    """

    def __init__(self,base_value):
        """    
        Parameters
        ----------
        base_value : float
            Base value of the correspconding model parameter, i.e. before LTP-induction.
        """
        self.base_value=base_value
        self.current_value=base_value
    
    def timecourse(self,function,pars):
        
        """
        Sets up time course of the parameter during LTP-induction. Requires a function and the parameters of that function to be passed.
        
        Parameters
        ----------
        function : function
            A function that describes the time course of the parameter during LTP-induction
        pars : array_like
            List of parameters for the function.
            
    
        Returns
        -------
        out: 
        """
            
        self.pars=np.array(pars,dtype=float)
        self.course=function
        
    def update(self,t):
        """
        Sets up time course of the parameter during LTP-induction. Requires a function and the parameters of that function to be passed.
        
        Parameters
        ----------
        t : float
            Time passed since the LTP induction stimulus.
    
        Returns
        -------
        out: Updates current value of the partameter.
        """
        self.current_value=self.base_value*self.course(t,*self.pars)
        


class Model_system():
    """This class contains the set of differential equations that describe AMPAR trafficking at spines.
    
    Attributes
    ----------

    Methods
    -------
    odes(Init,t,Vspine,P,kin,kout,kexo,kendo,kUB,kBU,Cooperativity,sLTP,kin_RE,kout_RE)
        Defines the ODEs that describe the synapse.
    """    
    def __init__(self):
        """    
        Parameters
        ----------
        container : array_like
            container that can be filled with any parameter one is interested in during numerical integration, e.g. t, kexo or Vspine etc..
        """
        self.container=[]
        
    def odes(self,Init,t,Vspine,P,kin,kout,kexo,kendo,kUB,kBU,Cooperativity,sLTP,kin_RE,kout_RE):
        
        """
        Defines the ODEs that describe the AMPAR dynamics at the spine.
        
        Parameters
        ----------
        Init : [U(0),B(0),Sexo(0)]
            Initial conditions for the three variables of the system.
        t : array_like
            Time.
        Vspine : float
            Describes the spine volume and spine volume change during E-LTP.
        P : float
            Number of binding site/slots at the PSD.
        kin : float
            Rate at which receptors hop from the dednritic membrane compartment onto the spine membrane compartment.
        kout : float 
            Rate at which receptors hop from the spine membrane compartment onto the dendritic membrane compartment.
        kexo : float
            Rate of exocytosis events occuring at the spine.
        kendo : float
            Rate at which receptors are endocytosed at the spine.
        kUB : float
            Rate at which AMPARs bind to PSD slots.
        kBU : float
            Rate at which AMPARs unbind from PSD slots.
        Cooperativity: 1,0
            States whether binding is cooperative (1) or not cooperative (0).
        sLTP : 1,0
            States whether slTP does occur (1) or does not occur (0).
        kin_RE : float
            Rate at which AMPAR containing endosomes enter the spine (dynamics of exocytosis event size Sexo).
        kout_RE : float
            Rate at which AMPAR containing endosomes leave the spine.
            
            
    
        Returns
        -------
        out: [dU,dB,dS_exo]
        """
    
        self.container.append([Vspine.current_value,t])
        
        U=Init[0]
        B=Init[1]
        S_exo=Init[2]
    
    
        def kUB_(B,kUB,Cooperativity,P):
            if Cooperativity==1:       
                m = 24.6/(12.5 + P);
                return kUB*(m*B**0.8 + 1);
            else:
                return kUB
            
        def kBU_(B,kBU,Cooperativity,P):
            if Cooperativity==1:
                Lambda = 0.9*P + 15.5;
                Beta = 0.6*P + 9.3;
                return kBU*(Lambda/(Beta + B)-0.5);
            else:
                return kBU
            
        kexo.update(t)
        kUB.update(t)
        
        if sLTP==1:
            Vspine.update(t)
        
        Aspine=4*np.pi*(3*Vspine.current_value/(4*np.pi))**(2/3)

        dS_exo=kin_RE-kout_RE*S_exo/Vspine.current_value
        dU=kexo.current_value*S_exo+kin+kBU_(B,kBU,Cooperativity,P)*B-kendo*U/Aspine-kout*U/Aspine-kUB_(B,kUB.current_value,Cooperativity,P)*(P-B)*U/Aspine
        dB=kUB_(B,kUB.current_value,Cooperativity,P)*(P-B)*U/Aspine-kBU_(B,kBU,Cooperativity,P)*B
    
        return [dU,dB,dS_exo]



def Stim_Resp(t,a0,a,b,c):
    """This function describes the response of a parameter to the LTP induction-stimulus.

    Parameters
    ----------
    t : float
        time passed since induction of LTP
    a0 : float
        baseline of the parameter (Should usually be set to 1). 
    a : float
        Parameter amplitude during LTP induction
    b : float
        time constant
    c : float
        time constant

    Returns
    -------
    float
        current factor of parameter change.
        
    Examples
    --------
    
    Import libraries:
        
    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> import mf as mf
    >>> import seaborn as sns
    
    Set parameters:

    >>> t=np.linspace(0,10*60,100)
    >>> kexo0=0.0018
    >>> exocytosis=True
    
    Plot:
        
    >>> plt.figure(figsize=(4,3), dpi=150)
    >>> plt.plot(t/60,kexo0*mf.Stim_Resp(t,1,5,25,60), linewidth=2)
    >>> plt.xlabel('Time (min.)')
    >>> plt.ylabel('$k_{exo}$ ($s^{-1}$)')
    >>> sns.despine()

    Output:
        
    .. image:: images/example_Stim_Resp.png
        :width: 70%
    """
    
    A=(b/c)**(c/(c-b))-(b/c)**(b/(c-b))
    return a0+a*(np.exp(-t/b)-np.exp(-t/c))/A



def DV(t,V0,exocytosis):
    """This function describes the evolution of the spine volume during sLTP.

    Parameters
    ----------
    t : float
        time passed since induction of LTP
    V0 : float
        initial spine volume. 
    exocytosis : True, False
        States whether exocytosis is blocked or not. Make sure that the boolean value is chosen in agreement with the definition of kexo.

    Returns
    -------
    float
        current factor of spine volume change. 
        
    Examples
    --------
    
    Import libraries:
        
    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> import mf as mf
    >>> import seaborn as sns
    
    Set parameters:

    >>> t=np.linspace(0,120*60,1000)
    >>> exocytosis=True
    >>> spineVolumes=[0.04,0.08,0.2]
    
    Plot:
        
    >>> plt.figure(figsize=(4,3), dpi=150)
    >>> for V0 in spineVolumes:
    >>>     plt.plot(t/60,mf.DV(t,V0,exocytosis)*100, linewidth=2, label='$V_{spine}^0=$'+'{:1.2f} $\mu m^3$'.format(V0))
    >>> plt.legend()
    >>> plt.xlabel('Time (min.)')
    >>> plt.ylabel('$V_{spine}$ (%)')
    >>> sns.despine()

    Output:
        
    .. image:: images/example_DV.png
        :width: 70%
    """
    
    def DV_long_(x):
        a=299
        b=0.102
        c=-0.9
        return (a*np.exp(-x/b)+c)/100
    
    if exocytosis==True:
        DVlong=DV_long_(V0)
    else:
        DVlong=0
        
    a0_1=1
    a_1=1+1.75*DVlong
    b_1=250
    c_1=5
    
    a0_2=0
    a_2=DVlong
    b_2=3000
    c_2=500
    
    return Stim_Resp(t,a0_1,a_1,b_1,c_1)+Stim_Resp(t,a0_2,a_2,b_2,c_2)
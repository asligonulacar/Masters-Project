#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 25 14:59:37 2021

@author: asligonulacar
"""

from FourVector import fv
from ME_squared import ME
from Generator import Gen
from Amplitude import amp
import numpy as np
from Channel_Basics import CB
import matplotlib.pyplot as plt

mx=0.
my=0.
m1=0.
m2=0.
ctmin=-1
ctmax=1
cR=cL=cR_p=cL_p=amp.e

class CS:
    def MC_int_Weighted(ch,N,mx,my,m1,m2,ctmin,ctmax,cR,cL,cR_p,cL_p):   
        sigma=[]
        cross=[]
        ct=[]
        E_track=[]
        for i in range(N):
            E=np.random.uniform(1,30)
            E_track.append(E)
            p=Gen.Iso_Gen(E_track[i], mx, my, m1, m2, ctmin, ctmax)
            Iso=Gen.Iso_Weight(p[0],p[1],p[2],p[3],ctmin,ctmax)
            f=(Gen.Flux(p[0], p[1], m1, m2))
            a=amp.Z_new(p[0], p[1], p[2], p[3], cR, cL, cR_p, cL_p, ctmin, ctmax,ch)
            M=ME.MESqr(a,ch,p[0],p[1])
            MF=(ME.Ref(mx, my, m1, m2, p[0], p[1], p[2], p[3], ch)[0]) #reference ME
            sigma.append(abs(M*Iso*f))
            cross.append((4*np.pi*((1/137)**2))/(3*E_track[i]**2)) #reference cross section       
            #error=np.std(sigma)/np.sqrt(N)
        
        sig_max=np.amax(sigma)
        sigma_tot=sum(sigma)/N
        plt.scatter(E_track,sigma)
        plt.ylabel('σ (cm)')
        plt.xlabel('√s (eV)')
        plt.xlim(0,30)
        plt.title('Total cross section versus CMS energy')
        
        #plt.colorbar()
        plt.show()
        cross=sum(cross)/N 
        return cross/sigma_tot
    
    
        
print(CS.MC_int_Weighted('s',1000,mx,my,m1,m2,ctmin,ctmax,cR,cL,cR_p,cL_p))  


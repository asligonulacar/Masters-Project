#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  5 16:50:42 2022

@author: asligonulacar
"""

from FourVector import fv
from Generator import Gen
from Amplitude import amp
import numpy as np
import matplotlib.pyplot as plt
import cmath

E=10
mx=0
my=0
m1=0
m2=0
ctmin=-1
ctmax=0.99
ch='s'
cL=cL_p=0.653
 
def MC_int_Weighted(E,ch,mx,my,m1,m2,ctmin,ctmax,cR,cL,cR_p,cL_p):   
        err=0
        #sigma=np.array([])
        #sigma2=[]
        #cross=np.array([])
        N=1000
        sigma=[]
        M=[]
        # Gm=[]
        # Ge=[]
        # Mref=[]
        A=[]
        while err==0 or err>10000:#0.05*(sum(sigma2)/len(ctl)):
            N=N+5000
            for i in range(N):
                if ((i/N)*100)%10==0:
                    print((str((i/N)*100)) +'%')
                p=Gen.Iso_Gen(E, mx, my, m1, m2, ctmin, ctmax)
                A.append(amp.Electron_Neutrino_amp(p[2], p[1], p[3], p[0], 0, cL, 0, cL_p))
                P=Gen.Iso_Weight(p[0], p[1], p[2], p[3], ctmin, ctmax)
                f=Gen.Flux_CMS(p[0], p[1], mx, my)
                #sigma.append(a*P*f*(2*np.pi**2))
                s=fv.Pro(fv(fv.Add(p[0],p[1])),fv(fv.Add(p[0],p[1])))
                M.append((0.653**4/(2*80**4))*s**2)
            err=1    
            #err=np.std(sigma2)/(N**0.5)
        #s_tot=((1/(8*np.pi))*(0.653**4)*((p[0].x[0]**2)/(80**4))*(1-((0.1)**2/(4*p[0].x[0]**2)))**2)
        #sigma_tot=sum(sigma)/N
        atot=sum(A)/N
        mtot=sum(M)/N
        
        return atot/mtot#s_tot/sigma_tot

print(MC_int_Weighted(E,ch,mx,my,m1,m2,ctmin,ctmax,0,cL,0,cL_p))
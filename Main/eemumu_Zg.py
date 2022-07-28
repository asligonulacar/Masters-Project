#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  5 19:55:25 2022

@author: asligonulacar
"""
from FourVector import fv
from Generator import Gen
from Amplitude import amp
import numpy as np
import matplotlib.pyplot as plt
import cmath

E=1
mx=0
my=0
m1=0
m2=0
ctmin=-1
ctmax=0.99
ch='s'
cL=cL_p=amp.e
 
def MC_int_Weighted(E,ch,mx,my,m1,m2,ctmin,ctmax,cR,cL,cR_p,cL_p):   
        err=0
        #sigma=np.array([])
        #sigma2=[]
        #cross=np.array([])
        N=1000
        sigma=[]
        S=[]
        # Gm=[]
        # Ge=[]
        # Mref=[]
        # A=[]
        while err==0 or err>10000:#0.05*(sum(sigma2)/len(ctl)):
            N=N+5000
            for i in range(N):
                if ((i/N)*100)%10==0:
                    print((str((i/N)*100)) +'%')
                p=Gen.Iso_Gen(E, mx, my, m1, m2, ctmin, ctmax)
                a=(amp.E_mu_amp(p[0], p[1], p[2], p[3]))
                P=Gen.Iso_Weight(p[0], p[1], p[2], p[3], ctmin, ctmax)
                f=Gen.Flux_CMS(p[0], p[1], mx, my)
                s=fv.Pro(fv(fv.Add(p[0],p[1])),fv(fv.Add(p[0],p[1])))
                sigma.append(a*P*f*2*np.pi**2)
                G=np.sqrt(2)*0.653**2/(8*80**2)
                A=1+(abs(((np.sqrt(2)*G*80**2)/(s-80**2))*(s/amp.e**2))**2)/16
                S.append((1/137)**2*4*np.pi*A/(3*s))
                
               
                
                #M.append(2*((0.653)**4/(80**4))*fv.Pro(p[0],p[1])*fv.Pro(p[2],p[3]))
            err=1    
            #err=np.std(sigma2)/(N**0.5)
        s_tot=sum(S)/N
        sigma_tot=sum(sigma)/N
        
        
        return s_tot/sigma_tot

print(MC_int_Weighted(E,ch,mx,my,m1,m2,ctmin,ctmax,0,cL,0,cL_p))
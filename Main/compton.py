#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  5 00:48:44 2022

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
        M=[]
        a=[]
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
                a.append(amp.Compton_amp(p[0], p[2], p[3], p[1], cR, cL, cR_p, cL_p))
                s=fv.Pro(fv(fv.Add(p[0],p[1])),fv(fv.Add(p[0],p[1])))
                u=fv.Pro(fv(fv.Sub(p[0],p[3])),fv(fv.Sub(p[0],p[3])))
                M.append(2*amp.e**4*(-u/s))
            err=1    
            #err=np.std(sigma2)/(N**0.5)
        M_tot=sum(M)/N
        a_tot=sum(a)/N
        
        return a_tot/M_tot

print(MC_int_Weighted(E,ch,mx,my,m1,m2,ctmin,ctmax,0,cL,0,cL_p))
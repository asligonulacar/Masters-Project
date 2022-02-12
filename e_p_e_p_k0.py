#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  3 15:34:11 2022

@author: asligonulacar
"""
from FourVector import fv
from ME_squared import ME
from Generator import Gen
from Amplitude import amp
import numpy as np
from Channel_Basics import CB 
import matplotlib.pyplot as plt
E=50
me=0
mp=0.938
ctmin=-1
ctmax=1-1e-7
ch='t'
cR=cR_p=amp.e
cL=cL_p=amp.e
def MC_int_Weighted(E,ch,mx,my,m1,m2,ctmin,ctmax,cR,cL,cR_p,cL_p):   
        err=0
        sigma=[]
        N=1000
        cross=[]
        cross2=[]
        ct=[]
        #En=[]
        for i in range(N):
            p=Gen.Iso_Gen(E, mx, my, m1, m2, ctmin, ctmax)
            ct.append(p[4])
            #convenience for ME calculations
            pm2=fv([p[1].x[0],-p[1].y[0],-p[1].y[1],-p[1].y[2]]) #-p2 
            pm3=fv([p[2].x[0],-p[2].y[0],-p[2].y[1],-p[2].y[2]]) #-p3
            Iso=Gen.Iso_Weight(p[0],p[1],p[2],p[3],ctmin,ctmax)
            f=(Gen.Flux(p[0], p[1], m1, m2))
            a=amp.Z_new(p[0], pm3, p[3], pm2, cR, cL, cR_p, cL_p, ctmin, ctmax,ch)
            M=ME.MESqr(a, ch, p[0], pm3, p[3], pm2)
            sigma.append(abs((M*Iso*f)))
            #s=fv.Pro(fv(fv.Add(p[0],pm2)),fv(fv.Add(p[0],pm2)))
            q=fv.Pro(fv(fv.Sub(p[0],pm3)),fv(fv.Sub(p[0],pm3)))
            Mref2=(8*amp.e**4/q**2)*(fv.Pro(pm3,p[3])*fv.Pro(p[0],pm2)+fv.Pro(pm3,pm2)*fv.Pro(p[0],p[3])+mp**2*fv.Pro(p[0],pm3))
            cross.append(Mref2/(16*np.pi*(pm2.x[0]+p[0].x[0])**2))
            #cross2.append((((1/137)*2/(16*p[0].x[0]**2*((1-ct[i])/2)**5/2))*(((1+ct[i])/2)-(q**2/(2*m2**2))*(1-ct[i])/2))*(1+((2*p[0].x[0]/m2)*((1-ct[i])/2)))**(-1))      
        #cross=sum(cross)/N
        stot=sum(sigma)/N
        ctot=sum(cross)/N
        ctot2=sum(cross2)/N
        return ctot/stot

print(MC_int_Weighted(E,ch,me,mp,me,mp,ctmin,ctmax,cR,cL,cR_p,cL_p))

#cross section for e mu e mu calculated successfully 




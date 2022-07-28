#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 30 16:25:22 2022

@author: asligonulacar
"""
from FourVector import fv
from ME_squared import ME
from Generator import Gen
from Amplitude import amp
import numpy as np
from Channel_Basics import CB 
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from numpy import arange
E=1
me=0
m=0
ctmin=-1
ctmax=1-1e-7
ch='s'
cR=cR_p=amp.e
cL=cL_p=amp.e
def objective(x, a, b, c):
	return a * x + b * x**2 +c

def MC_int_Weighted(E,ch,mx,my,m1,m2,ctmin,ctmax):   
        err=0
        sigma=[]
        N=10
        cross=[]
        #En=[]
        for i in range(N):
            p=Gen.Iso_Gen(E, mx, my, m1, m2, ctmin, ctmax)
            #convenience for ME calculations
            pm2=fv([p[1].x[0],-p[1].y[0],-p[1].y[1],-p[1].y[2]]) #-p2 
            pm3=fv([p[2].x[0],-p[2].y[0],-p[2].y[1],-p[2].y[2]]) #-p3
            Iso=Gen.Iso_Weight(p[0],p[1],p[2],p[3],ctmin,ctmax)
            f=(Gen.Flux_CMS(p[0], p[1], m1, m2))
            a=amp.Z(p[0], pm3, p[3], pm2, cR, cL, cR_p, cL_p)
            a=np.sum(a*np.conj(a))
            M=ME.MESqr(a, ch, p[0], pm3, p[3], pm2)
            sigma.append(abs((M*Iso*f)))
            #s=fv.Pro(fv(fv.Add(p[0],pm2)),fv(fv.Add(p[0],pm2)))
            q=fv.Pro(fv(fv.Add(p[0],pm3)),fv(fv.Add(p[0],pm3)))
            Mref2=(8*amp.e**4/q**2)*(fv.Pro(pm3,p[3])*fv.Pro(p[0],pm2)+fv.Pro(pm3,pm2)*fv.Pro(p[0],p[3])+m**2*fv.Pro(p[0],pm3))
            cross.append(Mref2/(16*np.pi*(pm3.x[0]+p[0].x[0])**2))
        er= np.var(sigma)
        stot=sum(sigma)/N
        ctot=sum(cross)/N
        return stot,ctot,er

#print(MC_int_Weighted(E,ch,me,m,me,m,ctmin,ctmax,cR,cL,cR_p,cL_p))

#cross section for e mu e mu calculated successfully 
def cross(Ex):
    Output= np.ones((len(Ex)))
    for i in range(len(Ex)):
        E=Ex[i]
        Output[i]=((4*np.pi*((1/137)**2)))/(3*E**2)
    return Output

def Energy(mx,my,m1,m2,ctmin,ctmax):
    En=[]
    sigma=[]
    m=[]
    err=[]
    for i in range(20):
        if ((i/10)*100)%10==0:
            print((str((i/10)*100)) +'%')
        E=np.random.uniform(1,5)
        En.append(E)
        sigma.append(MC_int_Weighted(E,ch,mx,my,m1,m2,ctmin,ctmax)[0])
        m.append(MC_int_Weighted(E,ch,mx,my,m1,m2,ctmin,ctmax)[1])
        err.append(MC_int_Weighted(E,ch,mx,my,m1,m2,ctmin,ctmax)[2])
   
    
    for i in range(len(sigma)):
            sigma[i]=0.27*10**5*sigma[i]
            m[i]=0.25*10**5*m[i]
    
    Ex=np.linspace(min(En),max(En),10)
   
    plt.grid(False)
    plt.xlabel('√s [GeV]')
    plt.ylabel('σ [nb]')
    plt.yscale('log')
    plt.xlim(0,6)
    #plt.ylim(0.5*10**-2,10)
    # x,y=En,m
    # popt, _ = curve_fit(objective, x, y)
    # a, b, c = popt
    # x_line = arange(min(x), max(x), 1)
    # y_line = objective(x_line, a, b, c)
    #plt.plot(x_line, y_line, color='black', label='analytical')
    plt.scatter(En,sigma,label='MC', color='red')
    plt.plot(Ex,0.4*10**5*cross(Ex),color='black',alpha=1,label='Analytical')
    #plt.scatter(En,m)
    plt.errorbar(En,sigma,yerr=err, ls='none')#,alpha=0.8,marker='^',markersize=8,label='Monte Carlo')
    plt.legend()
    return None


    
print(Energy(me,m,me,m,ctmin,ctmax))







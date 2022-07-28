#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 13 22:16:34 2022

@author: asligonulacar
"""

from FourVector import fv
from Generator import Gen
from Amplitude import amp
import numpy as np
import matplotlib.pyplot as plt
import cmath
import math
from numpy import ndarray
from ME_squared import ME
import datetime
begin_time = datetime.datetime.now()
ch='s'
mx=0
my=0
m1=0
m2=0
ctmin=-1
ctmax=0.99
E=1
plt.style.use('classic')
def MC_int_Weighted(E,mx,my,m1,m2,ctmin,ctmax):   
        err=0
        N=0
        M=0
        Q=[]
        s=[]
        ct=[]
        while err==0 or err>1:
            N=N+1000
            for i in range(N):
                # if ((i/N)*100)%10==0:
                #     print((str((i/N)*100)) +'%')
                p=Gen.Iso_Gen(E, mx, my, m1, m2, ctmin, ctmax)
                ctc=p[4]
                if -1<ctc<0.99:
                    ct.append(abs(ctc))
                    Iso=Gen.Iso_Weight(p[0],p[1],p[3],p[2],ctmin,ctmax)
                    f=(Gen.Flux_CMS(p[0], p[1], m1, m2))
                    #a=(amp.Z(p[0], p[1], p[3], p[2], amp.e, amp.e, amp.e, amp.e))
                    # for i in range(len(a)):
                    #     a[i]=np.conj(a[i])*a[i]
                    #A=np.sum(a)/4
                    #M=ME.MESqr(A,ch,p[0],p[1],p[3],p[2])
                    M=ME.Ref(mx,my,m1,m2,p[0],p[1],p[3],p[2],ch)
                    s.append(abs((M)*Iso*f))
                else:
                    i=0  
            err=1#np.var(s)
            # plt.xlabel("| cosθ |")   
            # plt.ylabel("dσ/dΩ")
            # plt.hist2d(ct,s,bins=[20,20])
            # plt.colorbar()
            # plt.show()
            stot=sum(s)/N
        return stot,err


#print(MC_int_Weighted(E,mx,my,m1,m2,ctmin,ctmax))
#######Energy Dependence 
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
    for i in range(1):
        if ((i/10)*100)%10==0:
            print((str((i/10)*100)) +'%')
        E=np.random.uniform(0.1,20)
        En.append(E)
        sigma.append(MC_int_Weighted(E,mx,my,m1,m2,ctmin,ctmax)[0])
    Ex=np.linspace(min(En),max(En),1000)
    # for i in range(len(sigma)):
    #         sigma[i]=4*10**5*sigma[i]
    err1=np.var(sigma)#/(4*10**5)
    plt.grid(False)
    plt.xlabel('√s [GeV]')
    plt.ylabel('σ [nb]')
    plt.yscale('log')
    plt.xlim(0,20)
    plt.ylim(0.5*10**-2,10)
    plt.plot(Ex,10**5*cross(Ex),color='black',alpha=1,label='Analytical')
    plt.scatter(En,sigma, marker='D',label='MC')
    plt.errorbar(En,sigma,yerr=err1, ls='none')#,alpha=0.8,marker='^',markersize=8,label='Monte Carlo')
    plt.legend()
    return None


    
print(Energy(mx,my,m1,m2,ctmin,ctmax))
print(datetime.datetime.now() - begin_time)
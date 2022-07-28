#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 31 16:14:28 2022

@author: asligonulacar
"""
from FourVector import fv
from Generator import Gen
from Amplitude import amp
import numpy as np
import matplotlib.pyplot as plt
import cmath
import datetime
begin_time = datetime.datetime.now()
E=1
mx=0
my=0
m1=0
m2=0
ctmin=-1
ctmax=0.99
ch='s'
cL=cL_p=0.653
 
def MC_int_Weighted(E,mx,my,m1,m2,ctmin,ctmax):   
        err=0
        #sigma=np.array([])
        #sigma2=[]
        #cross=np.array([])
        N=0
        sigma=[]
        M=[]
        # Gm=[]
        # Ge=[]
        # Mref=[]
        # A=[]
        while err==0 or err>10000:#0.05*(sum(sigma2)/len(ctl)):
            N=N+1000
            for i in range(N):
                # if ((i/N)*100)%10==0:
                #     print((str((i/N)*100)) +'%')
                p=Gen.Iso_Gen(E, mx, my, m1, m2, ctmin, ctmax)
                a=(amp.Inv_Mu_amp(p[2], p[0], p[3], p[1], 0, cL, 0, cL_p))
                P=Gen.Iso_Weight(p[0], p[1], p[2], p[3], ctmin, ctmax)
                f=Gen.Flux_CMS(p[0], p[1], mx, my)
                sigma.append(abs(a))
                
               
                
                M.append(2*((0.653)**4/(80**4))*fv.Pro(p[0],p[1])*fv.Pro(p[2],p[3]))
            err=1#np.var(sigma)
            #err=np.std(sigma2)/(N**0.5)
        #s_tot=((1/(8*np.pi))*(0.653**4)*((p[0].x[0]**2)/(80**4))*(1-((m1)**2/(4*p[0].x[0]**2)))**2)
        sigma_tot=sum(sigma)/N
        
        s_tot= sum(M)/N
        
        return s_tot/sigma_tot


print(MC_int_Weighted(1,mx,my,m1,m2,ctmin,ctmax))

def cross(Ex):
    Output= np.ones((len(Ex)))
    for i in range(len(Ex)):
        E=Ex[i]
        Output[i]=((1/(8*np.pi))*(0.653**4)*((E**2)/(80**4))*(1-((m1)**2/(4*E**2)))**2)
    return Output

def Energy(mx,my,m1,m2,ctmin,ctmax):
    En=[]
    sigma=[]
    m=[]
    err1=[]
    for i in range(1):
        if ((i/10)*100)%10==0:
            print((str((i/10)*100)) +'%')
        E=np.random.uniform(0.1,2)
        En.append(E/2)
        #sigma.append(MC_int_Weighted(E,mx,my,m1,m2,ctmin,ctmax)[1]*10**10)
        m.append(MC_int_Weighted(E,mx,my,m1,m2,ctmin,ctmax)[0])
        err1.append(MC_int_Weighted(E,mx,my,m1,m2,ctmin,ctmax)[1]*10**10)
        
    Ex=np.linspace(min(En),max(En),10)
    #plt.yscale('log')
    #plt.ylim(0,1.2*10**-38)
    plt.xlim(0,1.1)
    #err=np.var(sigma)
    plt.ylabel('σ [mb]')
    plt.xlabel('Ev [GeV]')
    #plt.ylim(0,1.3*max(sigma))
    plt.grid(False)
    plt.scatter(En,m, color='red')
    plt.plot(Ex,cross(Ex)*10**10,color='black',alpha=1,label='Analytical')
    #plt.plot(En,m, color='red')
    #plt.errorbar(En,sigma, yerr=err1, ls='none')
    
    
    plt.show()
    # plt.scatter(En,sigma)
    # plt.scatter(En,m)
    return 


    
#print(Energy(mx,my,m1,m2,ctmin,ctmax))
print(datetime.datetime.now() - begin_time)



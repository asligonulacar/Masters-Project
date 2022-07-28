#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 11 15:19:25 2021

@author: asligonulacar
"""
import matplotlib.pyplot as plt
from ME_squared import ME
from Generator import Gen
from Amplitude import amp
import numpy as np
E=1
mx=0
my=0
m1=0
m2=0
ctmin=-1
ctmax=1
cR=cL=cR_p=cL_p=amp.e

class CS:
    def MC_int_Iso(ch,E,mx,my,m1,m2,ctmin,ctmax,cR,cL,cR_p,cL_p):   
        err=0
        sigma=[]
        ct=[]
        dd=[0]
        cross=0
        N=1
        while err==0 or err>1e-8:
            N=N+1000
            for i in range(N):
                p=Gen.Iso_Gen(E, mx, my, m1, m2, ctmin, ctmax)
                ct.append((np.dot(p[1].y,p[2].y)/(np.linalg.norm(p[1].y)*np.linalg.norm(p[2].y))))
                Iso=Gen.Iso_Weight(p[0],p[1],p[2],p[3],ctmin,ctmax)
                f=(Gen.Flux(p[0], p[1], m1, m2))
                a=amp.Z_new(p[0], p[1], p[2], p[3], cR, cL, cR_p, cL_p, ctmin, ctmax,ch)
                M=ME.MESqr(a,ch,p[0],p[1])
                sigma.append(abs((M*Iso*f)))
                cross+=((4*np.pi*((1/137)**2))/(3*E**2))
            err=np.var(sigma)
        return sigma,ct,N,cross
    
    def Weighted_Event(ch,E,mx,my,m1,m2,ctmin,ctmax,cR,cL,cR_p,cL_p):
        pass
        
    
    def UnWeighted_Event(ch,E,mx,my,m1,m2,ctmin,ctmax,cR,cL,cR_p,cL_p):
        p_acc=[]
        sigma_tot=[]
        MC=CS.MC_int_Iso(ch,E,mx,my,m1,m2,ctmin,ctmax,cR,cL,cR_p,cL_p)
        sigma=MC[0]
        sigma_max=np.max(MC[0])
        N=MC[2]
        cross=MC[3]/N
        for i in range(N):
            R=np.random.uniform(0,1)
            p=sigma[i]/sigma_max
            sigtot=sum(sigma)/N
            if p<R:
                p_acc.append(p)
                sigma_tot.append(sigtot)
                
            else:
                p_acc.append(1-p)
                sigma_tot.append(0)
        efficiency=(sigtot)/sigma_max
        plt.hist2d(p_acc,sigma_tot,bins=[10,5])
        plt.colorbar()
        plt.show()
        return sigtot,efficiency
        
            
print(CS.UnWeighted_Event('s',E,mx,my,m1,m2,ctmin,ctmax,cR,cL,cR_p,cL_p))
        
        
        
         
    


#histogram, look for the angle between electron and muon for example
#fire 1 PS point (momentum generation) (lets call this j)
#p is the acceptance probability
#with p j add to histogram value of sigma total 
#with (1-p) acceptance with value 0.
#integrate it all up I should get sigma total back.
#if we normalise them correctly the two events should do the same.
#keep track of sigma max.
#expectation of value of p acceptance.
#anisotropic vs isotropic the unweighted efficiency should go up in t channel for aniso



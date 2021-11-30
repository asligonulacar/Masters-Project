#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 11 15:19:25 2021

@author: asligonulacar
"""
from FourVector import fv
from ME_squared import ME
from Generator import Gen
from Amplitude import amp
import numpy as np

# mx = 0
# my = 0
# m1=0
# m2=0
# E=100
# cR=amp.e   
# cL=amp.e   
# cR_p=amp.e       
# cL_p=amp.e   
# ctmin=-1
# ctmax=1
class CS:
    def MC_int(ch,N,E,mx,my,m1,m2,ctmin,ctmax,cR,cL,cR_p,cL_p):   
        sigma=[]
        cross=[]
        ct=[]
        theta=[]
        MF=[]
        for i in range(N):
            p=Gen.P_Gen(E, mx, my, m1, m2, ctmin, ctmax)
            #costheta=1-(fv.Pro(p[0], p[2])/(p[4]**2))
            ct.append(p[5])
            theta.append(p[6])
            Iso=Gen.Iso_Weight(p[0],p[1],p[2],p[3], ctmin, ctmax)
            f=Gen.Flux(p[0], p[1], m1, m2)
            a=amp.Z_new(p[0], p[1], p[2], p[3], cR, cL, cR_p, cL_p, ctmin, ctmax,ch)
            M=ME.MESqr(a,ch,p[0],p[1])
            MF.append(ME.Ref(mx, my, m1, m2, p[0], p[1], p[2], p[3], ch)[0]) #reference ME
            sigma.append((M*Iso*f))
            cross.append((4*np.pi*((1/137)**2))/(3*E**2)) #reference cross section       
            #error=np.std(sigma)/np.sqrt(N)
        sig_max=np.amax(sigma)
        sigma_tot=sum(sigma)/N
        cross=sum(cross)/N 

        #print('cross section         ','       error              ','   maximum of cross section')
        return sigma,ct,MF,theta
    
    # def ErrorMinimize(ch,N,E,mx,my,m1,m2,ctmin,ctmax,cR,cL,cR_p,cL_p):  
    #     SGM=CS.MC_int(ch, N, E, mx, my, m1, m2, ctmin, ctmax, cR, cL, cR_p, cL_p)
    #     error=SGM[1]
    #     while error>=1e-4:
    #         N=+100
    #         SGM
    #         error=SGM[1]
    #     return SGM[0],SGM[1],SGM[2],SGM[3],SGM[4]
        
#print(CS.MC_int('s',100,E,mx,my,m1,m2,ctmin,ctmax,cR,cL,cR_p,cL_p))  

# N=2
# def func(N):
#     for i in range(N):
#         i=N+1
#     return i

# def m(func,N):
#     i=func(N)
#     while i<=4:
#         N+=1
#         func(N)
#         i=func(N)
#     return i
        
#CS.ErrorMinimize('s', N, E, mx, my, m1, m2, ctmin, ctmax, cR, cL, cR_p, cL_p)
# print(CS.ErrorMinimize('s', N, E, mx, my, m1, m2, ctmin, ctmax, cR, cL, cR_p, cL_p)[1])
#check until the error is small.
    
#histogram, look for the angle between electron and muon for example
#fire 1 PS point (momentum generation) (lets call this j)
#p is the acceptance probability
#with p j add to histogram value of sigma total 
#with (1-p) acceptance with value 0.
#integrate it all up I should get sigma total back.
#if we normalise them correctky the two events should do the same.
#keep track of sigma max.
#expectation of value of p acceptance.
#anisotropic vs isotropic the unweighted efficiency should go up in t channel for aniso



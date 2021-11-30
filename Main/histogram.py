#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 21 18:28:21 2021

@author: asligonulacar
"""
from FourVector import fv
from Amplitude import amp
from CS_Integrator import CS
from Generator import Gen 
import numpy as np
import matplotlib.pyplot as plt

#WEIGHTED EVENTS
N=100
mx=0
my=0
m1=0
m2=0
E=1
cR=amp.e   
cL=amp.e   
cR_p=amp.e       
cL_p=amp.e   
ctmin=-1
ctmax=1
ct=np.array(CS.MC_int('s',N,E,mx,my,m1,m2,ctmin,ctmax,cR,cL,cR_p,cL_p)[1])
M=np.array(CS.MC_int('s',N,E,mx,my,m1,m2,ctmin,ctmax,cR,cL,cR_p,cL_p)[2])
theta=np.array(CS.MC_int('s',N,E,mx,my,m1,m2,ctmin,ctmax,cR,cL,cR_p,cL_p)[3])
for i in range(N):
    ct[i]=abs(ct[i])

dsigma=np.array(CS.MC_int('s',N,E,mx,my,m1,m2,ctmin,ctmax,cR,cL,cR_p,cL_p)[0])/(4*np.pi)
ref_dsigma=[]
for i in range(N):   
    ref_dsigma.append((2*(1/137)**2/((1-ct[i])**2))*(1+((1+ct[i])**2)/4))

# ref_dsigma=sum(ref_dsigma)/N
dsigma=sum(dsigma)/N
ref_ref_dsigma=[]
for i in range(N):
    ref_ref_dsigma.append(M[i]/(64*np.pi**2))
#print(ref_ref_dsigma/N)
# print(ref_dsigma)
plt.scatter(theta,ref_dsigma)
#confirm that dsigma/dcos is equal to the formula that I found. Produce a presentable graph





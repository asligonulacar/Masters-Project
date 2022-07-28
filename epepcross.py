#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 17 15:15:04 2022

@author: asligonulacar
"""

from FourVector import fv
from ME_squared import ME
from Generator import Gen
from Amplitude import amp
import numpy as np
from Channel_Basics import CB
import matplotlib.pyplot as plt
from histogram_test import histogram
import pylab
import cmath
from scipy import stats 
plt.style.use('seaborn-pastel')
E=1
mx=0
my=0.939
m1=0
m2=0.939#in GeV
ctmin=-1
ctmax=0.99
ch='s'
cR=cL=cR_p=cL_p=amp.e
 
def MC_int_Weighted(E,ch,mx,my,m1,m2,ctmin,ctmax,cR,cL,cR_p,cL_p):   
        err=0
        #sigma=np.array([])
        sigma2=[]
        #cross=np.array([])
        N=0
        ctl=[]
        RB=[]
        Q=[]
        #E3=[]
        #F2=[]
       # F1=[]
        # Gm=[]
        # Ge=[]
        Mref=[]
        # A=[]
        while err==0 or err>10000:#0.05*(sum(sigma2)/len(ctl)):
            N=N+5000
            for i in range(N):
                if ((i/N)*100)%10==0:
                    print((str((i/N)*100)) +'%')
                p=Gen.Lab(E, mx, my, m1, m2, ctmin, ctmax)
                ctc=p[4]
                if -1<ctc<0.99:
                    pm1=p[0]
                    pm2=p[1]
                    pm3=p[2]
                    pm4=p[3]
                    q=fv(fv.Sub(pm2,pm4))
                    Q2=abs(fv.Pro(q,q))
                    t=Q2/(4*my**2)
                    F1=(CB.Ff_Dipole(Q2, 1.79, my, 2.79)[0])
                    F2=(CB.Ff_Dipole(Q2, 1.79, my, 2.79)[1])
                    k=1.79
                    G_m=F1+k*F2
                    G_e=F1-t*F2
                    RBb=(((1/137)**2)/(4*pm1.x[0]**2*((1-ctc)/2)**2))*(pm1.x[0]/pm3.x[0])*(((G_e**2+t*G_m**2)/(1+t))*((1+ctc)/2)+ 2*t*G_m**2*((1-ctc)/2))
                    RB.append(RBb*10**2)
                    #Mm=((8*amp.e**4)/Q2**2)*((F1+k*F2)**2*((2*fv.Pro(p[0],p[1])*fv.Pro(p[2],p[1]))-fv.Pro(p[0],p[2])*(-fv.Pro(p[0],p[1])+fv.Pro(p[2],p[1])+my**2))-(((F1+k*F2)*F2)/2)*(4*(-my**2*fv.Pro(p[0],p[2])+2*fv.Pro(p[0],p[1])*fv.Pro(p[2],p[1])))+((k*F2)**2/(4*my**2))*(2*(fv.Pro(p[0],p[1])-fv.Pro(p[2],p[1])+2*my**2)*(-my**2*fv.Pro(p[0],p[2])+2*fv.Pro(p[0],p[1])*fv.Pro(p[2],p[1]))))                     
                    #Mref.append(Mm*10**-4)#((amp.e**4/(t**2))*((((G_e**2+t*G_m**2)/(1+t))*(lamb**2-t**2-t))+2*t**2*G_m**2))
                    ctl.append(ctc) 
                    Iso=Gen.Iso_Weight(p[0], p[1], p[2], p[3], ctmin, ctmax)
                    ff=Gen.Flux_CMS(p[0],p[1],m1,m2)
                    #Ps=Gen.PS(p[0], p[1], p[2], p[3],ctc)
                    a=amp.epAmp(pm1, pm2, pm3, pm4, cR, cL, cR_p, cL_p)[0]       
                    sigma2.append(abs(a*Iso*ff)*10**2/(4))#*p[2].x[0]/(p[0].x[0]**3))#              3.7*abs(Ps*a)
                else:
                     i=0
            
            err=np.var(sigma2)*0.9*10**-2/N**2
        plt.yscale('log')
        plt.ylabel('dσ/dΩ [10–30 cm2 /sr]')
        plt.xlabel('cos θ')
        plt.hist2d(ctl,sigma2)
        plt.colorbar()
        plt.plot(ctl,RB,'s',marker='.',markersize=2,color='black')
        plt.show()
        #plt.errorbar(ctl,sigma2,yerr=err,ls='none')
        sigma2_tot=sum(sigma2)/len(ctl)
        RB_tot=sum(RB)/len(ctl)
        return sigma2_tot/RB_tot

print(MC_int_Weighted(E,ch,mx,my,m1,m2,ctmin,ctmax,cR,cL,cR_p,cL_p))


def Energy(c,mx,my,m1,m2,ctmin,ctmax,cR,cL,cR_p,cL_p): 
        sigma=[]
        ct=[]
        RB=[]
        N=10
        for i in range(N):
            E=1
            sigma.append(MC_int_Weighted(E,ch,mx,my,m1,m2,ctmin,ctmax,cR,cL,cR_p,cL_p)[0])
            RB.append(MC_int_Weighted(E,ch,mx,my,m1,m2,ctmin,ctmax,cR,cL,cR_p,cL_p)[1])
            
        
       
        
        return None


    
#print(Energy(ch,mx,my,m1,m2,ctmin,ctmax,cR,cL,cR_p,cL_p))


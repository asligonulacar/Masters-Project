#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 21 18:28:21 2021

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
E=1
mx=0
my=0
m1=0
m2=0
ctmin=-1
ctmax=1
ch='s'
cR=cL=cR_p=cL_p=amp.e
def MC_int_Weighted(E,ch,mx,my,m1,m2,ctmin,ctmax,cR,cL,cR_p,cL_p):   
        err=0
        sigma=[]
        ct=[]
        cross=[]
        N=0
        #En=[]
        while err==0 or err>1e-9:
            N=N+100
            for i in range(N):
                # E=np.random.uniform(0,40*10**12)
                # En.append(E)
                p=Gen.Iso_Gen(E, mx, my, m1, m2, ctmin, ctmax)
                ct.append(abs(np.dot(p[1].y,p[2].y)/(np.linalg.norm(p[1].y)*np.linalg.norm(p[2].y))))
                Iso=Gen.Iso_Weight(p[0],p[1],p[2],p[3],ctmin,ctmax)
                f=(Gen.Flux(p[0], p[1], m1, m2))
                a=amp.Z_new(p[0], p[1], p[2], p[3], cR, cL, cR_p, cL_p, ctmin, ctmax,ch)
                M=ME.MESqr(a,ch,p[0],p[1],p[2],p[3])
                sigma.append(abs((M*Iso*f)))
                cross.append(((4*np.pi*((1/137)**2))/(3*E**2))*np.sqrt(1-(m1**2/(0.5)**2))*(1+(1/2)*(m1**2/(0.5)**2)))
            err=np.var(sigma)
            print(N)
            print(err)
        avg=np.mean(sigma)
        #err=np.std(sigma)/np.sqrt(N) 
        plt.hist2d(ct,sigma,bins=[20,20])
        plt.colorbar()
        plt.xlabel("| costheta |")
        plt.ylabel("Differential Cross Section")
        # plt.show()
        # plt.xlabel('E [MeV]' )
        # plt.ylabel('Cross Section [1/MeV^2]')
        # plt.scatter(En,sigma)
        # plt.xlabel('E [MeV]' )
        # plt.ylabel('Cross Section [1/MeV^2]')
        # plt.scatter(En,cross)
        
        
        sig_min=np.amin(sigma)
        sig_max=np.amax(sigma)
        cross=sum(cross)/N
        sigma_tot=sum(sigma)/N
        return cross/sigma_tot

print(MC_int_Weighted(E,ch,mx,my,m1,m2,ctmin,ctmax,cR,cL,cR_p,cL_p))
def MC_int_Unweighted(E,ch,mx,my,m1,m2,ctmin,ctmax,cR,cL,cR_p,cL_p):
        err=0
        sigma=[]
        sigma_tot=[]
        ct=[]
        N=10000
        p_acc=[]
        for i in range(N):
            p=Gen.Iso_Gen(E, mx, my, m1, m2, ctmin, ctmax)
            ct.append((np.dot(p[1].y,p[2].y)/(np.linalg.norm(p[1].y)*np.linalg.norm(p[2].y))))
            Iso=Gen.Iso_Weight(p[0],p[1],p[2],p[3],ctmin,ctmax)
            f=(Gen.Flux(p[0], p[1], m1, m2))
            a=amp.Z_new(p[0], p[1], p[2], p[3], cR, cL, cR_p, cL_p, ctmin, ctmax,ch)
            M=ME.MESqr(a,ch,p[0],p[1],p[2],p[3])
            sigma.append(abs((M*Iso*f)))            
        avg=np.mean(sigma)  
        sig_max=np.amax(sigma)
        for i in range(N):
            x=np.random.uniform(1,0)
            p=(sigma[i]/sig_max)
            if p<x:
                p_acc.append(p)
                sigma_tot.append(sigma[i]/N)
            else:
                p_acc.append(1-p)
                sigma_tot.append(0)
# =============================================================================
#         plt.hist2d(p_acc,sigma_tot,bins=[25,25])
#         plt.colorbar()
# =============================================================================
        err=np.var(p_acc)
        eff=np.average(p_acc)
        return err,eff
#print(MC_int_Unweighted(E,ch,mx,my,m1,m2,ctmin,ctmax,cR,cL,cR_p,cL_p))

#######Energy Dependence 
def Energy(ch,mx,my,m1,m2,ctmin,ctmax,cR,cL,cR_p,cL_p): 
        err=0
        sigma=[]
        ct=[]
        cross=[]
        N=0
        En=[]
        error=[]
        while err==0 or err>1e-7:
            N=N+100
            for i in range(N):
                E=np.random.uniform(1,40*10**12)
                En.append(E)
                p=Gen.Iso_Gen(E, mx, my, m1, m2, ctmin, ctmax)
                #ct.append(abs(np.dot(p[1].y,p[2].y)/(np.linalg.norm(p[1].y)*np.linalg.norm(p[2].y))))
                Iso=Gen.Iso_Weight(p[0],p[1],p[2],p[3],ctmin,ctmax)
                f=(Gen.Flux(p[0], p[1], m1, m2))
                a=amp.Z_new(p[0], p[1], p[2], p[3], cR, cL, cR_p, cL_p, ctmin, ctmax,ch)
                M=ME.MESqr(a,ch,p[0],p[1],p[2],p[3])
                sigma.append(abs((M*Iso*f)))
                cross.append(((4*np.pi*((1/137)**2))/(3*E**2))*np.sqrt(1-(m1**2/(E/2)**2))*(1+(1/2)*(m1**2/(E/2)**2)))
                error.append(np.var(sigma))
            err=np.var(sigma)
        #print(err)
        avg=np.mean(sigma)
        #err=np.std(sigma)/np.sqrt(N) 
        # plt.hist2d(ct,sigma,bins=[20,20])
        # plt.colorbar()
        #plt.xlabel("|)
        # plt.show()
        # plt.xlabel('E [MeV]' )
        # plt.ylabel('Cross Section [1/MeV^2]')
        # plt.scatter(En,sigma)
        # plt.xlabel('E [MeV]' )
        # plt.ylabel('Cross Section [1/MeV^2]')
        
        #plt.scatter(En,sigma)
        #plt.scatter(En,cross)
        sig_min=np.amin(sigma)
        sig_max=np.amax(sigma)
        cross=sum(cross)/N
        sigma_tot=sum(sigma)/N
        return cross/sigma_tot

    
#print(Energy(ch,mx,my,m1,m2,ctmin,ctmax,cR,cL,cR_p,cL_p))




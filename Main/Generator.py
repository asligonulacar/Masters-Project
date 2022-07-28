#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 11 15:00:59 2021

@author: asligonulacar
"""
from FourVector import fv  
from Channel_Basics import CB
import numpy as np
class Gen:
    def Iso_Gen(E,mx,my,m1,m2,ctmin,ctmax):
        #incoming momenta generation in the CMS frame
        p1x=(CB.SqLam(E**2,mx**2,my**2))**0.5/(2*E)
        #print(p1x)
        SphPolars=np.array([0,0,1])
        p1x_3=p1x*SphPolars
        p1y_3=-p1x_3
        Ex=(p1x**2+mx**2)**0.5
        Ey=(p1x**2+my**2)**0.5
        px=fv([Ex,p1x_3[0],p1x_3[1],p1x_3[2]])
        py=fv([Ey,p1y_3[0],p1y_3[1],p1y_3[2]])
        p = fv(fv.Add(px, py))
        p = fv(p.x)
        #outgoing momenta generation in the CMS frame 
        ran1 = np.random.uniform(0,1)
        ran2 = np.random.uniform(0,1)
        costheta=ctmin+(ctmax-ctmin)*ran1
        sintheta=np.sqrt(1.-costheta*costheta)
        phi=2.*np.pi*ran2
        new_p = fv.Boost_CMS(p, p)
        new_E = new_p.x[0]
        p1m = ((CB.SqLam(new_E**2, m1**2, m2**2))**0.5/(2*new_E))  
        #print(p1m)
        SphPolars = np.array([sintheta*np.cos(phi), sintheta*np.sin(phi), costheta])
        p1_3 = p1m*SphPolars  # three momentum
        p2_3 = -p1_3
        E1 = (p1m**2+m1**2)**0.5
        E2 = (p1m**2+m2**2)**0.5
        p1 = fv([E1, p1_3[0], p1_3[1], p1_3[2]])
        #print(p1.x)
        p2 = fv([E2, p2_3[0], p2_3[1], p2_3[2]])
        #print(p2.x)
        return px,py,p1,p2,costheta,sintheta,p1x,p1m
    def AnIso_Gen(E,mx,my,m1,m2,ctmin,ctmax,ctexp):
        #incoming momenta generation
        p1x=(CB.SqLam(E**2,mx**2,my**2))**0.5/(2*E)
        SphPolars=np.array([0,0,1])
        p1x_3=p1x*SphPolars
        p1y_3=-p1x_3
        Ex=(p1x**2+mx**2)**0.5
        Ey=(p1x**2+my**2)**0.5
        px=fv([Ex,p1x_3[0],p1x_3[1],p1x_3[2]],mx)
        py=fv([Ey,p1y_3[0],p1y_3[1],p1y_3[2]],my)
        p = fv(fv.Add(px, py))
        m = (p.x[0]**2-p.x[1]**2-p.x[2]**2-p.x[3]**2)**0.5  # calculates mass
        p = fv(p.x, m)
        #outgoing momenta generation
        ran1=np.random.uniform(0,1)
        ran2=np.random.uniform(0,1)
        s=CB.Abs2(p)
        pabs=np.sqrt(abs(s))
        e1=(s+m1-m2)/(pabs*2)
        p1m=pabs*(CB.SqLam(s,m1**2,m2**2))**0.5/2
        a=(p.x[0]*e1)/(m*p1m)
        if 0<=a<=1:
            a = 1.0000000001
        ct=CB.PeakedDist(a,ctexp,ctmin,ctmax,0,1,ran1)
        st=np.sqrt(1-ct*ct)
        phi=2*np.pi*ran2
        SphPolars=np.array([st*np.cos(phi), st*np.sin(phi), ct])
        p1_3 = p1m*SphPolars  # three momentum
        E1 = (p1m**2-m1**2)**0.5
        p1 = fv([E1, p1_3[0], p1_3[1], p1_3[2]], m1)
        #print(p1.x)
        p2 = fv(fv.Sub(p, p1))
        return px,py,p1,p2
    #this also feels wrong...
    def Iso_Weight(px,py,p1,p2,ctmin,ctmax):
        p = fv(fv.Add(px, py))
        pxc=fv.Boost_CMS(px, p)
        pyc=fv.Boost_CMS(py, p)
        s=(pxc.x[0]+pyc.x[0])**2
        p1c=fv.Boost_CMS(p1, p)
        p2c=fv.Boost_CMS(p2, p)
        r_2= CB.SqLam(s,CB.Abs2(p1c),CB.Abs2(p2c))**0.5*(ctmax-ctmin)/(16*np.pi*s)
        return r_2
        # cms_p=fv(fv.Add(px,py))
        # return 1/(4*(np.pi)*CB.SqLam(CB.Abs2(cms_p),CB.Abs2(p1),CB.Abs2(p2))**0.5*(ctmax-ctmin))
        
    def AnIso_Weight(px,py,p1,p2,m1,m2,ctmin,ctmax,ctexp):
        ran1=np.random.uniform(0,1)
        cms_p=fv(fv.Add(px,py))
        s=(px.x[0]+py.x[0])**2
        pabs=np.sqrt(abs(s))
        p1h_0=(s+m1-m2)/(2*pabs)
        a=(cms_p.x[0]*p1h_0)/(s*m1)
        if 0<=a<=1:
            a = 1.0000000001
        ct=(pabs*p1.x[0]-cms_p.x[0]*p1h_0)/(s*m1)
        if ct<ctmin or ct>ctmax:
            return 0
        wt=1/((np.pi*CB.SqLam(s,m1,m2)**0.5)/(4*(a+ct)**ctexp*CB.PeakedWeight(a,ctexp,ctmin,ctmax,ct,1,ran1)))    
        return wt
   #this is wrong somehow
    def Flux_CMS(px,py,m1,m2):
        p = fv(fv.Add(px, py))
        pxc=fv.Boost_CMS(px, p)
        pyc=fv.Boost_CMS(py, p)
        p1p2 = fv.Pro(pxc,pyc)
        #ff=8*np.pi**2*CB.SqLam((pxc.x[0]+pyc.x[0])**2, CB.Abs2(pxc), CB.Abs2(pyc))**0.5
        return  1/(4*((p1p2)**2-m1**2*m2**2)**0.5)
    
    def PS(p1,p2,p3,p4,ct):
        m=p2.m
        p1t=np.linalg.norm(p1.y)
        p3t=np.linalg.norm(p3.y)
        E1=p1.x[0]
        E3=p3.x[0]
        PS=((1/(64*np.pi**2*m*p1t))*(p3t**2/(((E1+m)*p3t)-p1t*E3*(ct))))
        #PS2=((1/(64*np.pi**2*m*p1t))*(p3t**2/(((E1+m)*p3t)-p1t*E3*(ct))))
        return PS
    
    def Flux_Lab(px,my):
        return 1/(4*my*np.sqrt(px.y[0]**2+px.y[1]**2+px.y[2]**2))
    
    def Lab(E,m1,m2,m3,m4,ctmin,ctmax):
        # p1x=(CB.SqLam(E**2,m1**2,m2**2))**0.5/(2*E)
        SphPolars = np.array([0,0,1])
        p1x_3=E*SphPolars
        p1=fv([E,p1x_3[0],p1x_3[1],p1x_3[2]])
        p2=fv([m2,0,0,0])
        p = fv(fv.Add(p1, p2))
        p1c=fv.Boost_CMS(p1, p)
        p2c=fv.Boost_CMS(p2, p)
        pcms=fv.Boost_CMS(p, p)
        ran1 = np.random.uniform(0,1)
        ran2 = np.random.uniform(0,1)
        costheta=ctmin+(ctmax-ctmin)*ran1
        sintheta=np.sqrt(1.-costheta*costheta)
        phi=2.*np.pi*ran2
        new_E = pcms.x[0]
        
        #s=new_E**2
        
        p3m = ((CB.SqLam(new_E**2, m3**2, m4**2))**0.5/(2*new_E))  
        SphPolars = np.array([sintheta*np.cos(phi),sintheta*np.sin(phi), costheta])
        p3_3 = p3m*SphPolars 
        p4_3 = -p3_3
        E3 = (p3m**2+m3**2)**0.5
        E4 = (p3m**2+m4**2)**0.5
        p3c = fv([E3, p3_3[0], p3_3[1], p3_3[2]])
        p4c = fv([E4, p4_3[0], p4_3[1], p4_3[2]])
        #t=fv.Pro(fv(fv.Sub(p1c,p3c)),fv(fv.Sub(p1c,p3c)))
        #u=fv.Pro(fv(fv.Sub(p1c,p4c)),fv(fv.Sub(p1c,p4c)))
        p3=fv.Boost_lab(p3c,p)
        p4=fv.Boost_lab(p4c,p)
        #s=fv.Pro(fv(fv.Add(p1,p2)),fv(fv.Add(p1,p2)))
        #t=fv.Pro(fv(fv.Sub(p1,p3)),fv(fv.Sub(p1,p3)))
        #u=fv.Pro(fv(fv.Sub(p1,p4)),fv(fv.Sub(p1,p4)))
        q=fv(fv.Sub(p1,p3))
        q2=(fv.Pro(q,q))
        ct=1-(fv.Pro(p1,p3)/(p1.x[0]*p3.x[0]))
        dq2=2*np.linalg.norm(p1c.y)*np.linalg.norm(p3c.y)*abs(costheta)
        
        #ct=((-2*m2**2*t)/((s-m2**2)*(u-m2**2)))+1
        #ctt=((s-m1**2-m2**2)*(m2**2+m3**2-u)+2*m2**2*(t-m1**2-m3**2))/((CB.SqLam(s,m1**2,m2**2))**0.5*(CB.SqLam(u,m2**2,m3**2))**0.5)
        return p1,p2,p3,p4,ct,dq2

# p1=fv([2,0,0,1])
# p2=fv([2,0,0,0])
# p=fv(fv.Add(p1,p2))
# p0=fv.Boost_CMS(p,p)
# print(p0.x,p0.m)
# p1c=fv.Boost_CMS(p1,p)
# print(p1c.x,p1c.m)
# p2c=fv.Boost_CMS(p2,p)
# print(p2c.x,p2c.m)
# print(p1c.x+p2c.x)
# p2l=fv.Boost_lab(p2c,p)
# p1l=fv.Boost_lab(p1c,p)
# print(p2l.x,p2l.m)
# print(p1l.x,p1l.m)

# E=1
# m1=m3=0
# m2=m4=0.71
# ctmin=-1
# ctmax=1
# print(Gen.Lab(E,m1,m2,m3,m4,ctmin,ctmax))
    
    
    
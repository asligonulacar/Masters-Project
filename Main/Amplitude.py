#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 10 18:09:31 2021

@author: asligonulacar
"""
import numpy as np
from Channel_Basics import CB
import cmath
from FourVector import fv
from Generator import Gen
class amp:  
    k0 = fv([1, 0, 1, 0])
    k1 = fv([0, 1, 0, 0])
    e = (4*np.pi/137)**0.5 
    def S_plus(p1, p2, k0, k1):
        X = (2*(fv.Pro(k0,p1)*fv.Pro(k1,p2)-fv.Pro(k0,p2)*fv.Pro(k1,p1)+complex(0, 1) *
            fv.eps_func(k0, k1, p1, p2)))/(fv.Eta(p1,k0)*fv.Eta(p2,k0))
                #np.einsum("ijkl,i,j,k,l->", fv.e_4, k0.x, k1.x, p1.x, p2.x)))/(fv.Eta(p1,k0)*fv.Eta(p2,k0))
        return X
    
    def S_minus(p1, p2, k0, k1):
        X = np.conj((2*(fv.Pro(k0,p2)*fv.Pro(k1,p1)-fv.Pro(k0,p1)*fv.Pro(k1,p2)+complex(0, 1) *
            fv.eps_func(k0, k1, p2, p1)))/(fv.Eta(p2,k0)*fv.Eta(p1,k0)))
                        #np.einsum("ijkl,i,j,k,l->", fv.e_4, k0.x, k1.x, p2.x, p1.x)))/(fv.Eta(p2,k0)*fv.Eta(p1,k0)))
        return X    
   
    def Helicity(A, initial_A):        
        initial_length = str(len(A))
        H = str(2**len(A)) 
        num1 = '0b1'  
        for index, value in enumerate(A):
            if value == -1:
                A[index] = 0
        A = bin(int(''.join(map(str, (A))), 2) << 1)  
        num = int(A, 2) 
        A = format(num, '0'+initial_length+'b')
        Helicities = [initial_A] 
        Binary = [A]
        for i in range(int(H)-1):
            A = bin(int(A, 2)+int(num1, 2))  
            num = int(A, 2) 
            A = format(num, '0'+initial_length+'b')  
            Binary.append(A)
            Y = [int(i)for i in A]  
            for index, value in enumerate(Y):
                if value == 0:
                    Y[index] = -1
            Helicities.append(Y)  
        return Binary
   
    def Z(px, py, p1, p2, cR, cL, cR_p, cL_p):
        A = [-1, -1, -1, -1]
        initial_A = [-1, -1, -1, -1]
        k0=amp.k0
        k1=amp.k1
        B = amp.Helicity(A,initial_A)
        E1 = fv.Eta(px,k0)  
        E2 = fv.Eta(py,k0)  
        E3 = fv.Eta(p1,k0)  
        E4 = fv.Eta(p2,k0)  
        Mu1 = px.m/E1#cmath.sqrt(fv.Pro(px,px))/fv.Eta(px,amp.k0)
        Mu2 = py.m/E2#cmath.sqrt(fv.Pro(py,py))/fv.Eta(py,amp.k0)
        Mu3 = p1.m/E3#cmath.sqrt(fv.Pro(p1,p1))/fv.Eta(p1,amp.k0)
        Mu4 = p2.m/E4#cmath.sqrt(fv.Pro(p2,p2))/fv.Eta(p2,amp.k0)
        AA=np.array([])
        for i in B:
            B[-1] = -2*(amp.S_plus(p1, px, k0, k1)*amp.S_minus(p2, py, k0, k1)*cR_p*cR-Mu1*Mu2*E3*E4*cR_p*cL-E1*E2*Mu3*Mu4*cL_p*cR)
            B[-2] = -2*E2*cR*(amp.S_plus(p2, px, k0, k1)*Mu3*cL_p-amp.S_plus(p1, px, k0, k1)*Mu4*cR_p) #changed second Sminus to Splus
            B[-3] = -2*E1*cR*(amp.S_minus(py, p1, k0, k1)*Mu4*cL_p-amp.S_minus(py, p2, k0, k1)*Mu3*cR_p)
            B[-4] = -2*(amp.S_plus(px, p2, k0, k1)*amp.S_minus(py, p1, k0, k1)*cL_p*cR-Mu1*Mu2*E3*E4*cL_p*cL-E1*E2*Mu3*Mu4*cR_p*cR)
            B[-5] = -2*E4*cR_p*(amp.S_plus(p1, px, k0, k1)*Mu2*cR-amp.S_plus(p1, py, k0, k1)*Mu1*cL) #changed first Sminus to Splus
            B[-6] = 0
            B[-7] = -2*(Mu1*Mu4*E2*E3*cL_p*cL+Mu2*Mu3*E1*E4*cR_p*cR -
                        Mu2*Mu4*E1*E3*cL_p*cR-Mu1*Mu3*E2*E4*cR_p*cL)
            B[-8] = -2*E3*cL_p*(amp.S_plus(py, p2, k0, k1)*Mu1*cL- amp.S_plus(px, p2, k0, k1)*Mu2*cR)
            B[-9] =  -2*E3*cR_p*(amp.S_minus(py, p2, k0, k1)*Mu1*cR-amp.S_minus(px, p2, k0, k1)*Mu2*cL) #changed all Splus to Sminus
            B[-10] = -2*(Mu1*Mu4*E2*E3*cR_p*cR+Mu2*Mu3*E1*E4*cL_p *
                         cL-Mu2*Mu4*E1*E3*cR_p*cL-Mu1*Mu3*E2*E4*cL_p*cR)
            B[-11] = 0
            B[-12] = -2*E4*cL_p*(amp.S_minus(p1, px, k0, k1)*Mu2*cL-amp.S_minus(p1, py, k0, k1)*Mu1*cR)
            B[-13] = -2*(amp.S_minus(px, p2, k0, k1)*amp.S_plus(py, p1, k0, k1)*cR_p*cL-Mu1*Mu2*E3*E4 *
                         cR_p*cR-E1*E2*Mu3*Mu4*cL_p*cL)
            B[-14] = -2*E1*cL*(amp.S_plus(py, p1, k0, k1)*Mu4*cR_p-amp.S_plus(py, p2, k0, k1)*Mu3*cL_p)
            B[-15] = -2*E2*cL*(amp.S_minus(p2, px, k0, k1)*Mu3*cR_p-amp.S_minus(p1, px, k0, k1)*Mu4*cL_p)
            B[-16] = -2*(amp.S_minus(p1, px, k0, k1)*amp.S_plus(p2, py, k0, k1)*cL_p*cL-Mu1*Mu2*E3*E4 *
                         cL_p*cR-E1*E2*Mu3*Mu4*cR_p*cL)
        AA=np.append(AA,B)
        # for i in range(0,16):
        #     AA[i]=AA[i]*np.conj(AA[i])
        return AA

    def Y(px, py, cR, cL):
        k0=amp.k0
        k1=amp.k1
        A=[-1,-1]
        initial_A=[-1,-1]
        B=amp.Helicity(A,initial_A)
        E1 = fv.Eta(px,k0)  
        E2 = fv.Eta(py,k0)  
        Mu1 = px.m/E1#cmath.sqrt(fv.Pro(px,px))/fv.Eta(px,amp.k0)/E1
        Mu2 = py.m/E2#cmath.sqrt(fv.Pro(py,py))/fv.Eta(py,amp.k0)/E2
        AA=np.array([])
        for i in B:
            B[0]=cR*Mu1*E2+cL*Mu2*E1
            B[1]=cL*amp.S_minus(px,py,k0,k1)
            B[2]=(cR*amp.S_plus(px,py,k0,k1))
            B[3]=(cL*Mu1*E2+cR*Mu2*E1)
        AA=np.append(AA,B)
        return AA
    
    def X(px, py, p1, cR, cL):
        k0=amp.k0
        k1=amp.k1
        A=[-1,-1]
        initial_A=[-1,-1]
        B=amp.Helicity(A,initial_A)
        E1 = fv.Eta(px,k0)
        E2 = fv.Eta(py,k0)
        E3 = fv.Eta(p1,k0) 
        Mu1 = px.m/E1#cmath.sqrt(fv.Pro(px,px))/fv.Eta(px,amp.k0)
        Mu2 = py.m/E2#cmath.sqrt(fv.Pro(py,py))/fv.Eta(py,amp.k0)
        Mu3 = p1.m/E3#cmath.sqrt(fv.Pro(p1,p1))/fv.Eta(p1,amp.k0)
        AA=np.array([])
        for i in B:
            B[-1]=cL*E2**2*Mu1*Mu3+cR*Mu2**2*E1*E3+cR*amp.S_plus(px,py,k0,k1)*amp.S_minus(py,p1,k0,k1) 
            B[-2]=E2*(cR*Mu1*amp.S_plus(py,p1,k0,k1)+cL*Mu3*amp.S_plus(px,py,k0,k1))
            B[-3]=E2*(cL*Mu1*amp.S_minus(py,p1,k0,k1)+cR*Mu3*amp.S_minus(px,py,k0,k1))
            B[-4]=(cR*E2**2*Mu1*Mu3+cL*Mu2**2*E1*E3+cL*amp.S_minus(px,py,k0,k1)*amp.S_plus(py,p1,k0,k1))
        AA=np.append(AA,B)
        return AA
    
    def epAmp(pe1, pp2, pe3, pp4, cR, cL, cR_p, cL_p):
        #with photon exchange
        q=fv(fv.Sub(pp2,pp4))
        ptilda=fv(fv.Add(pp4,pp2))
        Q=abs(fv.Pro(q,q))
        
        M=pp2.m#cmath.sqrt(fv.Pro(pp2,pp2))/fv.Eta(pp2,amp.k0)
        k=1.79
        # G_m=F1+k*F2#1/(1+(Q/0.71))**2
        # G_e=F1+(k*Q/(4*M**2))*F2#G_m/2.79
        F1=CB.Ff_Dipole(Q, k, M, 2.79)[0]
        F2=(CB.Ff_Dipole(Q, k, M, 2.79)[1])
        #f=(F1+k*F2)
        #g=-k*F2/(2*M) #(F1-((k*F2*Q/(4*M**2))))/Q 
        Z=amp.e**2*(F1+F2)*np.sum(amp.Z(pe3, pe1, pp4, pp2, 1, 1,1, 1))/(2*Q)
        XY=amp.e**2*(-F2)*np.sum(amp.X(pe3,ptilda,pe1,1,1))*(np.sum(amp.Y(pp4,pp2,1,1)))/(2*Q*M)
        
        #C=-f*np.sum(amp.X(pe3,q,pe1,cR,cL))*np.sum(amp.X(pp4,q,pp2,cR,cL))/(4*M**2)
        #D=-g*(fv.Pro(pp2,q)+fv.Pro(pp4,q))*np.sum(amp.X(pe3,q,pe1,cR,cL))*np.sum(amp.Y(pp4,pp2,cR,cL))/(4*M**2)
        
        ff=((Z+XY))
        ME2=ff*np.conj(ff)/(4)
        return ME2,F1,F2
    
    def Inv_Mu_amp(pvm, pe, pm, pve, cR, cL, cR_p, cL_p):
        Mw=80
        gw=0.653
        Z=np.sum((amp.Z(pvm, pe, pm, pve, 0, 2*gw/np.sqrt(8), 0, 2*gw/np.sqrt(8))))
        f=Z*((1/Mw)**2)
        ME=f*np.conj(f)/2
        return ME
    
    # def Compton_amp(pg1,pe2,pg3,pe4,cR,cL,cR_p,cL_p):
    #     s=fv.Pro(fv(fv.Add(pg1,pe2)),fv(fv.Add(pg1,pe2)))
    #     Z=np.sum(amp.Z(pg1, pe2, pg3, pe4, cR, cL,cR_p,cL_p))
    #     f=Z/s
    #     ME=f*np.conj(f)/4
    #     return ME
   
    def Electron_Neutrino_amp(pv1,pe2,pv3,pe4,cR,cL,cR_p,cL_p):
        Mw=80
        Z=np.sum((amp.Z(pv1, pe2, pv3, pe4, 0, cL, 0, cL_p)))
        f=Z*((1/Mw)**2)
        ME=f*np.conj(f)/8
        return ME
    
    def E_mu_amp(pe1,pe2,pm3,pm4):
        Mz=80
        G=np.sqrt(2)*0.653**2/(8*Mz**2)
        s=fv.Pro(fv(fv.Add(pe1,pe2)),fv(fv.Add(pe1,pe2)))
        Zg=np.sum(amp.Z(pe1, pe2, pm3, pm4, amp.e, amp.e, amp.e, amp.e))/s
        Zz=-np.sqrt(2)*G*Mz**2*(np.sum(amp.Z(pm3, pm4, pe1, pe2, 0.5, -0.5, 0.5, -0.5)))/(s-Mz**2)
        f=(Zg+Zz)
        ME=f*np.conj(f)/4
        return ME
    
    def Quasi_neutr_amp(pv,pn,pe,pp):
        q=fv(fv.Sub(pv,pe))
        Q2=abs(fv.Pro(q,q))
        ptilda=fv(fv.Add(pn,pp))
        
        Mw=80
        M2=pn.m
        M3=pp.m
        M=(M2+M3)/2
        
        #aw=0.653**2/(np.pi*4)
        Gf=0.653**2/(8*Mw**2)
        cc=0.975 #cos of cabibbo angle 
        #form factors
        F1=CB.FDipole(Q2,M)[0]
        F2=CB.FDipole(Q2,M)[1]
        FA=CB.FDipole(Q2,M)[2]
        
        Fp=0#(2*M**2/(0.134**2+Q2))*FA
        Z=np.sum(amp.Z((pe), (pv), pp, pn, 0 , 1, F1+F2+FA,F1+F2-FA))
        X1=np.sum(amp.X((pe), q, (pv), 0, 1))
        X2=np.sum(amp.X((pe), ptilda, (pv), 0 , 1))
        Y1=np.sum(amp.Y(pp, pn, -F2/(2*M),  -F2/(2*M)))
        Y2=np.sum(amp.Y(pp, pn, Fp/M, -Fp/M))
        f=(Z+Y1*X2+Y2*X1)*cc*Gf
        ME=f*np.conj(f)/(8*np.pi)
        
        # su=4*M*pv.x[0]-Q2
        # A=(Q2/(4*M**2))*((4+(Q2/(M**2)))*FA**2-(4+(-Q2/(M**2)))*F1**2+(Q2/M**2)*(1-(Q2/(4*M**2)))*F2**2+((4*Q2*F1*F2)/(M**2)))
        # B=(Q2/(M**2))*FA*(F1+F2)
        # C=(FA**2+F1**2+(Q2/(4*M**2))*F2**2)/4
        
        # Mref=Gf**2*M**4*cc**2*(A-(B*su/M**2)+(C*(su)**2/(M**4)))
        # Mref2=Gf**2*M**4*cc**2*(A+(B*su/M**2)+(C*(su)**2/(M**4)))
        #print(A,B,C)
        s=fv.Pro(fv(fv.Add(pv,pn)),fv(fv.Add(pv,pn)))
        t=fv.Pro(fv(fv.Sub(pv,pe)),fv(fv.Sub(pv,pe)))
        N=F1**2*(0.5*(2*(M**2-s)**2+2*s*t+t**2+pe.m**2*(-2*s-t))) + (F2**2/(4*M**2))*(0.25*(-2*t*(M**4-2*s*(2*M**2)+M**4+2*s**2)+2*t**2*(4*M**2-2*s)+pe.m**2*(t*(-4*M**2+4*s)+t**2)+pe.m**4*(-(4*M**2+t)))) + FA**2*(0.5*(2*(M**2-s)**2-t*((2*M)**2-2*s)+t**2+pe.m**2*(4*M**2--2*s-t))) + (Fp**2/(4*M**2))*(pe.m**2*(pe.m**2-t)*(-t))  + (F1*F2/(2*M))*(-(-t*2*M*t+pe.m**2*(M*t)+pe.m**4*M))   - F1*FA*(-(t*(2*M**2-2*s-t)+pe.m**2*(t))) - (F2*FA/(2*M))*(-(2*M)*(t*(2*M**2-2*s-t)+pe.m**2*(t)))+ (FA*Fp/(2*M))*(-2*pe.m**2*(pe.m**2*M-M*(s+t)+M*s))           
        Mref2=Gf**2*cc**2*N       
        return ME, Mref2
    
    def Quasi_antineutr_amp(pv,pn,pe,pp):
        q=fv(fv.Sub(pv,pe))
        Q2=abs(fv.Pro(q,q))
        ptilda=fv(fv.Add(pn,pp))
        
        Mw=80
        M2=pn.m
        M3=pp.m
        M=(M2+M3)/2
        
        #aw=0.653**2/(np.pi*4)
        Gf=0.653**2/(8*Mw**2)
        cc=0.975 #cos of cabibbo angle 
        #form factors
        F1=CB.FDipole(Q2,M)[0]
        F2=CB.FDipole(Q2,M)[1]
        FA=CB.FDipole(Q2,M)[2]
        
        Fp=0#(2*M**2/(0.134**2+Q2))*FA
        # Z=np.sum(amp.Z(fv.anti(pe), fv.anti(pv), pp, pn, 1 , 0, F1+F2+FA,F1+F2-FA))
        # X1=np.sum(amp.X(fv.anti(pe), q, fv.anti(pv), 1, 0))
        # X2=np.sum(amp.X(fv.anti(pe), ptilda, fv.anti(pv), 0 , 1))
        # Y1=np.sum(amp.Y(pp, pn, -F2/(2*M),  -F2/(2*M)))
        # Y2=np.sum(amp.Y(pp, pn, Fp/M, -Fp/M))
        # f=(Z+Y1*X2+Y2*X1)*cc*Gf
        # ME=f*np.conj(f)/(4*np.pi)
        
        su=4*M*pv.x[0]-Q2
        A=(Q2/(4*M**2))*((4+(Q2/(M**2)))*FA**2-(4+(-Q2/(M**2)))*F1**2+(Q2/M**2)*(1-(Q2/(4*M**2)))*F2**2+((4*Q2*F1*F2)/(M**2)))
        B=(Q2/(M**2))*FA*(F1+F2)
        C=(FA**2+F1**2+(Q2/(4*M**2))*F2**2)/4
        
        
        Mref2=Gf**2*M**4*cc**2*(A+(B*su/M**2)+(C*(su)**2/(M**4)))
        #print(A,B,C)
        # s=fv.Pro(fv(fv.Add(pv,pn)),fv(fv.Add(pv,pn)))
        # t=fv.Pro(fv(fv.Sub(pv,pe)),fv(fv.Sub(pv,pe)))
        # N=F1**2*(0.5*(2*(M**2-s)**2+2*s*t+t**2+pe.m**2*(-2*s-t))) + (F2**2/(4*M**2))*(0.25*(-2*t*(M**4-2*s*(2*M**2)+M**4+2*s**2)+2*t**2*(4*M**2-2*s)+pe.m**2*(t*(-4*M**2+4*s)+t**2)+pe.m**4*(-(4*M**2+t)))) + FA**2*(0.5*(2*(M**2-s)**2-t*((2*M)**2-2*s)+t**2+pe.m**2*(4*M**2--2*s-t))) + (Fp**2/(4*M**2))*(pe.m**2*(pe.m**2-t)*(-t))  + (F1*F2/(2*M))*(-(-t*2*M*t+pe.m**2*(M*t)+pe.m**4*M))   - F1*FA*(-(t*(2*M**2-2*s-t)+pe.m**2*(t))) - (F2*FA/(2*M))*(-(2*M)*(t*(2*M**2-2*s-t)+pe.m**2*(t)))+ (FA*Fp/(2*M))*(-2*pe.m**2*(pe.m**2*M-M*(s+t)+M*s))           
        # Mref2=Gf**2*cc**2*N       
        return Mref2
    
   



# cR=cR_p=amp.e
# cL=cL_p=amp.e
# m1=0
# m2=0.5
# m3=0
# m4=0.5
# ctmin=-1
# ctmax=1
# E=1

# p1=Gen.Lab(E,m1,m2,m3,m4,ctmin,ctmax)[0]
# p2=Gen.Lab(E,m1,m2,m3,m4,ctmin,ctmax)[1]
# p3=Gen.Lab(E,m1,m2,m3,m4,ctmin,ctmax)[2]
# p4=Gen.Lab(E,m1,m2,m3,m4,ctmin,ctmax)[3]

# #calculating this cross section in the lab frame
# ME=amp.epAmp(p1, p2, p3, p4, cR, cL, cR_p, cL_p, ctmin, ctmax)[0]
# W=Gen.Iso_Weight(p1, p2, p3, p4, ctmin, ctmax)
# F=Gen.Flux(p1, p2, m3, m4)
# Cross=(ME*W*F)

# #analytical
# G_m=amp.epAmp(p1, p2, p3, p4, cR, cL, cR_p, cL_p, ctmin, ctmax)[1]
# G_e=amp.epAmp(p1, p2, p3, p4, cR, cL, cR_p, cL_p, ctmin, ctmax)[2]
# t=-fv.Pro(fv(fv.Sub(p1, p3)),fv(fv.Sub(p1, p3)))#/(4*m4**2)
# ct=Gen.Lab(E,m1,m2,m3,m4,ctmin,ctmax)[4]
# RB=((1/137)**2/(4*p1.x[0]**2*((1-ct)/2)**2))*(p3.x[0]/p1.x[0])*(((G_e**2+t*G_m**2)/(1+t))*((1+ct)/2)+2*t*G_m**2*((1-ct)/2))






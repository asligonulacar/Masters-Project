#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 10 18:09:31 2021

@author: asligonulacar
"""
import numpy as np
from FourVector import fv
class amp:  
    k1 = fv([0, 1, 0, 0])
    k0 = fv([1, 0, 1, 0])
    e = (4*np.pi/137)**0.5 
    A = [-1, -1, -1, -1]
    initial_A = [-1, -1, -1, -1]
    def S_plus(p1, p2, k0, k1):
        X = (2*(fv.Pro(p1,k0)*fv.Pro(p2,k1)-fv.Pro(p1,k1)*fv.Pro(p2,k0)+complex(0, 1) *
            np.einsum("ijkl,i,j,k,l->", fv.e_4, k0.x, k1.x, p1.x, p2.x)))/(((2*fv.Pro(p1,k0))**0.5)*((2*fv.Pro(p2,k0))**0.5))
        return X
    
    def S_minus(p1, p2, k0, k1):
        X = np.conj((2*(fv.Pro(p2,k0)*fv.Pro(p1,k1)-fv.Pro(p2,k1)*fv.Pro(p1,k0)+complex(0, 1) *
            np.einsum("ijkl,i,j,k,l->", fv.e_4, k0.x, k1.x, p2.x, p1.x)))/(((2*fv.Pro(p2,k0))**0.5)*((2*fv.Pro(p1,k0))**0.5)))
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
   
    def Z_new(px, py, p1, p2, cR, cL, cR_p, cL_p, ctmin, ctmax,ch):
        k0=amp.k0
        k1=amp.k1
        B = amp.Helicity(amp.A, amp.initial_A)
        E1 = fv.Eta(px,k0)  
        E2 = fv.Eta(py,k0)  
        E3 = fv.Eta(p1,k0)  
        E4 = fv.Eta(p2,k0)  
        Mu1 = px.m/E1
        Mu2 = -py.m/E2
        Mu3 = p1.m/E3
        Mu4 = -p2.m/E4
        for i in B:
            B[-1] = -2*(amp.S_plus(p1, px, k0, k1)*amp.S_minus(p2, py, k0, k1)*cR_p*cR-Mu1*Mu2*E3*E4*cR_p*cL-E1*E2*Mu3*Mu4*cL_p*cR)
            B[-2] = -2*E2*cR*(amp.S_plus(p2, px, k0, k1)*Mu3*cL_p-amp.S_minus(p2, py, k0, k1)*Mu4*cR_p)
            B[-3] = -2*E1*cR*(amp.S_minus(py, p1, k0, k1)*Mu4*cL_p-amp.S_minus(py, p2, k0, k1)*Mu3*cR_p)
            B[-4] = -2*(amp.S_plus(px, p2, k0, k1)*amp.S_minus(py, p1, k0, k1)*cL_p*cR-Mu1*Mu2*E3*E4*cL_p*cL-E1*E2*Mu3*Mu4*cR_p*cR)
            B[-5] = -2*Mu4*cR_p*(amp.S_minus(p2, py, k0, k1)*Mu2*cR-amp.S_plus(p1, py, k0, k1)*Mu1*cL)
            B[-6] = 0
            B[-7] = -2*(Mu1*Mu4*E2*E3*cL_p*cL+Mu2*Mu3*E1*E4*cR_p*cR -
                        Mu2*Mu4*E1*E3*cL_p*cR-Mu1*Mu3*E2*E4*cR_p*cL)
            B[-8] = -2*E3*cL_p*(amp.S_plus(py, p2, k0, k1)*Mu1*cL- amp.S_plus(px, p2, k0, k1)*Mu2*cR)
            B[-9] = -2*E3*cR_p*(amp.S_plus(py, p2, k0, k1)*Mu1*cR-amp.S_plus(px, p2, k0, k1)*Mu2*cL)
            B[-10] = -2*(Mu1*Mu4*E2*E3*cR_p*cR+Mu2*Mu3*E1*E4*cL_p *
                         cL-Mu2*Mu4*E1*E3*cR_p*cL-Mu1*Mu3*E2*E4*cL_p*cR)
            B[-11] = 0
            B[-12] = -2*Mu4*cL_p*(amp.S_minus(p1, px, k0, k1)*Mu2*cL-amp.S_minus(p1, py, k0, k1)*Mu1*cR)
            B[-13] = -2*(amp.S_minus(px, p2, k0, k1)*amp.S_plus(py, p1, k0, k1)*cR_p*cL-Mu1*Mu2*E3*E4 *
                         cR_p*cR-E1*E2*Mu3*Mu4*cL_p*cL)
            B[-14] = -2*E1*cL*(amp.S_plus(py, p1, k0, k1)*Mu4*cR_p-amp.S_plus(py, p2, k0, k1)*Mu3*cL_p)
            B[-15] = -2*E2*cL*(amp.S_minus(p2, px, k0, k1)*Mu3*cR_p-amp.S_minus(p1, px, k0, k1)*Mu4*cL_p)
            B[-16] = -2*(amp.S_minus(p1, px, k0, k1)*amp.S_plus(p2, py, k0, k1)*cL_p*cL-Mu1*Mu2*E3*E4 *
                         cL_p*cR-E1*E2*Mu3*Mu4*cR_p*cL)
        A=[]
        for i in range(0,16):
            A.append(B[i]*np.conj(B[i]))   
        return sum(A)
         

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 21 21:18:38 2021

@author: asligonulacar
"""

import numpy as np
import FourVector as fv
class CB:
     def SqLam(E, m1, m2):  # SqLam function from Channel Basics
        arg1 = ((E-m1-m2)**2-4*m1*m2)
        if arg1 >= 0:
            return arg1
   
     def Abs2(self):
        return (self.x[0]**2-self.x[1]**2-self.x[2]**2-self.x[3]**2)
    
     def PeakedDist(a,cn,cxm,cxp,res,k,ran):
         ce=1-cn
         res=0
         if ce==0:
             res=k*((a+k*cxm)*pow((a+k*cxp)/(a+k*cxm),ran)-a)
         else:
             res=k*(pow(ran*pow(a+k*cxp,ce)+(1.-ran)*pow(a+k*cxm,ce),1/ce)-a)
         return res
    
     def PeakedWeight(a,cn,cxm,cxp,res,k,ran):         
         ce=1-cn
         if ce>0 or ce<0:
            amin=(a+k*cxm)**ce
            wt=(a+k*cxp)**ce-amin
            ran=((a+k*res)**ce-amin)/wt
            wt/=k*ce
         else:
            amin=(a+k*cxm)
            wt=np.log((a+k*cxp)/amin)
            ran =np.log((a+k*res)/amin)/wt;
            wt /= k
         return wt
    
     def Ff_Dipole(Q2,k,M,mu):
         #simple
         F1=((Q2/(4*M**2))+mu)*(1/(1+(Q2/0.71)))**(2)*(1/(1+(Q2/(4*M**2))))
         F2=(((1+(Q2/0.71))**(-2))-F1)/k
         return F1,F2
     
     def FDipole(Q2,M):
         Gev=(1+(Q2/0.71))**-2
         Gem=4.71*((1+(Q2/0.71))**-2)
         Fv1=(Gev+((Q2/(4*M**2)))*Gem)/(1+(Q2/(4*M**2)))
         cFv2=(Gem-Gev)/(1+(Q2/(4*M**2)))
         Fa=-1.23*(1+(Q2/1.79))**-2
         return Fv1,cFv2,Fa
        
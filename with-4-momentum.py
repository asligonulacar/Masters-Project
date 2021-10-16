#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 20 13:32:23 2021

@author: asligonulacar
"""

import numpy as np
import random

class FourVector:
    
    g=np.array([[1,0,0,0], #metric tensor
               [0,-1,0,0],
               [0,0,-1,0],
               [0,0,0,-1]])
    
    c=1 #natural units
    
    def __init__(self,X=None,m=0,gam=0): 
        
        self.x=np.array([X[0],X[1],X[2],X[3]]) #4 vector
        
        self.m=m
        
        self.y=np.array([X[1],X[2],X[3]])
     
    #def LorentzTransform(self): 
        #general Lorentz transform matrix for 
        
    def Pro(self,other): #as the name says...
        return self@FourVector.g@other  
        #return np.trace(self*FourVector.g*other)
    
    def Addition(self,other): #vector addition
        return np.array([self.x[0]+other.x[0],self.x[1]+other.x[1],self.x[2]+other.x[2],self.x[3]+other.x[3]])  

    def Subtraction(self,other): #vector subtraction because why not
        return np.array([self.x[0]-other.x[0],self.x[1]-other.x[1],self.x[2]-other.x[2],self.x[3]-other.x[3]])  
    
    def EtaForZ(self,other):
        return (2*FourVector.FourDotProduct(self,other))**0.5


# def NewP1(new_E,m1,m2,SqLam):
#     phi=random.uniform(0,1)*2*np.pi
#     n=np.array([0,0,1,1])
#     l=np.array([0,1,0,0])
#     ran1=random.uniform(0,np.pi)
#     p1m=(SqLam(new_E**2,m1**2,m2**2))**0.5/(2*new_E) #from the notes, this should give abs(3-momentum) squared
#     E1=(new_E/2)*(1+(m1**2/new_E**2)-(m2**2/new_E**2))
#     v=np.array([E1,p1m*ran1*new_p.y[0],p1m*ran1*new_p.y[1],p1m*ran1*new_p.y[2]])
#     v+=p1m*(np.cos(phi)*n+np.sin(phi)*l)
#     return FourVector(v,m1)
            

#----------------------------_MOMENTA_---------------------------------------------------------
p=FourVector([0.7,0.1,0.2,0.3]) #input four monenta here
m_beta=p
norm=(p.x[0]**2+p.x[1]**2+p.x[2]**2+p.x[3]**2)**0.5

m=(p.x[0]**2-p.x[1]**2-p.x[2]**2-p.x[3]**2)**0.5 #calculates mass 
#print(m)
p=FourVector(p.x,m)


def SqLam(E,m1,m2): #SqLam function from Channel Basics
  arg1=((E-m1-m2)**2-4*m1*m2)
  if arg1>=0:
      return arg1 

def Boost(p):
    p0=(p.x[0]*m_beta.x[0]-np.dot(m_beta.y,p.y))/p.m
    c1=(p.x[0]+p0)/(p.m+m_beta.x[0])
    c2=p.y-c1*m_beta.y
    p=FourVector([p0,c2[0],c2[1],c2[2]],p0)
    return p

new_p=Boost(p)


m1=0.1
m2=0.1
new_E=new_p.x[0]
#print(new_p.x)
def NewP1(new_E,m1,m2,SqLam):
    #CMS frame
    p1m=(SqLam(new_E**2,m1**2,m2**2))**0.5/(2*new_E) 
    theta=random.uniform(0,2)*np.pi
    phi=random.uniform(0,1)*np.pi
    SphPolars=np.array([np.sin(theta)*np.cos(phi),np.sin(theta)*np.sin(phi),np.cos(theta)])
    p1_3=p1m*SphPolars #three momentum
    E1=(np.dot(p1_3,p1_3)+m1**2)**0.5
    p1=FourVector([E1,p1_3[0],p1_3[1],p1_3[2]])
    #LAB frame
    # if m1 and m2>0:
    #     ELAB=(s-m**2-m**2)/(2*m)
    #     print(ELAB)
    #     pLAB=(ELAB**2-m**2)**0.5
        
    #     p=FourVector([ELAB,pLAB,0,0])
    # else:
    #     return 0
    return p1

p1=NewP1(new_E,m1,m2,SqLam) #p1 calculated in the CMS frame.

#Calcaulating p2 in the CMS frame.
p2_3=-p1.y
E2=(np.dot(p2_3,p2_3)+m2**2)**0.5
p2=FourVector([E2,p2_3[0],p2_3[1],p2_3[2]])
#print(new_p.x)
print(p1.x)

print(p2.x)

print(p1.x+p2.x)

#alpha is the fine structure constant 
#e is the strength of the coupling 
#costheta max is deflection angle. costheta max closer to one  the more the cross section explodes
#costheta min =-1 the electron goes back to where it comes from.



#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 10 18:31:45 2021

@author: asligonulacar
"""

import numpy as np
import cmath
class fv:
    #metric tensor
    g = np.array([[1, 0, 0, 0],  
                  [0, -1, 0, 0],
                  [0, 0, -1, 0],
                  [0, 0, 0, -1]])
    #natural units
    c = 1  
    #epsilon four
    e_4 = ([[[[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]],
                  [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 1], [0, 0, -1, 0]],
                  [[0, 0, 0, 0], [0, 0, 0, -1], [0, 0, 0, 0], [0, 1, 0, 0]],
                  [[0, 0, 0, 0], [0, 0, 1, 0], [0, -1, 0, 0], [0, 0, 0, 0]]],
                 [[[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, -1], [0, 0, 1, 0]],
                  [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]],
                  [[0, 0, 0, 1], [0, 0, 0, 0], [0, 0, 0, 0], [-1, 0, 0, 0]],
                  [[0, 0, -1, 0], [0, 0, 0, 0], [1, 0, 0, 0], [0, 0, 0, 0]]],
                 [[[0, 0, 0, 0], [0, 0, 0, 1], [0, 0, 0, 0], [0, -1, 0, 0]],
                  [[0, 0, 0, -1], [0, 0, 0, 0], [0, 0, 0, 0], [1, 0, 0, 0]],
                  [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]],
                  [[0, 1, 0, 0], [-1, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]],
                 [[[0, 0, 0, 0], [0, 0, -1, 0], [0, 1, 0, 0], [0, 0, 0, 0]],
                  [[0, 0, 1, 0], [0, 0, 0, 0], [-1, 0, 0, 0], [0, 0, 0, 0]],
                  [[0, -1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]],
                  [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]]])

    
    def __init__(self, X=None,m=0):
        self.x = np.array([X[0], X[1], X[2], X[3]])  #4 vector
        self.y = np.array([X[1], X[2], X[3]]) #spatial component    
        self.m = cmath.sqrt((X[0]**2-X[1]**2-X[2]**2-X[3]**2)) #mass
        
                      
    def Pro(self, other):  #four vector product
        return self.x@fv.g@other.x

    def Add(self, other):  #four vector addition
        return np.array([self.x[0]+other.x[0], self.x[1]+other.x[1], self.x[2]+other.x[2], self.x[3]+other.x[3]])

    def Sub(self, other):  #four vector subtraction
        return np.array([self.x[0]-other.x[0], self.x[1]-other.x[1], self.x[2]-other.x[2], self.x[3]-other.x[3]])

    def Mul(self,other):
        return fv([other*self.x[0], other*self.x[1],other*self.x[2],other*self.x[3]])
    
    def Eta(self, other): #for the z functions
        return cmath.sqrt((2*fv.Pro(self, other)))
    
    def Boost_CMS(plab,p):
        beta=p.y/p.x[0]
        beta2=np.dot(p.y,p.y)/(p.x[0]**2)
        gamma=(1-beta2)**(-0.5)
        Ecms=gamma*(plab.x[0]-np.dot(beta,plab.y))
        pcms=-gamma*(beta*plab.x[0]-(plab.y))
        p_new=fv([Ecms,pcms[0],pcms[1],pcms[2]])
        # p0 = (p.x[0]*m_beta.x[0]-(m_beta.y[0]*p.y[0]+m_beta.y[1]*p.y[1]+m_beta.y[2]*p.y[2]))/m_beta.m
        # c1 = (p.x[0]+p0)/(p.m+m_beta.x[0])
        # c2 = p.y-c1*m_beta.y
        # p_new = fv([p0, c2[0], c2[1], c2[2]])
        return p_new
    
    
    def Boost_lab(pcms,p):
        beta=p.y/p.x[0]
        beta2=np.dot(p.y,p.y)/(p.x[0]**2)
        gamma=(1-beta2)**(-0.5)
        Elab=gamma*(pcms.x[0]+np.sqrt(beta2)*pcms.y[2])
        plabz=gamma*(np.sqrt(beta2)*pcms.x[0]+(pcms.y[2]))
        p_new=fv([Elab,pcms.y[0],pcms.y[1],plabz])
        # p0 = (p.x[0]*m_beta.x[0]+(m_beta.y[0]*p.y[0]+m_beta.y[1]*p.y[1]+m_beta.y[2]*p.y[2]))/p.m
        # c1 = (p.x[0]+p0)/(p.m+m_beta.x[0])
        # c2 = p.y+c1*m_beta.y
        # p_new = fv([p0, c2[0], c2[1], c2[2]])
        return p_new
    
    def eps_func(p1,p2,p3,p4):
        X=(p1.x[0]*(p2.x[1]*(p3.x[2]*p4.x[3]-p3.x[3]*p4.x[2])
	           +p2.x[2]*(p3.x[3]*p4.x[1]-p3.x[1]*p4.x[3])
	           +p2.x[3]*(p3.x[1]*p4.x[2]-p3.x[2]*p4.x[1]))
        +p1.x[1]*(p2.x[0]*(p3.x[3]*p4.x[2]-p3.x[2]*p4.x[3])
	           +p2.x[2]*(p3.x[0]*p4.x[3]-p3.x[3]*p4.x[0])
	           +p2.x[3]*(p3.x[2]*p4.x[0]-p3.x[0]*p4.x[2]))
        +p1.x[2]*(p2.x[0]*(p3.x[1]*p4.x[3]-p3.x[3]*p4.x[1])
	           +p2.x[1]*(p3.x[3]*p4.x[0]-p3.x[0]*p4.x[3])
	           +p2.x[3]*(p3.x[0]*p4.x[1]-p3.x[1]*p4.x[0]))
        +p1.x[3]*(p2.x[0]*(p3.x[2]*p4.x[1]-p3.x[1]*p4.x[2])
	           +p2.x[1]*(p3.x[0]*p4.x[2]-p3.x[2]*p4.x[0])
	           +p2.x[2]*(p3.x[1]*p4.x[0]-p3.x[0]*p4.x[1])))
        return X
        

    def anti(self):
        self.m=-self.m
        self.x=self.x
        self.y=self.y
        return self 
        

    
    
    
    
    
    

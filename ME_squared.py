#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 10 19:49:09 2021

@author: asligonulacar
"""
#rename this module
import numpy as np
from FourVector import fv
class ME:
    e = (4*np.pi/137)**0.5 
    def MESqr(M,ch,px,py,p1,p2):
        if ch =='s':     
            prop=(fv.Pro(fv(fv.Add(py,px)),fv(fv.Add(py,px))))**2
        if ch=='t':
            prop=(fv.Pro(fv(fv.Sub(px,p1)),fv(fv.Sub(px,p1))))**2
        return (1/4)*M/prop
    
    def Ref(mx,my,m1,m2,px,py,p1,p2,ch):
        if ch =='s':           
            prop=(fv.Pro(fv(fv.Add(py,px)),fv(fv.Add(py,px))))**2
        if ch=='t':
            prop=(fv.Pro(fv(fv.Sub(px,p1)),fv(fv.Sub(px,p1))))**2
        M_squared = (8*ME.e**4/prop)*(2*mx*my*m1*m2+mx*my*fv.Pro(p1,p2)+m1*m2*fv.Pro(px,py) +
                              fv.Pro(px,p2)*fv.Pro(py, p1)+fv.Pro(px, p1)*fv.Pro(py, p2))
        return M_squared

            
        

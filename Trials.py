#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 11 16:16:17 2021

@author: asligonulacar
"""
from FourVector import fv
import numpy as np
#test momentum generator 
ctmax = 1
ctmin = -1
mx = 0.5
my = 0.5
m1=0
m2=0
E=2
ch='s'          
from Generator import Gen 
p=Gen.P_Gen(E, mx, my, m1, m2, ctmin, ctmax)
print(p[0].x)
print(p[1].x)
print(p[2].x)
print(p[3].x)

#test passed
px=(p[0])
py=(p[1])
p1=(p[2])
p2=(p[3])

#test matrix element
from Amplitude import amp
from ME_squared import ME
M=amp.Z_new(p[0], p[1], p[2], p[3], amp.e, amp.e, amp.e, amp.e, ctmin, ctmax, ch)
HA=ME.MESqr(M,ch,p[0],p[1]) #calculated ME
MF=ME.Ref(mx, my, m1, m2, p[0], p[1], p[2], p[3], ch)[0] #reference ME
print(HA)
print(MF)
#test passed



#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 25 11:46:02 2022

@author: asligonulacar"""

from FourVector import fv
from ME_squared import ME
from Generator import Gen
from Amplitude import amp
import numpy as np
from Channel_Basics import CB 
E=1
me=0
mmu=0.01
ctmin=-1
ctmax=1-1e-7
ch='t'
cR=cR_p=amp.e
cL=cL_p=amp.e
p=Gen.Iso_Gen(E, me, mmu, me, mmu, ctmin, ctmax)

p1=p[0] #electron
p2=p[1] #anti-muon
p3=p[2] #electron
p4=p[3] #anti-muon

#momenta for ME calculation
pm1=p1
pm2=fv([p[1].x[0],-p[1].y[0],-p[1].y[1],-p[1].y[2]]) #-p2 
pm3=fv([p[2].x[0],-p[2].y[0],-p[2].y[1],-p[2].y[2]]) #-p3
pm4=p4

#now I can calculate the matrix element using the pmx momenta, taking advantage 
#of the crossing symmetry. I can now have the reaction ee mumu.
s=fv.Pro(fv(fv.Add(pm1,pm2)),fv(fv.Add(pm1,pm2)))
k=p1.x[0]
Ek=p2.x[0]
ct=1-(fv.Pro(p1,p3)/(k**2))
q=fv.Pro(fv(fv.Sub(p1,p3)),fv(fv.Sub(p1,p3)))#-2*fv.Pro(p1,p3)
a=amp.Z_new(pm1, pm3, pm4, pm2, cR, cL, cR_p, cL_p, ctmin, ctmax,ch)
a2=amp.Z_new(p1, p3, p4, p2, cR, cL, cR_p, cL_p, ctmin, ctmax,ch)
M=(1/4)*a/q**2
M2=(1/4)*a2/q**2
ct=((fv.Pro(p3,p2)/p1.x[0])-p2.x[0])/p1.x[0]
#print(ct/ct1) #ct confirmed


Mref=(8*amp.e**4/q**2)*(fv.Pro(pm1,pm4)*fv.Pro(pm3,pm2)+fv.Pro(pm1,pm2)*fv.Pro(pm3,pm4)+mmu**2*fv.Pro(pm1,pm3))
Mref2=(8*amp.e**4/q**2)*(fv.Pro(pm3,p4)*fv.Pro(p1,pm2)+fv.Pro(pm3,pm2)*fv.Pro(pm1,pm4)+mmu**2*fv.Pro(pm1,pm3))
Mref3=((2*amp.e**4)/(k**2*(1-ct)**2))*((Ek+k)**2+(Ek+k*ct)**2-mmu**2*(1-ct))
print(Mref2/M)

    


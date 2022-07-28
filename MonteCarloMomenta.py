#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 14 19:05:12 2021

@author: asligonulacar
"""
from Generator import Gen
import numpy as np
from FourVector import fv
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d 
#test momentum generator 
ctmax = 1
ctmin = -1
mx = 0
my = 0
m1=0
m2=0
#means p1m is 1
E=2
ch='s'          


N=10000

#plots random points of one of the output momenta
x=[]
y=[]
z=[]
for i in range(N):
    p=Gen.Iso_Gen(E, mx, my, m1, m2, ctmin, ctmax)
    p1=p[2].y
    x.append(p1[0])
    y.append(p1[1])
    z.append(p1[2])

fig = plt.figure(figsize=(3,3))
ax=plt.axes(projection='3d')
# ax=fig.gca(projection='3d')
# u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
# h = np.cos(u)*np.sin(v)
# e = np.sin(u)*np.sin(v)
# f = np.cos(v)
# ax.plot_wireframe(h, e, f)
# ax.tick_params(axis='both', which='major', labelsize=8)
# ax.tick_params(axis='both', which='minor', labelsize=7.5)
ax.set_title('Generated Momenta')
#ax.set_title('Unit Sphere')

ax.scatter3D(x,y,z,c=z,cmap=plt.cm.cool)

# plt.plot(y,z)
# plt.show()

#monte carlo integration of the volume of the sphere that the function covers.
a=[]
b=[]
c=[]
V=0
r=0
for i in range(N):
    p=Gen.Iso_Gen(E, mx, my, m1, m2, ctmin, ctmax)
    p1=p[2].y
    a.append(abs(p1[0]))
    b.append(abs(p1[1]))
    c.append(abs(p1[2]))
    r+=np.sqrt(a[i]**2+b[i]**2+c[i]**2)
r/=N    
V=4*np.pi*r**3/3

R=1
sph=4*np.pi*R**3/3
print(sph/V)
    


    

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  8 19:19:24 2021

@author: asligonulacar
"""
#histogram test using the function x**2
#this is practice to calculate dsigma/dcostheta 
import numpy as np
import matplotlib.pyplot as plt

class histogram:
                
    def diff(x,y):
        diff=[1]
        for i in range(1,len(y)) :
            diff.append((y[i]-y[i-1])/(x[i]-x[i-1]))
        return diff

# def main():
#     x=[]
#     y=[]
#     err=0
#     N=100
#     while err==0 or err>1e-1:
#         N+=100
#         for i in range(N):
#             x.append(np.random.uniform(0,1))
#             y.append(x[i]**2)
#             err=np.var(y)
#             y.sort()
#             x.sort()
#             d=histogram.diff(x,y)
#     plt.hist2d(x,y,bins=[25,25],density=True)
#     plt.colorbar()
#     plt.show()
#     plt.hist2d(x,d,bins=[25,25],density=True)
#     plt.colorbar()
#     return N,err
            
# print(main())

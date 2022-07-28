#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  8 19:19:24 2021

@author: asligonulacar
"""
#histogram test using the function x**2
import numpy as np
import matplotlib.pyplot as plt
def main():
    x=[]
    y=[]
    err=0
    N=100
    while err==0 or err>1e-1:
        N+=100
        for i in range(N):
            x.append(np.random.uniform(0,1))
            y.append(x[i]**2)
            err=np.var(y)
    plt.hist2d(x,y)
    plt.colorbar()
    return N,err
            
print(main())
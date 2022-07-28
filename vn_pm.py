#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  6 23:49:55 2022

@author: asligonulacar
"""

from FourVector import fv
from Generator import Gen
from Amplitude import amp
import numpy as np
import matplotlib.pyplot as plt
import cmath
import math
from numpy import ndarray
from scipy.optimize import curve_fit
from numpy import arange
from scipy.signal import savgol_filter
import datetime
begin_time = datetime.datetime.now()
E=1
mx=0
my=0.000009
m1=0
m2=0.00009
ctmin=-1
ctmax=1

# def smooth_data_np_cumsum_my_average(arr, span):
#     cumsum_vec = np.cumsum(arr)
#     moving_average = (cumsum_vec[2 * span:] - cumsum_vec[:-2 * span]) / (2 * span)
#     front, back = [np.average(arr[:span])], []
#     for i in range(1, span):
#         front.append(np.average(arr[:i + span]))
#         back.insert(0, np.average(arr[-i - span:]))
#     back.insert(0, np.average(arr[-2 * span:]))
#     return np.concatenate((front, moving_average, back))

def objective(x, a, b, c, d, e, f):
 	return a * x + b * x**2 + c*x**3+d*x**5+e*x**6*f

def MC_int_Weighted(E,mx,my,m1,m2,ctmin,ctmax):   
        err=0
        N=0
        M=[]
        Q=[]
        s=[]
        C1=[]
        while err==0 or err>1:
            N=N+1000
            for i in range(N):
                if ((i/N)*100)%10==0:
                    print((str((i/N)*100)) +'%')
                p=Gen.Lab(E, mx, my, m1, m2, ctmin, ctmax)
                ctc=p[4]
                if -1<ctc<0.99:
                    # q=fv(fv.Sub(p[0],p[2]))
                    # Q2=fv.Pro(q,q)
                    # Q.append(abs(Q2))
                    dq2=p[5]
                    Me=(amp.Quasi_neutr_amp(p[0], p[1], p[2], p[3]))[0]
                    m=(amp.Quasi_neutr_amp(p[0], p[1], p[2], p[3])[1])
                    #C=(amp.Quasi_neutr_amp(p[0], p[1], p[2], p[3]))[2]
                    M.append((abs(Me*dq2))/(8*np.pi*m2**2*E**2))
                    s.append(abs(m*dq2)/(8*np.pi*m2**2*E**2))
                    
                else:
                    i=0  
            err=np.var(M)
            plt.xlabel('E_ν [GeV]')
            plt.ylabel('σ [10^-38 cm^2]')
            #plt.xlim(0,1)
            #plt.ylim(0,2)
            #plt.scatter(Q,s,color='blue')
            #print(min(s))
            #plt.scatter(Q,M, color='red')
            
            #plt.errorbar(Q,s,xerr=None,yerr=err,ls='none',color='yellow',alpha=0.5)
            Mtot=(sum(M)/N)
            stot=sum(s)/N
        return stot,Mtot,(stot/Mtot), err

print(MC_int_Weighted(E,mx,my,m1,m2,ctmin,ctmax))
# def Energy(mx,my,m1,m2,ctmin,ctmax):
#     En=[]
#     sigma=[]
#     m=[]
#     #Q=[]
#     N=1
#     for i in range(N):
#         if ((i/N)*100)%10==0:
#             print((str((i/N)*100)) +'%')
#         E=np.random.uniform(0.00001,6) 
#         En.append(E)
#         #sigma.append(MC_int_Weighted(E,mx,my,m1,m2,ctmin,ctmax)[1])
#         m.append(MC_int_Weighted(E,mx,my,m1,m2,ctmin,ctmax)[0])
#     plt.ylim(0,1.4)
#     plt.xlim(0,6)
#     err=np.var(m)
#     #err=np.var(sigma)*50
#     # plt.ylim(0,1*10**-38)
#     plt.grid(False)
#     #print(min(sigma))
#     plt.scatter(En, m)
#     #plt.scatter(En, sigma, color='magenta')
#     # plt.plot(En,m,'s', marker='^', color='magenta')
#     plt.errorbar(En,m, yerr=err, ls='none')
    
#     x,y=En,m

    
#     #yhat = smooth_data_np_cumsum_my_average(y, 10) # window size 51, polynomial order 3
#     # popt, _ = curve_fit(objective, x, y)
#     # a, b, c, d, e, f = popt
#     # x_line = arange(min(x), max(x), 1)
#     # y_line = objective(x_line, a, b, c, d, e, f)
    
#     # plt.plot(x_line, y_line, color='black')
#     # plt.show()
#     # plt.scatter(En,sigma)
#     # plt.scatter(En,m)
#     return None 
# print(Energy(mx,my,m1,m2,ctmin,ctmax))
# print(datetime.datetime.now() - begin_time)

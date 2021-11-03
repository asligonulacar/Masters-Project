#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 14 19:10:44 2021

@author: asligonulacar
"""

import numpy as np

class FourVector:
    g = np.array([[1, 0, 0, 0],  # metric tensor
                  [0, -1, 0, 0],
                  [0, 0, -1, 0],
                  [0, 0, 0, -1]])

    c = 1  # natural units

    def __init__(self, X=None, m=0, gam=0):
        self.x = np.array([X[0], X[1], X[2], X[3]])  # 4 vector
        self.m = m
        self.y = np.array([X[1], X[2], X[3]])

    def Pro(self, other):  # as the name says...
        return self@FourVector.g@other

    def Addition(self, other):  # vector addition
        return np.array([self.x[0]+other.x[0], self.x[1]+other.x[1], self.x[2]+other.x[2], self.x[3]+other.x[3]])

    def Subtraction(self, other):  # vector subtraction because why not
        return np.array([self.x[0]-other.x[0], self.x[1]-other.x[1], self.x[2]-other.x[2], self.x[3]-other.x[3]])

    def EtaForZ(self, other):
        return (2*FourVector.Pro(self, other))**0.5


def SqLam(E, m1, m2):  # SqLam function from Channel Basics
    arg1 = ((E-m1-m2)**2-4*m1*m2)
    if arg1 >= 0:
        return arg1


def Boost(p, m_beta):
    p0 = (p.x[0]*m_beta.x[0]-np.dot(m_beta.y, p.y))/p.m
    c1 = (p.x[0]+p0)/(p.m+m_beta.x[0])
    c2 = p.y-c1*m_beta.y
    p = FourVector([p0, c2[0], c2[1], c2[2]], p0)
    return p

# -------------------------------ATTEMPTING THE LEVI-CIVITA TENSOR------------------------------------------------------
# levi-civita symbol in 4-dimensions-->need to make a 4D tensor.
# The tensor would have to be updated if the user wants to change the metric tensor defined in the FourVector class.
# This tensor definition has been taken from the file 'Tensor Class from github.py'.

Epsilon_Four = ([[[[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]],
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

# -------------------------------THE S FUNCTION------------------------------------------------------------------

def S_plus(p1, p2, k0, k1):
    X = (2*(FourVector.Pro(p1.x, k0.x)*FourVector.Pro(p2.x, k1.x)-FourVector.Pro(p1.x, k1.x)*FourVector.Pro(p2.x, k0.x)+complex(0, 1) *
            np.einsum("ijkl,i,j,k,l->", Epsilon_Four, k0.x, k1.x, p1.x, p2.x)))/(((2*FourVector.Pro(p1.x, k0.x))**0.5)*((2*FourVector.Pro(p2.x, k0.x))**0.5))
    return X


def S_minus(p1, p2, k0, k1):
    X = np.conj((2*(FourVector.Pro(p2.x, k0.x)*FourVector.Pro(p1.x, k1.x)-FourVector.Pro(p2.x, k1.x)*FourVector.Pro(p1.x, k0.x)+complex(0, 1) *
                    np.einsum("ijkl,i,j,k,l->", Epsilon_Four, k0.x, k1.x, p2.x, p1.x)))/(((2*FourVector.Pro(p2.x, k0.x))**0.5)*((2*FourVector.Pro(p1.x, k0.x))**0.5)))
    return X


# -----------------------------------GENERATE MOMENTA-------------------------------------------------------------------------
ctmax = 1
ctmin = -1

mx = 1
my = 0.6
# combined energy of the incoming particles
Ex=3
Ey=1

px=FourVector([Ex,0,0,(Ex**2-mx**2)**0.5],mx)
py=FourVector([Ey,0,0,(Ey**2-my**2)**0.5],my)
p = FourVector(FourVector.Addition(px, py))
m = (p.x[0]**2-p.x[1]**2-p.x[2]**2-p.x[3]**2)**0.5  # calculates mass
p = FourVector(p.x, m)

def P_Generate(px, py, ctmin, ctmax):
    # random numbers ran1 and ran2
    ran1 = np.random.uniform(0, 1)
    ran2 = np.random.uniform(0, 1)
    phi = 2*np.pi*ran2
    costheta = ctmin+(ctmax-ctmin)*ran1
    sintheta = np.sqrt(1-costheta*costheta) 
    new_p = Boost(p, p)  # p in the CMS frame
    new_E = new_p.x[0]
    # momenta of outgoing particles
    m1 = mx 
    m2 = my 
    p1m = ((SqLam(new_E**2, m1**2, m2**2))**0.5/(2*new_E))  
    SphPolars = np.array([sintheta*np.cos(phi), sintheta*np.sin(phi), costheta])
    p1_3 = p1m*SphPolars  # three momentum
    p2_3 = -p1_3
    E1 = (p1m**2+m1**2)**0.5
    E2 = (p1m**2+m2**2)**0.5
    p1 = FourVector([E1, p1_3[0], p1_3[1], p1_3[2]], m1)
    p2 = FourVector([E2, p2_3[0], p2_3[1], p2_3[2]], m2)

    return p1, p2


#show if momenta fills a unit sphere
# dont randomize incoming (tick)
# ------------------------------------NEW HELICITY FUNCTION---------------------------------------------------------------------------------------------------------------------------------------
A = [-1, -1, -1, -1]
initial_A = [-1, -1, -1, -1]


def Helicity(A, initial_A):
    initial_length = str(len(A))
    H = str(2**len(A))  # amount of helicities to be summed over
    #print('need to sum over '+H+' helicities')
    num1 = '0b1'  # preparing to add bin(1) to get the sequence

    # switch -1's with zeroes to prepare to convert to binary
    for index, value in enumerate(A):
        if value == -1:
            A[index] = 0

    A = bin(int(''.join(map(str, (A))), 2) << 1)  # convert A to binary
    num = int(A, 2)  # get the decimal value of A
    # pad it with appropriate amount of zeroes depending on the initial length of A
    A = format(num, '0'+initial_length+'b')

    Helicities = [initial_A]  # initialising the list of heliciies
    # so I can also return the binary of all the helicities cause why not
    Binary = [A]
    for i in range(int(H)-1):
        A = bin(int(A, 2)+int(num1, 2))  # perform the bit addition
        num = int(A, 2)  # calculate the new decimal value of A
        A = format(num, '0'+initial_length+'b')  # re-format it
        Binary.append(A)
        Y = [int(i)for i in A]  # prepare to convert A back to a list
        # convert A back to list form, change 0's back to -1's
        for index, value in enumerate(Y):
            if value == 0:
                Y[index] = -1
        Helicities.append(Y)  # adding all the new helicities to the list

    return Binary

# add return Helicities if you want a list output.


# ----------------------------------------Z FUNCTION-----------------------------------------
k0 = FourVector([1, 0, 0, 1])  # gauge vectors
k1 = FourVector([0, 1, 0, 0])


def Z_new(px, py, p1, p2, cR, cL, cR_p, cL_p, ctmin, ctmax):
    B = Helicity(A, initial_A)
    E1 = FourVector.EtaForZ(px.x, k0.x)  # eta 1
    E2 = FourVector.EtaForZ(py.x, k0.x)  # eta 2
    E3 = FourVector.EtaForZ(p1.x, k0.x)  # eta 3
    E4 = FourVector.EtaForZ(p2.x, k0.x)  # eta 4
    # Mu is negative if its being calculated for an anti-particle.
    Mu1 = px.m/E1
    Mu2 = -py.m/E2
    Mu3 = p1.m/E3
    Mu4 = -p2.m/E4
    S1 = S_plus(p1, px, k0, k1)
    S13 = S_minus(p1, px, k0, k1)
    S2 = S_minus(p2, py, k0, k1)
    S14 = S_plus(p2, py, k0, k1)
    S3 = S_plus(p2, px, k0, k1)
    S15 = S_minus(p2, px, k0, k1)
    S5 = S_minus(py, p1, k0, k1)
    S17 = S_plus(py, p1, k0, k1)
    S6 = S_minus(py, p2, k0, k1)
    S11 = S_plus(py, p2, k0, k1)
    S7 = S_plus(px, p2, k0, k1)
    S10 = S_plus(p1, py, k0, k1)
    S21 = S_minus(p1, py, k0, k1)
    S19 = S_minus(px, p2, k0, k1)
    S23 = S_plus(px, p2, k0, k1)

    for i in B:
        B[-1] = -2*(S1*S2*cR_p*cR-Mu1*Mu2*E3*E4*cR_p*cL-E1*E2*Mu3*Mu4*cL_p*cR)
        B[-2] = -2*E2*cR*(S3*Mu3*cL_p-S1*Mu4*cR_p)
        B[-3] = -2*E1*cR*(S5*Mu4*cL_p-S6*Mu3*cR_p)
        B[-4] = -2*(S7*S5*cL_p*cR-Mu1*Mu2*E3*E4*cL_p*cL-E1*E2*Mu3*Mu4*cR_p*cR)
        B[-5] = -2*Mu4*cR_p*(S1*Mu2*cR-S10*Mu1*cL)
        B[-6] = 0
        B[-7] = -2*(Mu1*Mu4*E2*E3*cL_p*cL+Mu2*Mu3*E1*E4*cR_p*cR -
                    Mu2*Mu4*E1*E3*cL_p*cR-Mu1*Mu3*E2*E4*cR_p*cL)
        B[-8] = -2*E3*cL_p*(S11*Mu1*cL-S7*Mu2*cR)
        B[-9] = -2*E3*cR_p*(S11*Mu1*cR-S23*Mu2*cL)
        B[-10] = -2*(Mu1*Mu4*E2*E3*cR_p*cR+Mu2*Mu3*E1*E4*cL_p *
                     cL-Mu2*Mu4*E1*E3*cR_p*cL-Mu1*Mu3*E2*E4*cL_p*cR)
        B[-11] = 0
        B[-12] = -2*Mu4*cL_p*(S13*Mu2*cL-S21*Mu1*cR)
        B[-13] = -2*(S19*S17*cR_p*cL-Mu1*Mu2*E3*E4 *
                     cR_p*cR-E1*E2*Mu3*Mu4*cL_p*cL)
        B[-14] = -2*E1*cL*(S17*Mu4*cR_p-S11*Mu3*cL_p)
        B[-15] = -2*E2*cL*(S15*Mu3*cR_p-S13*Mu4*cL_p)
        B[-16] = -2*(S13*S14*cL_p*cL-Mu1*Mu2*E3*E4 *
                     cL_p*cR-E1*E2*Mu3*Mu4*cR_p*cL)
    C = []
    for i in range(0, 16):
        C.append(B[i]*np.conj(B[i]))
    C = sum(C)
    # t=FourVector.Subtraction(px,py)
    # t_4=FourVector.Pro(t,t)*FourVector.Pro(t,t)
    q = FourVector.Addition(py, px)
    q_4 = FourVector.Pro(q, q)*FourVector.Pro(q, q)
    e = 1  # (4*np.pi/137)**0.5 #from the fine structure constant
    M = 0.25*e**4*C/q_4
    m1 = p1.m
    m2 = p2.m
    M_squared = (8*e**4/q_4)*(2*mx*my*m1*m2+mx*my*FourVector.Pro(p1.x, p2.x)+m1*m2*FourVector.Pro(px.x, py.x) +
                              FourVector.Pro(px.x, p2.x)*FourVector.Pro(py.x, p1.x)+FourVector.Pro(px.x, p1.x)*FourVector.Pro(py.x, p2.x))
    # I've tried to implement isotropic weight but I dont know if its right
    Iso_Weight = 4*np.pi*(SqLam(m**2, mx**2, my**2))/(ctmax-ctmin)
    An_Iso_Weight = 0  # haven't implemented this yet.
    return C, M, Iso_Weight, M_squared

# Z takes those inputs and runs a loop over all helicities.
# construct a skeleton first and put the Z functions together by hand.

# --------------------CALCULATING THE CROSS SECTION-MONTE CARLO------------------------------------------------------------------------------------------------------------------------------------------------
def main():
    N = 1
    c = []
    c_check = []
    for i in range(N):
        p = P_Generate(px, py, ctmin, ctmax)
        z = Z_new(px, py, p[0], p[1], 1, 1, 1, 1, ctmin, ctmax)
        z1 = z[1]  # this gives me the matrix element squared from helicity amplitudes
        z2 = z[2]  # this gives me the isotropic weight
        z3 = z[3]  # this gives me the reference matrix element
        # this is the array I want to sum over to get integral result
        c.append(z1*z2)
        c_check.append(z3*z2)

    c = sum(c)/N  # divide by number of samples to get cross section
    c_check = sum(c_check)/N
    return c, c_check

# I need to produce error.
#virtual classes and give specific implementations names.

#print(main())





































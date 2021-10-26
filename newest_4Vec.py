#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 14 19:10:44 2021

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
        return (2*FourVector.Pro(self,other))**0.5
        

#---------------------------DEFINING GAMMA MATRICES----------------------------------------------------------------    
#Pauli Matrices
Pauli1=np.array([[0,1],[1,0]])
Pauli2=np.array([[0,-complex(0,1)],[complex(0,1),0]])
Pauli3=np.array([[1,0],[0,-1]])

#Gamma Matrices constructed from Pauli Matrices using np.block()
Gamma0=np.block([[np.identity(2),np.zeros((2,2))],[np.zeros((2,2)),-np.identity(2)]])
Gamma1=np.block([[np.zeros((2,2)),Pauli1],[-Pauli1,np.zeros((2,2))]])
Gamma2=np.block([[np.zeros((2,2)),Pauli2],[-Pauli2,np.zeros((2,2))]])
Gamma3=np.block([[np.zeros((2,2)),Pauli3],[-Pauli3,np.zeros((2,2))]])
Gamma5=complex(0,1)*(Gamma0@Gamma1@Gamma2@Gamma3)


#-------------------------------ATTEMPTING THE LEVI-CIVITA TENSOR------------------------------------------------------
#levi-civita symbol in 4-dimensions-->need to make a 4D tensor.
#The tensor would have to be updated if the user wants to change the metric tensor defined in the FourVector class.
#This tensor definition has been taken from the file 'Tensor Class from github.py'. 

Epsilon_Four=([[[[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]], 
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
    
#-------------------------------ATTEMPTING THE S FUNCTION------------------------------------------------------------------
#testing out the S function

def S_plus(p1,p2,k0,k1): 
    X=(2*(FourVector.Pro(p1.x,k0.x)*FourVector.Pro(p2.x,k1.x)-FourVector.Pro(p1.x,k1.x)*FourVector.Pro(p2.x,k0.x)+complex(0,1)*np.einsum("ijkl,i,j,k,l->",Epsilon_Four,k0.x,k1.x,p1.x,p2.x)))/(((2*FourVector.Pro(p1.x,k0.x))**0.5)*((2*FourVector.Pro(p2.x,k0.x))**0.5))
    return X

def S_minus(p1,p2,k0,k1):
    X=np.conj((2*(FourVector.Pro(p2.x,k0.x)*FourVector.Pro(p1.x,k1.x)-FourVector.Pro(p2.x,k1.x)*FourVector.Pro(p1.x,k0.x)+complex(0,1)*np.einsum("ijkl,i,j,k,l->",Epsilon_Four,k0.x,k1.x,p2.x,p1.x)))/(((2*FourVector.Pro(p2.x,k0.x))**0.5)*((2*FourVector.Pro(p1.x,k0.x))**0.5)))
    return X


k0=FourVector([1,0,0,1]) #gauge vectors
k1=FourVector([0,1,0,0])

#-----------------------------------GENERATE MOMENTA-------------------------------------------------------------------------
#NOTE: Its important to check that we have px=-py with px,py being three momenta.
#UPDATE: The IncomingMomentumX makes sure that this is satisfied now.

def SqLam(E,m1,m2): #SqLam function from Channel Basics
  arg1=((E-m1-m2)**2-4*m1*m2)
  if arg1>=0:
      return arg1 

def Boost(p,m_beta):
    p0=(p.x[0]*m_beta.x[0]-np.dot(m_beta.y,p.y))/p.m
    c1=(p.x[0]+p0)/(p.m+m_beta.x[0])
    c2=p.y-c1*m_beta.y
    p=FourVector([p0,c2[0],c2[1],c2[2]],p0)
    return p

#masses of incoming  particles
mx=0.003
my=0.004
mxy=mx+my

#combined energy of the incoming particles
E=0.5

#generating random incoming momenta using random numbers for the angles
def IncomingMomentumX(E,mx,my,SqLam):
       p1x=(SqLam(E**2,mx**2,my**2))**0.5/(2*E)
       ran1=np.random.uniform(0,1)
       ran2=np.random.uniform(0,1)
       costheta=1-2*ran1
       sintheta=2*(ran1*(1-ran1))**0.5
       phi=2*np.pi*ran2
       SphPolars=np.array([sintheta*np.cos(phi),sintheta*np.sin(phi),costheta])
       p1x_3=p1x*SphPolars
       Ex=(np.dot(p1x_3,p1x_3)+mx**2)**0.5
       px=FourVector([Ex,p1x_3[0],p1x_3[1],p1x_3[2]],mx)
       return px

#px and py represent the incoming momenta       
px=IncomingMomentumX(E,mx,my,SqLam)
py_3=-px.y
Ey=(np.dot(py_3,py_3)+my**2)**0.5
py=FourVector([Ey,py_3[0],py_3[1],py_3[2]],my)
p=FourVector(FourVector.Addition(px,py))
m_beta=p #this keeps the original value of p
m=(p.x[0]**2-p.x[1]**2-p.x[2]**2-p.x[3]**2)**0.5 #calculates mass 
p=FourVector(p.x,m)


new_p=Boost(p,m_beta) #p in the CMS frame 

#masses of outgoing particles
m1=(mx/mxy)*p.m
m2=(my/mxy)*p.m
new_E=new_p.x[0]
#print(new_p.x)
def NewP1(new_E,m1,m2,SqLam):
    #CMS frame
    p1m=(SqLam(new_E**2,m1**2,m2**2))**0.5/(2*new_E) 
    theta=np.random.uniform(0,2)*np.pi
    phi=np.random.uniform(0,1)*np.pi 
    #rho1 and rho2 will be important when it comes to the Monte Carlo integration
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
En1=(np.dot(p1.y,p1.y)+m1**2)**0.5
#Calcaulating p2 in the CMS frame.
p2_3=-p1.y
En2=(np.dot(p2_3,p2_3)+m2**2)**0.5
p2=FourVector([En2,p2_3[0],p2_3[1],p2_3[2]])

# print(px.x+py.x)

# print(p1.x)

#print(((p1.x[0]**2-p1.x[1]**2-p1.x[2]**2-p1.x[3]**2)**0.5))
#alpha is the fine structure constant 
#e is the strength of the coupling 
#costheta max is deflection angle. costheta max closer to one  the more the cross section explodes
#costheta min =-1 the electron goes back to where it comes from.


#------------------------------------NEW HELICITY FUNCTION---------------------------------------------------------------------------------------------------------------------------------------
A=[-1,-1,-1,-1] 
initial_A=[-1,-1,-1,-1] 

def Helicity(A,initial_A):
    initial_length=str(len(A))
    H=str(2**len(A)) #amount of helicities to be summed over
    #print('need to sum over '+H+' helicities')
    num1='0b1' #preparing to add bin(1) to get the sequence
    
    for index,value in enumerate(A): #switch -1's with zeroes to prepare to convert to binary
        if value==-1:
           A[index]=0
           
    A=bin(int(''.join(map(str, (A))), 2) << 1) #convert A to binary 
    num=int(A,2) #get the decimal value of A
    A=format(num,'0'+initial_length+'b')  #pad it with appropriate amount of zeroes depending on the initial length of A
    
    Helicities=[initial_A] #initialising the list of heliciies 
    Binary=[A] #so I can also return the binary of all the helicities cause why not
    for i in range(int(H)-1):
        A=bin(int(A,2)+int(num1,2)) #perform the bit addition 
        num=int(A,2) #calculate the new decimal value of A
        A=format(num,'0'+initial_length+'b') #re-format it 
        Binary.append(A)
        Y=[int(i)for i in A] #prepare to convert A back to a list
        for index,value in enumerate(Y): #convert A back to list form, change 0's back to -1's
            if value==0:
                Y[index]=-1
        Helicities.append(Y) #adding all the new helicities to the list 
                
    return Binary
#add return Helicities if you want a list output.

#----------------------------------------Z FUNCTION-----------------------------------------
E1=FourVector.EtaForZ(px.x,k0.x) #eta 1
E2=FourVector.EtaForZ(py.x,k0.x) #eta 2
E3=FourVector.EtaForZ(p1.x,k0.x) #eta 3
E4=FourVector.EtaForZ(p2.x,k0.x) #eta 4

#Mu is negative if its being calculated for an anti-particle.
Mu1=mx/E1   
Mu2=-my/E2
Mu3=m1/E3
Mu4=-m2/E4


#setting the S function values beforehand to make Z function calculations more efficient.
S1=S_plus(p1,px,k0,k1)
S13=S_minus(p1,px,k0,k1)
S2=S_minus(p2,py,k0,k1)
S14=S_plus(p2,py,k0,k1)
S3=S_plus(p2,px,k0,k1)
S15=S_minus(p2,px,k0,k1)
S5=S_minus(py,p1,k0,k1)
S17=S_plus(py,p1,k0,k1)
S6=S_minus(py,p2,k0,k1)
S11=S_plus(py,p2,k0,k1)
S7=S_plus(px,p2,k0,k1)
S10=S_plus(p1,py,k0,k1)
S21=S_minus(p1,py,k0,k1)
S19=S_minus(px,p2,k0,k1)
S23=S_plus(px,p2,k0,k1)

def Z_new(Helicity,px,py,p1,p2,cR,cL,cR_p,cL_p):
    B=Helicity(A,initial_A)
    for i in B:
        B[-1]=-2*(S1*S2*cR_p*cR-Mu1*Mu2*E3*E4*cR_p*cL-E1*E2*Mu3*Mu4*cL_p*cR)
        B[-2]=-2*E2*cR*(S3*Mu3*cL_p-S1*Mu4*cR_p)
        B[-3]=-2*E1*cR*(S5*Mu4*cL_p-S6*Mu3*cR_p)
        B[-4]=-2*(S7*S5*cL_p*cR-Mu1*Mu2*E3*E4*cL_p*cL-E1*E2*Mu3*Mu4*cR_p*cR)
        B[-5]=-2*Mu4*cR_p*(S1*Mu2*cR-S10*Mu1*cL)
        B[-6]=0
        B[-7]=-2*(Mu1*Mu4*E2*E3*cL_p*cL+Mu2*Mu3*E1*E4*cR_p*cR-Mu2*Mu4*E1*E3*cL_p*cR-Mu1*Mu3*E2*E4*cR_p*cL)
        B[-8]=-2*E3*cL_p*(S11*Mu1*cL-S7*Mu2*cR)
        B[-9]=-2*E3*cR_p*(S11*Mu1*cR-S23*Mu2*cL)
        B[-10]=-2*(Mu1*Mu4*E2*E3*cR_p*cR+Mu2*Mu3*E1*E4*cL_p*cL-Mu2*Mu4*E1*E3*cR_p*cL-Mu1*Mu3*E2*E4*cL_p*cR)
        B[-11]=0
        B[-12]=-2*Mu4*cL_p*(S13*Mu2*cL-S21*Mu1*cR)
        B[-13]=-2*(S19*S17*cR_p*cL-Mu1*Mu2*E3*E4*cR_p*cR-E1*E2*Mu3*Mu4*cL_p*cL)
        B[-14]=-2*E1*cL*(S17*Mu4*cR_p-S11*Mu3*cL_p)
        B[-15]=-2*E2*cL*(S15*Mu3*cR_p-S13*Mu4*cL_p)
        B[-16]=-2*(S13*S14*cL_p*cL-Mu1*Mu2*E3*E4*cL_p*cR-E1*E2*Mu3*Mu4*cR_p*cL)
    C=[]
    for i in range(0,16):
        C.append(B[i]*np.conj(B[i]))
    return sum(C)

t=FourVector.Subtraction(px,py)
t_4=FourVector.Pro(t,t)*FourVector.Pro(t,t)
q=FourVector.Addition(py,px)
q_4=FourVector.Pro(q,q)*FourVector.Pro(q,q)


e=1#(4*np.pi/137)**0.5 #from the fine structure constant
def M_sqr(Z_new):
    M=0.25*e**4*Z_new(Helicity,px,py,p1,p2,1,1,1,1)/t_4
    return M

print(M_sqr(Z_new))
M_squared=(8*e**4/t_4)*(2*mx*my*m1*m2+mx*my*FourVector.Pro(p1.x,p2.x)+m1*m2*FourVector.Pro(px.x, py.x)+FourVector.Pro(px.x,p2.x)*FourVector.Pro(py.x,p1.x)+FourVector.Pro(px.x,p1.x)*FourVector.Pro(py.x,p2.x))
print(abs(M_squared))
#Z takes those inputs and runs a loop over all helicities.
#construct a skeleton first and put the Z functions together by hand.

#---------------------FROM MATRIX ELEMENT TO CROSS SECTION-----------------------------------------------------------------------------------------------------------------------------------------------------
#this is d(rho), I dont even need this rip.
def Diff_Cross(M_sqr,dcostheta,dphi):
    s=4*En1*En2
    d_omega=dcostheta*dphi
    pi=(px.y[0]**2+px.y[1]**2+px.y[2]**2)**0.5
    pf=(p1.y[0]**2+p1.y[1]**2+p1.y[2]**2)**0.5
    return(1/(8*np.pi)**2)*(1/s)*(pf/pi)*M_sqr*d_omega
    

#--------------------CALCULATING THE CROSS SECTION-MONTE CARLO------------------------------------------------------------------------------------------------------------------------------------------------
#This should be the monte carlo integration bit.

#number of samples 
N=1000
#function to Monte Carlo
J=4*np.pi

#the general idea is here however have to modify line 327 to fit my matrix element function 
def Monte_Carlo(f,a,b,N):
    b=1
    a=0
    subsets=np.arange(0,N+1,N/100)
    steps=N/100
    u=np.zeros(N)
    for i in range(100):
        start=int(subsets[i])
        end=int(subsets[i+1])
        u[start:end]=np.random.uniform(low=i/100,high=(i+1)/100,size=end-start)
    np.random.shuffle(u)
    u_f=f(a+(b-a)*u)
    s=((b-a)/N)*u_f.sum()
    return s
    
print(Monte_Carlo((M_sqr(Z_new)),0,1,N))





















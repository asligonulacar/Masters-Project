#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import numpy as np

class FourVector:
    
    g=np.array([[1,0,0,0], #metric tensor
               [0,-1,0,0],
               [0,0,-1,0],
               [0,0,0,-1]])
    
    c=1 #natural units
    
    def __init__(self,X=None,m=0,gam=0): 
        
        self.x=np.array([X[0],X[1],X[2],X[3]]) #4 vector
        
        self.mass=m
     
    #def LorentzTransform(self): 
        #general Lorentz transform matrix for 
        
    def FourDotProduct(self,other): #as the name says...
        return self@FourVector.g@other  
        #return np.trace(self*FourVector.g*other)
    
    def Addition(self,other): #vector addition
        return np.array([self.x[0]+other.x[0],self.x[1]+other.x[1],self.x[2]+other.x[2],self.x[3]+other.x[3]])  

    def Subtraction(self,other): #vector subtraction because why not
        return np.array([self.x[0]-other.x[0],self.x[1]-other.x[1],self.x[2]-other.x[2],self.x[3]-other.x[3]])  
    
    def EtaForZ(self,other):
        return (2*FourVector.FourDotProduct(self,other))**0.5

            
            

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
    X=(2*(FourVector.FourDotProduct(p1.x,k0.x)*FourVector.FourDotProduct(p2.x,k1.x)-FourVector.FourDotProduct(p1.x,k1.x)*FourVector.FourDotProduct(p2.x,k0.x)+complex(0,1)*np.einsum("ijkl,i,j,k,l->",Epsilon_Four,k0.x,k1.x,p1.x,p2.x)))/(((2*FourVector.FourDotProduct(p1.x,k0.x))**0.5)*((2*FourVector.FourDotProduct(p2.x,k0.x))**0.5))
    return X

def S_minus(p1,p2,k0,k1):
    X=complex(0,1)*(2*(FourVector.FourDotProduct(p2.x,k0.x)*FourVector.FourDotProduct(p1.x,k1.x)-FourVector.FourDotProduct(p2.x,k1.x)*FourVector.FourDotProduct(p1.x,k0.x)+complex(0,1)*np.einsum("ijkl,i,j,k,l->",Epsilon_Four,k0.x,k1.x,p2.x,p1.x)))/(((2*FourVector.FourDotProduct(p2.x,k0.x))**0.5)*((2*FourVector.FourDotProduct(p1.x,k0.x))**0.5))
    return X


k0=FourVector([1,0,0,1]) #gauge vectors
k1=FourVector([0,1,0,0])

#four momentum conservation 
p1=FourVector([1,0,-0.1,0]) #momentum vectors, have to be on their mass shell. take energy momentum equals same. cosphi and phi the two outgoing momentum have to be back to back 
p2=FourVector([1,0,0.7395917,0])
p3=FourVector([1,0,-0.1,0])
p4=FourVector([1,0,0.7395917,0])

E1=FourVector.EtaForZ(p1.x,k0.x) #eta 1
E2=FourVector.EtaForZ(p2.x,k0.x) #eta 2
E3=FourVector.EtaForZ(p3.x,k0.x) #eta 3
E4=FourVector.EtaForZ(p4.x,k0.x) #eta 4

Mu1=p1.mass/E1   #all mus are zero for now.
Mu2=p2.mass/E2
Mu3=p3.mass/E3
Mu4=p4.mass/E4


#-----------------------------ATTEMPTING THE Z FUNCTIONS-----------------------------------------------------------------------
#helicity summer -- I think I should define a separate function for this
#list-dictionary of spins
#loop over loops 
#bit addition -binary 2**N-1, sum over them in binary space
#simple loop you have one number and decompose it into bits.
#Equations to code: A27, set all couplings to one and mu's to zero.


#trying out a very simple helicity dictionary method for the 8 spins given in the paper. 

Helicities={'Helicity_1':np.array([1,1,1,1]),
            'Helicity_2':np.array([1,1,1,-1]),
            'Helicity_3':np.array([1,1,-1,1]),
            'Helicity_4':np.array([1,1,-1,-1]),
            'Helicity_5':np.array([1,-1,1,1]),
            'Helicity_6':np.array([1,-1,1,-1]),
            'Helicity_7':np.array([1,-1,-1,1]),  
            'Helicity_8':np.array([1,-1,-1,-1]),
            'Helicity_9':np.array([-1,-1,-1,-1]),   
            'Helicity_10':np.array([-1,-1,-1,1]),
            'Helicity_11':np.array([-1,-1,1,-1]),
            'Helicity_12':np.array([-1,-1,1,1]),
            'Helicity_13':np.array([-1,1,-1,-1]),
            'Helicity_14':np.array([-1,1,-1,1]),
            'Helicity_15':np.array([-1,1,1,-1]),
            'Helicity_16':np.array([-1,1,1,1])}


#will figure out a more efficient way to do this.

def Z(Helicities,p1,p2,p3,p4,cR,cL,cR_p,cL_p): #p in cR_p and cL_p means prime. The c's of course give the coupling.
    
    if Helicities=='Helicity_1':
        return -2*(S_plus(p3,p1,k0,k1)*(S_minus(p4,p2,k0,k1))*cR_p*cR-(Mu1*Mu2*E3*E4*cR_p*cL)-(E1*E2*Mu3*Mu4*cL_p*cR))
    elif Helicities=='Helicity_2':
        return -2*E2*cR*(S_plus(p4,p1,k0,k1)*Mu3*cL_p-S_plus(p3,p1,k0,k1)*Mu4*cR_p)
    elif Helicities=='Helicity_3':
        return -2*E1*cR*(S_minus(p2,p3,k0,k1)*Mu4*cL_p-S_minus(p2,p4,k0,k1)*Mu3*cR_p)
    elif Helicities=='Helicity_4':
        return -2*(S_plus(p1,p4,k0,k1)*S_minus(p2,p3,k0,k1)*cL_p*cR-Mu1*Mu2*E3*E4*cL_p*cL-E1*E2*Mu3*Mu4*cR_p*cR)
    elif Helicities=='Helicity_5':
        return -2*Mu4*cR_p*(S_plus(p3,p1,k0,k1)*Mu2*cR-S_plus(p3,p2,k0,k1)*Mu1*cL)
    elif Helicities=='Helicity_6': #this one is actually zero, even without mu=0.
        return 0
    elif Helicities=='Helicity_7':
        return -2*(Mu1*Mu4*E2*E3*cL_p*cL+Mu2*Mu3*E1*E4*cR_p*cR-Mu2*Mu4*E1*E3*cL_p*cR-Mu1*Mu3*E2*E4*cR_p*cL)
    elif Helicities=='Helicity_8':
        return -2*E3*cL_p*(S_plus(p2,p4,k0,k1)*Mu1*cL-S_plus(p1,p4,k0,k1)*Mu2*cR)
    elif Helicities=='Helicity_9':
        return -2*(S_minus(p3,p1,k0,k1)*(S_plus(p4,p2,k0,k1))*cL_p*cL-(Mu1*Mu2*E3*E4*cL_p*cR)-(E1*E2*Mu3*Mu4*cR_p*cL))
    elif Helicities=='Helicity_10':
        return -2*E2*cL*(S_minus(p4,p1,k0,k1)*Mu3*cR_p-S_minus(p3,p1,k0,k1)*Mu4*cL_p)
    elif Helicities=='Helicity_11':
        return -2*E1*cL*(S_plus(p2,p3,k0,k1)*Mu4*cR_p-S_plus(p2,p4,k0,k1)*Mu3*cL_p)
    elif Helicities=='Helicity_12':
        return -2*(S_minus(p1,p4,k0,k1)*S_plus(p2,p3,k0,k1)*cR_p*cL-Mu1*Mu2*E3*E4*cR_p*cR-E1*E2*Mu3*Mu4*cL_p*cL)
    elif Helicities=='Helicity_13':
        return -2*Mu4*cL_p*(S_minus(p3,p1,k0,k1)*Mu2*cL-S_minus(p3,p2,k0,k1)*Mu1*cR)
    elif Helicities=='Helicity_14':
        return 0
    elif Helicities=='Helicity_15':
        return -2*(Mu1*Mu4*E2*E3*cR_p*cR+Mu2*Mu3*E1*E4*cL_p*cL-Mu2*Mu4*E1*E3*cR_p*cL-Mu1*Mu3*E2*E4*cL_p*cR)
    elif Helicities=='Helicity_16':
        return -2*E3*cR_p*(S_plus(p2,p4,k0,k1)*Mu1*cR-S_plus(p1,p4,k0,k1)*Mu2*cL)
    


print((Z('Helicity_12',p1,p2,p3,p4,1,1,1,1)))

Z_Sum=(abs(Z('Helicity_1',p1,p2,p3,p4,1,1,1,1))**2+abs(Z('Helicity_4',p1,p2,p3,p4,1,1,1,1))**2+abs(Z('Helicity_9',p1,p2,p3,p4,1,1,1,1))**2+abs(Z('Helicity_12',p1,p2,p3,p4,1,1,1,1))**2)


print(Z_Sum)

q=FourVector.Addition(p1,p2)

q_4=FourVector.FourDotProduct(q,q)*FourVector.FourDotProduct(q,q)

e=1

Z_squared=(32*e**4/q_4)*(FourVector.FourDotProduct(p1.x,p4.x)*FourVector.FourDotProduct(p2.x,p3.x)+FourVector.FourDotProduct(p1.x,p3.x)*FourVector.FourDotProduct(p2.x,p4.x))

print(abs(Z_squared))


#------------------------------------NEW HELICITY FUNCTION---------------------------------------------------------------------------------------------------------------------------------------
A=[-1,-1,-1,-1] 
initial_A=[-1,-1,-1,-1] 

def Helicities(A,initial_A):
    initial_length=str(len(A))
    H=str(2**len(A)) #amount of helicities to be summed over
    print('need to sum over '+H+' helicities')
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
                
    return Helicities,Binary
    
#-----------------------------------NEW Z FUNCTION--------------------------------------------------------------------------------------













    

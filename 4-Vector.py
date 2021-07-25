import numpy as np

class FourVector:
    
    g=np.array([[1,0,0,0], #metric tensor
               [0,-1,0,0],
               [0,0,-1,0],
               [0,0,0,-1]])
    
    
    m=complex(0.1,0) #trying out invariant mass, will set this as an incident variable later
    
    def __init__(self,t=0,R=None,V=None): 
        
        self.x=np.array([t,R[0],R[1],R[2]]) #space-time 4 vector
        
        if R and V is not None:
            
            self.R=np.array([R[0],R[1],R[2]]) #spatial part of the vector
            
            self.V=np.array([V[0],V[1],V[2]]) #velocity is given in c's, natural units
            
            self.p=np.array(FourVector.m*self.V) #3-momentum 
            
            self.gam=1/np.sqrt(1-np.linalg.norm(V)**2) #gamma factor of the vector
            
            self.normR=np.array([0,2]) if np.linalg.norm(R)==0 else R/np.linalg.norm(R) #calculates the norm of 3-position
            
            self.normV=np.array([0,2]) if np.linalg.norm(V)==0 else V/np.linalg.norm(V) #calculates the norm of 3-velocity
     
    def SwitchToCoM(self,other): #this is currently for only two particles, will update later.
        for i in range(0,2):
            CoMR=np.array([(self.R[i-1]+other.R[i-1])/2,(self.R[i]+other.R[i])/2,(self.R[i+1]+other.R[i+1])/2])            
        return CoMR
    
    def LorentzTransform(self): 
        #general Lorentz transform matrix for velocity
        self.TransfromMatrix=np.array([[self.gam,-self.gam*self.V[0],-self.gam*self.V[1],-self.gam*self.V[2]],
                                 [-self.gam*self.V[0],1+(self.gam-1)*self.normV[0]**2,(self.gam-1)*self.normV[0]*self.normV[1],(self.gam-1)*self.normV[0]*self.normV[2]],
                                 [-self.gam*self.V[1],(self.gam-1)*self.normV[0]*self.normV[1],1+(self.gam-1)*self.normV[1]**2,(self.gam-1)*self.normV[1]*self.normV[2]],
                                 [-self.gam*self.V[2],(self.gam-1)*self.normV[2]*self.normV[0],(self.gam-1)*self.normV[2]*self.normV[1],1+(self.gam-1)*self.normV[2]**2]])
    
        return np.dot(self.TransformMatrix,self.x)
    
    def FourDotProduct(self,other): #as the name says...
        return self.x@FourVector.g@other.x
        #return self.x[0]*other.x[0]-self.x[1]*other.x[1]-self.x[2]*other.x[2]-self.x[3]*other.x[3] 
        #return np.trace(self.x*FourVector.g*other.x), this one also seems to be accurate but not as straightforward as the one above.
    
    def FourVelocity(self): #get 4-velocity from self.x type of input
        return np.array([self.gam*FourVector.c,self.gam*self.V[0],self.gam*self.V[1],self.gam*self.V[2]])
    
    def FourMomentum(self): #get4-momentum from 4-velocity
        return FourVector.m*FourVector.FourVelocity(self)
    
    def Addition(self,other): #vector addition
        return np.array([self.x[0]+other.x[0],self.x[1]+other.x[1],self.x[2]+other.x[2],self.x[3]+other.x[3]])  

    def Subtraction(self,other): #vector subtraction because why not
        return np.array([self.x[0]-other.x[0],self.x[1]-other.x[1],self.x[2]-other.x[2],self.x[3]-other.x[3]])  
    
    def Energy(self): #energy from E**2 = p**2 + m**2, would need to tweak this later.        
        E=((FourVector.m)**2+np.dot(self.p,self.p))**0.5
        return E 
      
        

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
N=FourVector(0,[2,5,1],[0.01,0.02,0.03])

#--------------------------ASSIGNING REPRESENTATIONS TO SPINORS----------------------------------------------------

#Defining Pauli Spinors
A=np.array([0,1])
B=np.array([1,0])

X1=A.reshape(-1,1) #gives chi1, (0,1)
X2=B.reshape(-1,1) #gives chi2, (1,0)


#Introducing Dirac Spinors-for a specific vector N for now
C1=np.array([N.p[2]/(FourVector.Energy(N)+FourVector.m),complex(N.p[0],N.p[1])/(FourVector.Energy(N)+FourVector.m)])
D1=C1.reshape(-1,1)

C2=np.array([complex(N.p[0],-N.p[1])/(FourVector.Energy(N)+FourVector.m),-N.p[2]/(FourVector.Energy(N)+FourVector.m)])
D2=C2.reshape(-1,1)

eta=(FourVector.Energy(N)+FourVector.m)**0.5

u1=eta*np.concatenate([X2,D1]) #u and v dirac spinors for a given mass, momentum and energy 
u2=eta*np.concatenate([X1,D2])
v1=eta*np.concatenate([D1,X2])
v2=eta*np.concatenate([D2,X1])

#-------------------------------ATTEMPTING THE S FUNCTION----------------------------------------------------------
#Would it be accurate to use np.dot() to multiply a k vector and a p vector?












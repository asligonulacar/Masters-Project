#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np

A=[-1,-1,-1,-1] 

def Swap(A): #replacing -1's with a 0 to get ready to convert the list into binary.
    for index,value in enumerate(A):
        if value==-1:
           A[index]=0
    return list(A)


def Convert(Swap,A): #converting into binary numbers.     
    A= bin(int(''.join(map(str, Swap(A))), 2) << 1)  
    return A
    

#print(Convert(Swap,A))


def main(A,Convert): #addition of binary numbers to obtain sequence 
    num1='0b001'
    A='0b000'
    X=bin(int(A,2)+int(num1,2))
    return X


#print(main(A,Convert(Swap, A)))


def Invert(A): #inverts the binary value to list 
    A='0b0001'
    X=[int(i)for i in A[2:]]
    for index,value in enumerate(X):
        if value==0:
            X[index]=-1
    return list(X)

#print(Invert(A))

#Notes: Need to manually update the value of A at the moment which is pretty shit.
#Need to change that so its automated somehow 
#bitwise comparison



#------------------------ADDING ALL THE FUNCTIONS TOGETHER-----------------------------------------------------------------

A=[-1,-1,-1,1] 
def alg(A):
    initial_length=str(len(A))
    num1='0b1' #preparing to add bin(1) to get the sequence
    for index,value in enumerate(A): #switch -1's with zeroes to prepare to convert to binary
        if value==-1:
           A[index]=0
           A=A
    A=bin(int(''.join(map(str, (A))), 2) << 1) #convert A to binary 
    num=int(A,2) #get the decimal value of A
    A=format(num,'0'+initial_length+'b')  #pad it with appropriate amount of zeroes depending on the initial length of A
    
    A=bin(int(A,2)+int(num1,2)) #perform the bit addition 
    num=int(A,2) #calculate the new decimal value of A
    A=format(num,'0'+initial_length+'b') #re-format it
    
    Y=[int(i)for i in A] #prepare to convert A back to a list
    for index,value in enumerate(Y): #convert A back to list form, change 0's back to -1's
        if value==0:
            Y[index]=-1
    return (Y)

#print(alg(A))

#now I need to input how manyt helicities to sum over into the function.






#----------------------HELICITY TRY OUTS-------------------------------------------------------------------------------

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
    



print((Helicities(A,initial_A)))



#----------------------Z FUNCTION TRY OUTS-------------------------------------------------------------------------------




























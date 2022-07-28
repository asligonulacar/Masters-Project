#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 21 18:28:21 2021

@author: asligonulacar
"""
from FourVector import fv
from ME_squared import ME
from Generator import Gen
from Amplitude import amp
import numpy as np
from Channel_Basics import CB
import matplotlib.pyplot as plt
from histogram_test import histogram
import pylab
import cmath
from scipy import stats 
E=1
mx=0
my=0.939
m1=0
m2=0.939#in GeV
ctmin=-1
ctmax=0.99
ch='s'
cR=cL=cR_p=cL_p=amp.e
 
def MC_int_Weighted(E,ch,mx,my,m1,m2,ctmin,ctmax,cR,cL,cR_p,cL_p):   
        err=0
        #sigma=np.array([])
        sigma2=[]
        #cross=np.array([])
        N=0
        ctl=np.array([])
        RB=[]
        Q=[]
        #E3=[]
        #F2=[]
       # F1=[]
        # Gm=[]
        # Ge=[]
        Mref=[]
        # A=[]
        while err==0 or err>10000:#0.05*(sum(sigma2)/len(ctl)):
            N=N+500
            for i in range(N):
                if ((i/N)*100)%10==0:
                    print((str((i/N)*100)) +'%')
                p=Gen.Lab(E, mx, my, m1, m2, ctmin, ctmax)
                ctc=p[4]
                if -1<ctc<0.99:
                    
                    pm1=p[0]
                    pm2=p[1]
                    pm3=p[2]
                    pm4=p[3]
                    q=fv(fv.Sub(pm2,pm4))
                    Q2=abs(fv.Pro(q,q))
                    
                    
                    #Q.append(Q2)
                    t=Q2/(4*my**2)
                    #p2=fv([p[1].x[0]/2,p[1].x[1]/2,p[1].x[2]/2,p[1].x[3]/2])
                    #p4=fv([p[3].x[0]/2,p[3].x[1]/2,p[3].x[2]/2,p[3].x[3]/2])
                    #pp=fv(fv.Add(p2,p4))
                    #p1=fv([p[0].x[0]/2,p[0].x[1]/2,p[0].x[2]/2,p[0].x[3]/2])
                    #p3=fv([p[2].x[0]/2,p[2].x[1]/2,p[2].x[2]/2,p[2].x[3]/2])
                    #k=fv(fv.Add(p1,p3))
                    F1=(CB.Ff_Dipole(Q2, 1.79, my, 2.79)[0])
                    F2=(CB.Ff_Dipole(Q2, 1.79, my, 2.79)[1])
                    k=1.79
                    G_m=F1+k*F2
                    G_e=F1-t*F2
                    RBb=(((1/137)**2)/(4*pm1.x[0]**2*((1-ctc)/2)**2))*(pm1.x[0]/pm3.x[0])*(((G_e**2+t*G_m**2)/(1+t))*((1+ctc)/2)+ 2*t*G_m**2*((1-ctc)/2))
                    RB.append(RBb*10**2)
                    #lamb=fv.Pro(pp,k)/my**2#((my**2+(Q/4))**0.5)*(Q/4)**0.5/(my**2)
                    #G_e=F1-(1.79)*t*F2
                    #G_m=F1+1.79*F2
                    Mm=((8*amp.e**4)/Q2**2)*((F1+k*F2)**2*((2*fv.Pro(p[0],p[1])*fv.Pro(p[2],p[1]))-fv.Pro(p[0],p[2])*(-fv.Pro(p[0],p[1])+fv.Pro(p[2],p[1])+my**2))-(((F1+k*F2)*F2)/2)*(4*(-my**2*fv.Pro(p[0],p[2])+2*fv.Pro(p[0],p[1])*fv.Pro(p[2],p[1])))+((k*F2)**2/(4*my**2))*(2*(fv.Pro(p[0],p[1])-fv.Pro(p[2],p[1])+2*my**2)*(-my**2*fv.Pro(p[0],p[2])+2*fv.Pro(p[0],p[1])*fv.Pro(p[2],p[1]))))                     
                    Mref.append(Mm*10**-4)#((amp.e**4/(t**2))*((((G_e**2+t*G_m**2)/(1+t))*(lamb**2-t**2-t))+2*t**2*G_m**2))
                    ctl=np.append(ctl,ctc) 
                    Iso=Gen.Iso_Weight(p[0], p[1], p[2], p[3], ctmin, ctmax)
                    ff=Gen.Flux_CMS(p[0],p[1],m1,m2)
                    Ps=Gen.PS(p[0], p[1], p[2], p[3],ctc)
                    a=amp.epAmp(pm1, pm2, pm3, pm4, cR, cL, cR_p, cL_p)[0]       
                    #print(F1)               
                    #print(F2)
                    sigma2.append(abs(a*Iso*ff)*10**2/(2*np.pi))#*p[2].x[0]/(p[0].x[0]**3))#              3.7*abs(Ps*a)
                    #rb=(amp.e**2/(4*p[0].x[0]**2*((1-ctc)/2)**2))*(p[2].x[0]/p[0].x[0])*(((F1**2-(1.79**2*Q2/(4*my**2))*F2**2)*(1+ctc)/2)-(Q2/(2*my**2))*(F1+1.79*F2)**2*(1-ctc)/2)
                    #RB.append(rb)
                    #E3.append(p[2].x[0])
                else:
                     i=0
            
            err=np.var(sigma2)*10**-7
            
            #RB=((1/137)**2*((1-(ctc))*t*G_m**2+(((1-(ctc))*(G_e**2+t*G_m))/(2*(1+t))))*p[1].x[0])/(p[0].x[0]**3*((1-(ctc))/2)**2)
            #RB=((1/137)**2/(4*p[0].x[0]**2*((1-ctc)/2)**2))*(p[2].x[0]/p[0].x[0])*(((G_e**2+t*G_m**2)/(1+t))*((1+ctc)/2)+2*t*G_m**2*((1-ctc)/2))
            #cross=np.append(cross,RB)
            #print(f*Iso/(1/(32*np.pi*my**2*p[0].x[0]**2))*(p[2].x[0]**2))
        #plt.hist(ctl)
        
        plt.xlabel('cos θ')
        plt.xlim(0,1)
        plt.ylim(0,1)
        plt.ylabel('dσ/dΩ [10–30 cm2 /sr]')
        plt.hist2d(ctl,sigma2)
        plt.colorbar()
        plt.yscale('log')
        plt.scatter(ctl,sigma2)
        plt.errorbar(ctl,sigma2,yerr=err,ls='none',color='magenta')
        plt.scatter(ctl,RB, color='white')
        # plt.yscale('log')   
        # plt.plot(Q,G_m,'s')
        # plt.show()
        # avg=np.mean(sigma)
        #err=np.std(sigma)/np.sqrt(N) 
        # plt.hist2d(ct,sigma,bins=[20,20])
        # plt.colorbar()
        # plt.xlabel("| costheta |")   
        # plt.ylabel("Differential Cross Section")
        # plt.show()
        # plt.xlabel('E [MeV]' )
        # plt.ylabel('Cross Section [1/MeV^2]')
        # plt.scatter(En,sigma)
        # plt.xlabel('E [MeV]' )
        # plt.ylabel('Cross Section [1/MeV^2]')
        # plt.scatter(En,cross)
        # RB_max=((1/137)**2*((1-max(ctl))*t*G_m**2+(((1-max(ctl))*(G_e**2+t*G_m))/(2*(1+t))))*p[1].x[0])/(p[0].x[0]**3*((1-max(ctl))/2)**2)
        # RB_min=((1/137)**2*((1-min(ctl))*t*G_m**2+(((1-min(ctl))*(G_e**2+t*G_m))/(2*(1+t))))*p[1].x[0])/(p[0].x[0]**3*((1-min(ctl))/2)**2)
        # sig_min=np.amin(sigma)
        # sig_max=np.amax(sigma)
        #E3=sum(E3)/len(ctl)
        # Q=sum(Q)/len(ctl)
        
        #print(Q)
        # Ge=sum(Ge)/len(ctl)
        # Gm=sum(Gm)/len(ctl)
        #E1=p[0].x[0]
        #F1=sum(F1)/len(ctl)
        #F2=sum(F2)/len(ctl)
        # Mref=sum(Mref)/len(ctl)
        # A=sum(A)/len(ctl)
        sigma2_tot=sum(sigma2)/len(ctl)
        #rb_tot=sum(RB)/len(ctl)
        
        Mref_tot=sum(Mref)/len(ctl)
        
        #ctlmax=ctmax#max(ctl)
        #ctlmin=ctmin#min(ctl)
        #k=1.79
       # A=((-k*Q/(2*my**2))*F2**2-2*F1**2)
        #B=(F1**2+((k*Q/(4*my**2))*F2**2)-(Q/(2*my**2))*(F1+k*F2)**2)
        ##rbmax=(1/137)**2*(E3/E1**3)*((A/(2*(1-ctlmax)))-(B*cmath.log(1-ctlmax))/2)
        #rbmin=(1/137)**2*(E3/E1**3)*((A/(2*(1-ctlmin)))-(B*cmath.log(1-ctlmin))/2)
        #print((E3/E1**3))
        #rbmax2=(E3/(75076*(-1+ctlmax)*E1**3*(1+t)))*((-3-2*ctlmax+ctlmax**2)*(Ge**2+Gm**2*t)+4*(-1+ctlmax)*(Ge**2-Gm**2*t**2)*cmath.log(-1+ctlmax))#(E3/(75076*E1**3))*((-4/(ctlmax-1))+ctlmax-4*(1+t)*cmath.log(ctlmax-1))
        #rbmin2=(E3/(75076*(-1+ctlmin)*E1**3*(1+t)))*((-3-2*ctlmin+ctlmin**2)*(Ge**2+Gm**2*t)+4*(-1+ctlmin)*(Ge**2-Gm**2*t**2)*cmath.log(-1+ctlmin))#(E3/(75076*E1**3))*((-4/(ctlmin-1))+ctlmin-4*(1+t)*cmath.log(ctlmin-1))
        
        # rbmax=2*np.pi*(-1*((2*(Ge**2 + Gm**2*t))/(-1 +ctlmax)) + ((2*E1*(Ge**2 + Gm**2*t) + m*(Ge**2 - Gm**2*t*(1+2*t)))*cmath.log(-1 + ctlmax))/m - ((2*E1*(Ge**2 + Gm**2*t) + m*(Ge**2 - Gm**2*t*(1 + 2*t)))*cmath.log((-1 + ctlmax)*E1 - m))/m)/(37538*E1**2*(1 + t))
        # rbmin=2*np.pi*(-1*((2*(Ge**2 + Gm**2*t))/(-1 +ctlmin)) + ((2*E1*(Ge**2 + Gm**2*t) + m*(Ge**2 - Gm**2*t*(1+2*t)))*cmath.log(-1 + ctlmin))/m - ((2*E1*(Ge**2 + Gm**2*t) + m*(Ge**2 - Gm**2*t*(1 + 2*t)))*cmath.log((-1 + ctlmin)*E1 - m))/m)/(37538*E1**2*(1 + t))
        #rb_tot=abs(rbmax-rbmin)
        return #sigma2_tot/Mref_tot#sigma2_tot/(rb_tot),sigma2_tot,rb_tot

print(MC_int_Weighted(E,ch,mx,my,m1,m2,ctmin,ctmax,cR,cL,cR_p,cL_p))

def MC_int_Unweighted(E,ch,mx,my,m1,m2,ctmin,ctmax,cR,cL,cR_p,cL_p):
        err=0
        sigma=[]
        sigma_tot=[]
        ct=[]
        N=10000
        p_acc=[]
        for i in range(N):
            p=Gen.Iso_Gen(E, mx, my, m1, m2, ctmin, ctmax)
            ct.append((np.dot(p[1].y,p[2].y)/(np.linalg.norm(p[1].y)*np.linalg.norm(p[2].y))))
            Iso=Gen.Iso_Weight(p[0],p[1],p[2],p[3],ctmin,ctmax)
            f=(Gen.Flux(p[0], p[1], m1, m2))
            a=amp.Z_new(p[0], p[1], p[2], p[3], cR, cL, cR_p, cL_p, ctmin, ctmax,ch)
            M=ME.MESqr(a,ch,p[0],p[1],p[2],p[3])
            sigma.append(abs((M*Iso*f)))            
        avg=np.mean(sigma)  
        sig_max=np.amax(sigma)
        for i in range(N):
            x=np.random.uniform(1,0)
            p=(sigma[i]/sig_max)
            if p<x:
                p_acc.append(p)
                sigma_tot.append(sigma[i]/N)
            else:
                p_acc.append(1-p)
                sigma_tot.append(0)
# =============================================================================
#         plt.hist2d(p_acc,sigma_tot,bins=[25,25])
#         plt.colorbar()
# =============================================================================
        err=np.var(p_acc)
        eff=np.average(p_acc)
        return err,eff
#print(MC_int_Unweighted(E,ch,mx,my,m1,m2,ctmin,ctmax,cR,cL,cR_p,cL_p))

#######Energy Dependence 
def cross(Ex):
    Output= np.ones((len(Ex)))
    for i in range(len(Ex)):
        E=Ex[i]
        Output[i]=((4*np.pi*((1/137)**2)))/(3*E**2)#*np.sqrt(1-(m1**2/(E/2)**2))*(1+(1/2)*(m1**2/(E/2)**2))
    
    return Output
def Energy(ch,cross,mx,my,m1,m2,ctmin,ctmax,cR,cL,cR_p,cL_p): 
        err=0
        sigma=[]
        ct=[]
        #cross=[]
        N=0
        En=[]
        error=[]
        while err==0 or err>1e-7:
            N=N+500
            for i in range(N):
                E=np.random.uniform(1,40)
                En.append(E)
                p=Gen.Iso_Gen(E, mx, my, m1, m2, ctmin, ctmax)
                #ct.append(abs(np.dot(p[1].y,p[2].y)/(np.linalg.norm(p[1].y)*np.linalg.norm(p[2].y))))
                Iso=Gen.Iso_Weight(p[0],p[1],p[2],p[3],ctmin,ctmax)
                f=(Gen.Flux_CMS(p[0], p[1], m1, m2))
                a=np.sum(amp.Z(p[0], p[1], p[2], p[3], cR, cL, cR_p, cL_p))
                A=a*np.conj(a)/4
                M=ME.MESqr(A,ch,p[0],p[1],p[2],p[3])
                sigma.append(abs((M*Iso*f*4*np.pi**2)))
                #cross.append(((4*np.pi*((1/137)**2))/(3*E**2))*np.sqrt(1-(m1**2/(E/2)**2))*(1+(1/2)*(m1**2/(E/2)**2)))
                error.append(np.var(sigma))
            err=np.var(sigma)
        #print(err)
        avg=np.mean(sigma)
        Ex=np.linspace(min(En),max(En),1000)
        #err=np.std(sigma)/np.sqrt(N) 
        # plt.hist2d(ct,sigma,bins=[20,20])
        # plt.colorbar()
        #plt.xlabel("|)
        # plt.show()
        for i in range(len(sigma)):
            sigma[i]=1.2*10**5*sigma[i]
        plt.xlabel('E_cm [GeV]')
        plt.ylabel('σ [nb]')
        plt.yscale('log')
        plt.ylim(10**-2,10)
        plt.plot(Ex,10**5*cross(Ex),color='red',alpha=1,label='Analytical')
        plt.plot(En,sigma,'s',alpha=0.6,marker='+',markersize=8,label='Monte Carlo')
        plt.legend()
        # plt.ylabel('Cross Section [1/MeV^2]')
        # plt.scatter(En,sigma)
        # plt.xlabel('E [MeV]' )
        # plt.plot(Ex,cross(Ex),color='red',alpha=0.6)
        # plt.scatter(En,sigma,marker='+')
        sig_min=np.amin(sigma)
        sig_max=np.amax(sigma)
        # cross=sum(cross)/N
        # sigma_tot=sum(sigma)/N
        return None


    
#print(Energy(ch,cross,mx,my,m1,m2,ctmin,ctmax,cR,cL,cR_p,cL_p))


# -*- coding: utf-8 -*-
"""
Created on Thu Jun 17 14:17:42 2021

@author: pa267340
"""

import numpy as np
from scipy.integrate import quad
from scipy.integrate import odeint
from scipy.special import gamma as Gammaf
import matplotlib.pyplot as plt
from matplotlib import cm
from constants import k,h,Te,cl,sT,Ne,me,qe,mu0, kTe, Tp, kTp
from plotter import setFigureParameters
from plotter import noLogsetFigureParameters
from scipy.special import kve

"""
beta = v/cl
gamma = 1/np.sqrt(1-beta**2)
Pe = gamma.me.nu
p = Pe/me/cl
"""   


"""
constants (like temperature or electron density) set in the constants.py file
"""

def generNu():
    """
    generation of the frequency mesh
    """
    global nu
    
    M=100
    numin=1e12
    numax=1e16
    nu = np.logspace(np.log10(numin),np.log10(numax),M)

def setParameters(B):
    """
    initialization of general parameters
    """
    
    #Larmor frequency
    nuL = qe*B/me/2/np.pi
    Ub=B**2/2/mu0
    theta = k*Te/me/cl**2
    K2 = kve(2, 1/theta)
    
    return nuL, Ub, theta, K2

    
def Ps(w,p):
    """
    function used in the calculation of j(nu,p)
    """

    fact=27*np.sqrt(3.0)/16.0/np.pi
    b=.02
    c=1.75
    d=3.45
    R0 = 2**(1./3)/5*Gammaf(1./3)**2
    g = np.sqrt(p*p+1)
    f = 2/(0.01*p+p**2) * (b + b*c*p + d*p*p/3.0) / ( 1.0 + c*p + d*p*p)
    j = np.zeros_like(w)
    X = f*(w-1./(g+p)) 
    j[X>0] = X[X>0]**(1./3.) / np.sqrt(1+(2*R0/np.pi)**2*X[X>0]**(2./3.)) * np.exp(-X[X>0])
    j[X>0] = f * R0 * p*p * j[X>0]
    j[X>0] = fact * j[X>0]
    
    return j

def j_nu_p(nu, Ub, nuL, p):
    """
    calculation of j(nu,p)
    """
    
    return (4/3) * cl * (sT*Ub / (4*np.pi*nuL)) * Ne * Ps(nu/nuL,p)



def f_p(p, K2, theta):
    """
    calculation of the maxwell-boltzmann distribution f(p,theta) defined by the 
    temperature theta
    """
    
    gamma = np.sqrt(p*p+1)
    
    return 1 / (theta * K2) * p*p * np.exp(-p*p/theta/(gamma+1))

        
def jMaxBol(xp, nu, Ub, nuL, K2, theta):
    """
    calculus of j(nu,p)*f(p,theta)
    """
    
    return j_nu_p(nu, Ub, nuL, xp) * f_p(xp, K2, theta)


def J_nu_theta(B, p, nuL, Ub, theta, K2):
    """
    Calculus and plotting of the synchrotron emission J(nu,theta) for a Maxwell-Boltzmann distribution of electrons
    Depends on the temperature theta which defines the profile of the distribution
    """
    
    # initialization of the array
    J_theta = np.zeros_like(nu)

    # for each frequency, evaluate J, the integral of j(nu,p)*f(p,theta) over p
    # as we can't take inf for the sup limit, we choose a mulitple of the characteristics border of the spectrum
    sup = 10*np.sqrt(theta *(1+theta))
    for i,nuu in enumerate(nu):
        J_theta[i], err = quad(jMaxBol,0, sup, args=(nuu, Ub, nuL, K2, theta))
        
    plt.plot(nu, J_theta)
    
    setFigureParameters(r"Tracé de l'émissivité synchrotron $J_\nu(\nu,\theta)$ pour un gaz d'électrons chauffé à {:} keV".format(kTe),
                        r'$J_\nu(\nu,\theta)$',r'$\nu$',5e5,5e8,1e2*nu[0],nu[-1])      
    
    return J_theta


def B_nu(theta):
    """
    Planck's black body law
    """
    
    B_nu = ((2*h)/cl**2) * nu**3 *  (1 / (np.exp((h*nu)/(k*Te)) - 1 )) 
    #print(B_nu)
    
    plt.plot(nu, B_nu)
    
    setFigureParameters(r"$B_\nu$".format(kTe)
                        ,r'$B_\nu(\theta)$',r'$\nu$',5e3,5e10,nu[0],nu[-1])      
    
    return B_nu


def alpha_nu_theta(B, p, nuL, Ub, theta, K2):
    """
    Synchrotron self absorption
    """
    
    J = J_nu_theta(B, p, nuL, Ub, theta, K2)
    B = B_nu(theta)
    
    alpha = J / B

    plt.plot(nu, B, label = r'$B_\nu(\theta)$')

    plt.plot(nu, J, label = r'$J_\nu(\nu,\theta)$')

    plt.plot(nu, alpha,label=r'$\alpha_\nu(\nu,\theta)$')
    
    #setFigureParameters(r"Tracé de l'absorption synchrotron $\alpha_\nu(\nu,\theta)$ pour un gaz d'électrons chauffé à {:} keV".format(kTe)
    #                    ,r'$\alpha_\nu(\nu,\theta)$',r'$\nu$',5e-5,5e30,nu[0],nu[-1])      
    
    setFigureParameters(r"Tracés de $J_\nu(\nu,\theta)$, $B_\nu(\theta)$ et $\alpha_\nu(\nu,\theta)$ pour un gaz d'électrons chauffé à {:} keV".format(kTe)
                        ,r'$f_\nu(\nu,\theta)$',r'$\nu$',5e-3,5e7,nu[0],nu[-1])      
    
    print(J)

    #print(alpha)
    return alpha 


def f(I, t, j, alpha):
    
    return cl * (j - alpha*I)
    

def init_planckian(Tp):
    """
    initialization of the photons distribution with a planckian profile
    """
    
    u0 = (1 /  (np.exp((h*nu)/(k*Te)) - 1 ))
    return u0

def solveTransfertEq(B, p, nuL, Ub, theta, K2):
    
    j = J_nu_theta(B, p, nuL, Ub, theta, K2)
    alpha = j / B_nu(theta)
        
    t = np.linspace(0,10,100)

    I0 = init_planckian(Te)
    I = odeint(f,I0,t, args = (j, alpha))


    print(I)

    plt.plot(t,I)
    setFigureParameters(r"Evolution de l'intensité spécifique de l'émission synchrotron pour un gaz d'électrons chauffé à {:} keV".format(kTe)
                        ,r'$I_\nu(t)$',r'$t$',5e0,5e30,0,10)      

#######################################################################
###################  TEST FUNCTIONS  #############################
######################################################################

def integJTh(p,Ub):
    """
    theoritical integration of j(nu,p)
    """
    #print(Ub)
    return 4/3*p*p*Ub*Ne*cl*sT/4/np.pi


def plotfp(p, B, K2, theta):
    """
    plotting of the maxwell-boltzmann distribution f(p, theta)
    """
           
    fp = np.zeros_like(p)
    for i, pp in enumerate(p):
        fp[i] = f_p(pp, K2, theta)
    
    plt.plot(p,fp,label = "theta = {:.2E}".format(theta))
    
    setFigureParameters("Gaz d'électrons chauffé à 100 kev",'f_p(p)','p',5e-4,5e2,1e-3,p[-1])

def testPs(nu,nuL):
    """
    tests of the Ps function
    """  

    p = [1e-3, 1e-2, 1.e-1, 1, 1.e1]
    col = [ cm.jet(x) for x in np.linspace(0, 0.3 , len(p))]
    w=nu/nuL

    for i, pp in enumerate(p):
        PS=Ps(w,p[i])
        gamma = np.sqrt(pp*pp+1)
        plt.plot(w,PS, color=col[i], label="p = {:.3f}".format(pp))
        print("p = ",p[i], " | p^2 = ",p[i]**2, " | integral_Ps = ", quad(Ps,1/(pp+gamma),10*3*gamma**2/2,args=(p[i],))[0])
        
    setFigureParameters('','Ps(w,p)','nu/nuL',1e-2,5e0,1e-2,3e3)
    
def testj_same_p(p, Ub, nuL):
    """
    calculus and plotting of the synchrotron radiation for an uniform distribution of electron
    (same p for all electrons)
    we use this function only for test purposes, so B is variating for tests
    """
    
    B=np.linspace(1e9,1e10,10)
    for b in B:
    
        #plus B est faible (respt fort), plus nu doit couvrir des freq basses (respt hautes)
        #B=1.e10
     
        gamma = np.sqrt(p*p+1)

        j=j_nu_p(nu, Ub, nuL, p)
        J = quad(j_nu_p,nuL/(p+gamma),nuL*10*3*gamma**2/2,args=(Ub,nuL,p))
        Jth = integJTh(p, Ub)
        
        print("Jth = {:.4E} (pour B = {:.1E})".format(Jth,b))
        print("J = {:.4E} (pour B = {:.1E})".format(J[0],b))
        #print(j)
        plt.plot(nu,j, label='B={:.1E} T'.format(b))

    setFigureParameters("Spectre de l'émissivité synchrotron pour p={:.1E} et Ne={:.1E}".format(p,Ne),'j(nu)','nu',5e3,5e7,nu[0],nu[-1])

def testJ_nu_theta(B, p, nuL, Ub, theta, K2):
    
    global temp
    temps = np.logspace(0, 5, 6)
    for temp in temps:
        # Te in K
        Te = temp * 1.602e-16 / k
        theta = k*Te/me/cl**2
        K2 = kve(2, 1/theta)

       # initialization of the array
        J_theta = np.zeros_like(nu)
    
        # for each frequency, evaluate J, the integral of j(nu,p)*f(p,theta) over p
        # as we can't take inf for the sup limit, we choose a mulitple of the characteristics border of the spectrum
        sup = 10*np.sqrt(theta *(1+theta))
        for i,nuu in enumerate(nu):
            J_theta[i], err = quad(jMaxBol,0, sup, args=(nuu, Ub, nuL, K2, theta))
            
        plt.plot(nu, J_theta, label="Temp = {:.0E} kev".format(temp))
        
    setFigureParameters(r"Tracé de l'émissivité synchrotron $J_\nu(\nu,\theta)$ pour un gaz d'électrons chauffé à différentes températures"
                        ,r'$J_\nu(\nu,\theta)$',r'$\nu$',5e5,5e8,nu[0],nu[-1])      
        

####################################################################
###################################################################


def synch():
    

    generNu()
    
    B=1e3
    p = np.linspace(1e-3, 1e1, 100)
    
    nuL, Ub, theta, K2 = setParameters(B)

    
    #plotfp(p, B, K2, theta)
    #testPs(nu,nuL)
    # tests pour différentes valeurs de p (le même pour tous les électrons)
    """
    p = [1e-2, 1.e-1, 1, 1.e1, 5e1]
    for pp in p:
        print("\n******  p = {:.1E}".format(pp))
        testj_same_p(pp, Ub, nuL)
    """    
    #testJ_nu_theta(B, p, nuL, Ub, theta, K2)
    #J_nu_theta(B, p, nuL, Ub, theta, K2)
    #alpha_nu_theta(B, p, nuL, Ub, theta, K2)
    solveTransfertEq(B, p, nuL, Ub, theta, K2)
    
    
synch()

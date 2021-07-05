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
from scipy.special import kve

from constants import k,h,Te,cl,sT,Ne,me,qe,mu0, kTe, B, R
from plotter import setFigureParameters
from plotter import plotIntensity
import changcooper as cc
from initDistrib import init_gaussian
from initDistrib import init_planckian

"""
beta = v/cl
gamma = 1/np.sqrt(1-beta**2)
Pe = gamma.me.nu
p = Pe/me/cl
"""   


"""
constants (like temperature or electron density) set in the constants.py file
"""

def setParameters(B):
    """
    initialization of general parameters
    """
    p = np.linspace(1e-3, 1e1, 100)

    #Larmor frequency
    nuL = qe*B/me/2/np.pi
    Ub=B**2/2/mu0
    theta = k*Te/me/cl**2
    K2 = kve(2, 1/theta)
    
    return p, nuL, Ub, theta, K2

    
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

def j_nu_p(nui, Ub, nuL, p):
    """
    calculation of j(nu,p) (for one specific frequency nui)
    """
    
    return (4/3) * cl * (sT*Ub / (4*np.pi*nuL)) * Ne * Ps(nui/nuL,p)



def f_p(p, K2, theta):
    """
    calculation of the maxwell-boltzmann distribution f(p,theta) defined by the 
    temperature theta
    """
    
    gamma = np.sqrt(p*p+1)
    
    return 1 / (theta * K2) * p*p * np.exp(-p*p/theta/(gamma+1))

        
def jManuol(xp, nui, Ub, nuL, K2, theta):
    """
    calculus of j(nu,p)*f(p,theta) (for one specific frequency nui)
    """
    
    return j_nu_p(nui, Ub, nuL, xp) * f_p(xp, K2, theta)


def J_nu_theta(nu, M, B, p, nuL, Ub, theta, K2):
    """
    Calculus and plotting of the synchrotron emission J(nu,theta) for a Maxwell-Boltzmann distribution of electrons
    Depends on the temperature theta which defines the profile of the distribution
    """
    
    # initialization of the array
    J_theta = np.zeros(M)

    # for each frequency, evaluate J, the integral of j(nu,p)*f(p,theta) over p
    # as we can't take inf for the sup limit, we choose a mulitple of the characteristics border of the spectrum
    sup = 100*np.sqrt(theta *(1+theta))
    for i,nuu in enumerate(nu[:-1]):
        J_theta[i], err = quad(jManuol,0, sup, args=(nuu, Ub, nuL, K2, theta))
        
    #plt.plot(nu, J_theta)
    
    #setFigureParameters(r"Tracé de l'émissivité synchrotron $J_\nu(\nu,\theta)$ pour un gaz d'électrons chauffé à {:} keV".format(kTe),
     #                   r'$J_\nu(\nu,\theta)$',r'$\nu$',5e5,5e8,1e2*nu[0],nu[-1])      
    
    return J_theta


def B_nu(nu, theta):
    """
    Planck's black body law
    """
    
    #B_nu = ((2*h)/cl**2) * nu**3 *  (1 / (np.exp((h*nu)/(k*Te)) - 1 )) 
    a = h*nu/k/Te
    """
    Bnu = np.zeros(len(nu))
    for i in range(len(nu)):
        
        if a[i] > 1e-8:
            Bnu[i] = ((2*h)/cl**2) * nu[i]**3 *  (np.exp(-a[i]) / (1 - np.exp(-a[i])))
        else:
            Bnu[i] = ((2*h)/cl**2) * nu[i]**2 *  k * Te
"""
    
    #plt.plot(nu, B_nu)
    # reshaping of the expression to avoid overflows
    Bnu = ((2*h)/cl**2) * nu**3 *  (np.exp(-a) / (1 - np.exp(-a)))

    #setFigureParameters(r"$B_\nu$".format(kTe)
     #                   ,r'$B_\nu(\theta)$',r'$\nu$',5e3,5e10,nu[0],nu[-1])      
    
    return Bnu



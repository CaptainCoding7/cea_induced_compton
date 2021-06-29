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
    sup = 10*np.sqrt(theta *(1+theta))
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
    Bnu = np.zeros(len(nu))
    for i in range(len(nu)):
        
        if a[i] > 1e-8:
            Bnu[i] = ((2*h)/cl**2) * nu[i]**3 *  (np.exp(-a[i]) / (1 - np.exp(-a[i])))
        else:
            Bnu[i] = ((2*h)/cl**2) * nu[i]**2 *  k * Te

    
    #plt.plot(nu, B_nu)
    # reshaping of the expression to avoid overflows
    Bnu = ((2*h)/cl**2) * nu**3 *  (np.exp(-a) / (1 - np.exp(-a)))

    #setFigureParameters(r"$B_\nu$".format(kTe)
     #                   ,r'$B_\nu(\theta)$',r'$\nu$',5e3,5e10,nu[0],nu[-1])      
    
    return Bnu[:-1]


def alpha_nu_theta(Jnu, Bnu):
    """
    Synchrotron self absorption
    """
    
    alpha = Jnu / Bnu
    
    """
    plt.plot(nu[1:], B, label = r'$B_\nu(\theta)$')

    plt.plot(nu[1:], J, label = r'$J_\nu(\nu,\theta)$')

    plt.plot(nu[1:], alpha,label=r'$\alpha_\nu(\nu,\theta)$')
    
    #setFigureParameters(r"Tracé de l'absorption synchrotron $\alpha_\nu(\nu,\theta)$ pour un gaz d'électrons chauffé à {:} keV".format(kTe)
    #                    ,r'$\alpha_\nu(\nu,\theta)$',r'$\nu$',5e-5,5e30,nu[0],nu[-1])      
    
    setFigureParameters(r"Tracés de $J_\nu(\nu,\theta)$, $B_\nu(\theta)$ et $\alpha_\nu(\nu,\theta)$ pour un gaz d'électrons chauffé à {:} keV".format(kTe)
                        ,r'$f_\nu(\nu,\theta)$',r'$\nu$',5e-3,5e7,nu[0],nu[-1])      
    """
    
    #print(alpha)
    return alpha 


def meshGeneration():
    """
    generation of the time and frequency meshes
    """
    
    #########################  MESH GENERATION  #######################
    
    #### time mesh  ################################
    
    to = ((k*Te) / (me*cl**2)) * Ne * cl * sT
    #to = cl / R

    # observation instants
    # the last instant is taken as the max time in the simulation
    tobs = [1e-7, 3e-7, 4e-7]
    tobs = [3e-9] # stabilité atteinte    
    
    tmax = tobs[-1]
    
    dt=1e-10
    dto = dt
    dto =  to * dt
    #print("dt;dto: ",dt,dto)
    N=int(tmax/dt)
    
    #### frequency mesh  ################################

    M=100
    numin=1e12
    numax=1e21
    #valeurs au bords des bins (M+1 valeurs)
    nu = np.logspace(np.log10(numin),np.log10(numax),M+1)
    
    return nu, tobs,tmax,dt,dto, M


def solveTransfertEq(nu,tobs,tmax,dt,dto,M,B, p, nuL, Ub, theta, K2):
    """
    solve the Transfert equation with the chang-cooper scheme
    call on of the plotting function
    """
        
    j = J_nu_theta(nu, M, B, p, nuL, Ub, theta, K2)
    # as the first values of j are 0, we give random values to avoid zero divisions
    #j[0:12] = 1e-3
    Bnu = B_nu(nu, theta)
    alpha = j / Bnu
    
    ################# TRANSFERT EQUATION  PARAMETERS  #################
    # we use the chang-cooper scheme without taking account of the A,B,C vectors
    
    # emission term
    emSyn = cl*j
    
    # absorption term
    absSyn = 1/cl/alpha
    
    # photons injection at 0.5 keV, gaussian width of 0.03keV
    xinj = 1e14
    width = 1e2
    I0=init_gaussian(nu, xinj, width)
    epho = nu * h / 1.602e-16
    #I0 = init_planckian(epho[1:])
    
    # find the solution of the transfert equation with the cc scheme for each instant
    # specified in tobs
    N=int(tobs[-1]/dt)
    I = cc.changCooper2(nu, N, I0,absSyn, emSyn, dt)

    #t=np.logspace(np.log10(1e-14),np.log10(1e-10),100)
    t = tobs[-1]
    Ith = I0 * np.exp(-cl*alpha*t) + Bnu*(1-np.exp(-cl*alpha*t))
    #print(alpha)
    
    intensity = ((2*h)/cl**2) * (I*nu[1:]**3) 

    plt.plot(nu[1:],Ith, '+', label = 'Ith')
    plt.plot(nu[1:],I0, label = 'I0')
    plt.plot(nu[1:],I, label = 'I_nu')
    plt.plot(nu[1:],Bnu, label = 'B_nu')
    setFigureParameters(r"Intensité spécifique de l'émission synchrotron pour un gaz d'électrons chauffé à {:} keV".format(kTe)
                        ,r'$I_\nu(t)$',r'$\nu$',1e-6,1e4,nu[0],nu[-1])      
    
    return I

def synch():
        
    p, nuL, Ub, theta, K2 = setParameters(B)

    
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
    
    nu,tobs,tmax,dt,dto, M = meshGeneration()
    #alpha_nu_theta(B, p, nuL, Ub, theta, K2)
    
    solveTransfertEq(nu, tobs, tmax, dt, dto, M, B, p, nuL, Ub, theta, K2)    
    
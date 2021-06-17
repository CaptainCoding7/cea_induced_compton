# -*- coding: utf-8 -*-
"""
Created on Thu Jun  3 13:44:28 2021

@author: Paul
"""

# Script to find a solution for the equation of Kompaneets with the Chang-Cooper scheme

import numpy as np
from scipy.stats import norm
from scipy.signal import unit_impulse
from scipy.stats import lognorm
from scipy.integrate import quad
from scipy.special import zeta
from scipy.special import gamma as Gammaf
import tridiagonal_matrix_solver as tms
from plotter import linplot
from plotter import logplot
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.optimize import fsolve



"""
x: array of the energy x=(h.nu)/(k.T) with nu the frequency
u: solutions of the equation for different instants

"""


def Ps(w,p):
    
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

def computej(w, Ub, nuL, p):
    
    return (4/3) * (sT*Ub / (4*np.pi*nuL)) * Ne * Ps(w,p)

   
def integJTh(p,Ub):
       print(Ub)
       return 4*p*p*Ub*Ne*cl*sT/4/np.pi


def synch():

    ## calculus and plotting of the synchrotron radiation    

    M=10000
    numin=1e18
    numax=1e24
    nu = np.logspace(np.log10(numin),np.log10(numax),M)
    #print("{:.1E}".format(nu[0]), nu[len(nu)-1])
    
    B=np.linspace(1e9,1e10,10)
    
    for b in B:
    
        #plus B est faible (respt fort), plus nu doit couvrir des freq basses (respt hautes)
        #B=1.e10
        
        """
        beta = v/cl
        gamma = 1/np.sqrt(1-beta**2)
        Pe = gamma.me.nu
        p = Pe/me/cl
        """
        p = 10
        
        #Larmor frequency
        nuL = qe*b/me/2/np.pi
        Ub=b**2/2/mu0
        w=nu/nuL
        
        j=computej(w, Ub, nuL, p)
        J = quad(computej,0,np.inf,args=(Ub,nuL,p))
        Jth = integJTh(p, Ub)
        
        print("Jth = {:.4E} (pour B = {:.1E})".format(Jth,b))
        print("J = {:.4E} (pour B = {:.1E})".format(J[0],b))
        #print(j)
        plt.plot(nu,j, label='B={:.1E} T'.format(b))
    
    setFigureParameters("Spectre de l'émissivité synchrotron pour p={:.1E} et Ne={:.1E}".format(p,Ne),'j(nu)','nu',5e-4,5e0,numin,numax)


def testPs(w):
        
    ############## tests fonction Ps    

    p = [1e-3, 1e-2, 1.e-1, 1, 1.e1]
    col = [ cm.jet(x) for x in np.linspace(0, 0.3 , len(p))]

    for i, pp in enumerate(p):
        PS=Ps(w,p[i])
        plt.plot(w,PS, color=col[i], label="p = {:.3f}".format(pp))
        print("p = ",p[i], " | p^2 = ",p[i]**2, " | integral_Ps = ", quad(Ps,0,np.inf,args=(p[i],))[0])
        
    setFigureParameters('Ps(w,p)','nu/nuL',1e-2,5e0,1e-2,3e3)
    
    

def init_gaussian(xinj, width):

    # on identifie le bin correspondant a xinj
    iinj = np.argwhere(xb<xinj)[-1]
    # on recentre xinj sur le bin en question
    xinj = xa[iinj]

    u0 = norm.pdf(xa , loc = xinj , scale = width )

    return u0


def init_dirac(xinj):
        
    #xinj = np.sqrt(nbPhotons)
    
    # on identifie le bin correspondant a xinj
    iinj = np.argwhere(xb<xinj)[-1]
    # on recentre xinj sur le bin en question
    xinj = xa[iinj]
    # on initialise avec des zeros partout sauf dans le bin en question
    # ou on met 1/dx de maniere avec avoir une integrale=1
    #u0 = np.zeros_like(xa)
    #u0[iinj] = 1.0/dxa[iinj]

    nbPhotons = xinj**2

    uinj = nbPhotons / ((xinj**2) * dxa[iinj])
    
    u0 = uinj * unit_impulse(len(xa), int(iinj))

    return u0

def init_planckian(Tp):
    
    #u0 = (1 /  (np.exp((h*f_photons)/(k*Tp)) - 1 ))
    u0 = (1 /  (np.exp(e_photons/kTp) - 1 ))
    return u0


def setFigureParameters(title, ylabel, xlabel, ymin, ymax,  xmin, xmax):
    
    # Reglages affichage
    plt.xscale('log')
    plt.yscale('log')
    leg = plt.legend(loc="upper right",prop={'size': 7}, bbox_to_anchor=[1, 1],
                     ncol=1, shadow=True, title="Legend", fancybox=True)
    leg.get_title().set_color("black")
    plt.title(title)
    plt.ylim(ymin,ymax)
    plt.xlim(xmin,xmax)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.show()
    
def noLogsetFigureParameters(title, ylabel, xlabel, ymin, ymax,  xmin, xmax):
    
    # Reglages affichage
    #plt.xscale('log')
    leg = plt.legend(loc="upper right",prop={'size': 7}, bbox_to_anchor=[1, 1],
                     ncol=1, shadow=True, title="Legend", fancybox=True)
    leg.get_title().set_color("black")
    plt.ylim(ymin,ymax)
    plt.xlim(xmin,xmax)
     #   plt.xticks([1, 2, 3, 4, 5])
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.show()


# solve the equation with the chang-cooper scheme
def changCooper(pech, Q, dT):
            
    #print(dT)
    
    #dT=dt

    for n in range(N-1):
        
       #we redefine B as it depends on u    
        B = 0.5*(xa[1:]**4+xa[:-1]**4) * (u[n][1:]+1) # sur les bords      (M-1)
        w = (B/C)*dxb
        lwl = np.abs(w)
        W = lwl * np.exp(-lwl)/(1.0-np.exp(-lwl))
        Wp = W*np.exp(+w/2)
        Wm = W*np.exp(-w/2)
        
        
        # sous-diagonale (taille M-1)
        a = 1/A[1:]/dxa[1:] * C/dxb * Wm
    
        # sur-diagonale (taille M-1)
        c = 1/A[:-1]/dxa[:-1] * C/dxb * Wp
    
        # diagonale (taille M)
        b = np.zeros_like(xa)
        b[0]    = 1/dT + 1/A[   0]/dxa[   0] * (                         0 + C[0 ]/dxb[ 0]*Wm[ 0] ) + pech[0]
        b[1:-1] = 1/dT + 1/A[1:-1]/dxa[1:-1] * ( C[:-1]/dxb[:-1] * Wp[:-1] + C[1:]/dxb[1:]*Wm[1:] ) + pech[1:-1]
        b[-1]   = 1/dT + 1/A[  -1]/dxa[  -1] * ( C[ -1]/dxb[ -1] * Wp[ -1] +                    0 ) + pech[-1]
        
        r = Q + u[n]/dT
                
        u[n+1] = tms.tridag(-a,b,-c,r)
        
        #print("---- u = ",u)
            
    return u


def plotOccupationRate():
    
    print("Plotting occupation rate...")
        
    # plot the occupation rate
    for i,tt in enumerate(tobs):
        print("Tracé ",i+1)
        index = int(tt / dt) - 1
        plt.plot(xa,u[index],color=col[i],label='t={:.2f}s '.format(index*tmax/5))    
    
    setFigureParameters('u(x,t)','x')
    
    
def plotDensityFreq():
    
    print("Plotting spectral density...")
                  
    phoDensTh = ((8*np.pi)/cl**3) * f_photons**2 / (findCst()*np.exp((h*f_photons)/(k*Te)) - 1 ) 
    plt.plot(f_photons, phoDensTh, "+",color='red',label='theoritical solution')
    
    # Calculus of the final density (after 30s)
    # Must be constant (same as the initial one)
    Ndens0 = 8*np.pi * ((k*Te)/(cl*h))**3 * np.sum(u[N-1]*dxa*xa**2)
    #print("Densité finale: ",Ndens)

    for i,tt in enumerate(tobs):
        print("Tracé ",i+1)
        index = int(tt / dt) - 1
        if index==-1:
            index=0
        
        # plot the photon density rate
        phoDens = ((8*np.pi)/cl**3) * (u[index]*f_photons**2)

        # multiply by Enorm as we represent the energy in abscisse
        plt.plot(f_photons, phoDens, color = col[i],label='t={:.2f}s '.format(tt))
    
    setFigureParameters('Spectral Density $ (nbPhotons . m^{-3} . Hz^{-1})$','Frequency (Hz)',1e2,1e13,1e13,1e21)
    
    
def plotDensity():
    
    print("Plotting density...")
                 
    phoDensTh = ((8*np.pi)/cl**3) * f_photons**2 * (1 / (findCst()*np.exp((h*f_photons)/(k*Te)) - 1 )) *1.602e-16 / h 
    plt.plot(e_photons, phoDensTh, "+",color='red',label='theoritical solution')

    # Calculus of the final density (after 30s)
    # Must be constant (same as the initial one)
    Ndens = np.sum(u[N-1]*xa*xa*dxa)
    #print("Densité finale: ",Ndens)

    for i,tt in enumerate(tobs):
        print("Tracé ",i+1)
        index = int(tt / dt) - 1
        if index==-1:
            index=0
        print(index)
        
        # plot the photon density rate
        """
        two to convert the density from nbPho/m3/Hz into nbPho/m3/kev
        - 1kev = 2.4e17
        - E=h.nu ==> divide by h (and  multiply by 1.602e-16)
        """
        #phoDens = ((8*np.pi)/cl**3) * (u[index]*(f_photons**2)) * 2.4e17
        phoDens = ((8*np.pi)/cl**3) * (u[index]*(f_photons**2)*1.602e-16) / h 
        #phoDens = u[index]*e_photons**2
        plt.plot(e_photons, phoDens, color = col[i],label='t={:.1E}s'.format(tt))
    
    setFigureParameters('Photon Number Density $ (nbPhotons . m^{-3} . keV^{-1})$','Energy (keV)',1e20,1e30,1e-4,1e3)


def plotIntensity():
    
    print("Plotting intensity...")
    
    #C=findCst()
    #intensityTh = ((2*h)/cl**2) * f_photons**3 * ( (1 /  (C*np.exp((h*f_photons)/(k*Te)) - 1 ))) 
    #plt.plot(e_photons, intensityTh, "+", color='red',label='theoritical solution')
    
    # Calculus of the final density (after 30s)
    # Must be constant (same as the initial one)
    Ndens = 8*np.pi * ((k*Te)/(cl*h)) * np.sum(u[N-1]*xa**2*dxa)
    print("Densité finale: ",Ndens)

    for i,tt in enumerate(tobs):
        print("Tracé ",i+1)
        index = int(tt / dt) - 1
        if index==-1:
            index=0
        print(index)

        # plot the intensity of the photon field
        intensity = ((2*h)/cl**2) * (u[index]*f_photons**3)
        
        # multiply by Enorm as we represent the energy in abscisse
        plt.plot(e_photons, intensity, color = col[i],label='t={:.1E}s'.format(tt))
    

    setFigureParameters('Spectral Radiance $(keV.m^{-2}.s^{-1}.Hz^{-1}.str^{-1})$','Energy (keV)',1e-1,5e7,1e-4,1e3)


def plotEnergyDensity():
    
    print("Plotting Energy density... ")
    
    t = np.linspace(dt, tmax/5., 25., endpoint=True)
    Er = np.zeros(25)
    i=0
    for tt in t:
        index = int(tt / dt) - 1
        print(index)
        # Calculus of the photon energy density
        Er[i] = ( ((2*(k*Te)**4) / ((h**3)*(cl**3) )) * np.sum(u[index]*4*np.pi*xa**3))
        #Er[i] = (8*np.pi**5*k**4/(15*h**3*cl**3))*Te**4
        print(Er[i])
        i+=1
    Ee = (3/2)*Ne*k*np.full_like(Er,Te)
    print(Ee)
    plt.plot(t, Ee,label='Ee')    
    plt.plot(t, Er,label='Er')
        
    setFigureParameters('Photon Energy Density $ (erg . m^{-3})$','Time (s)',1e12, 1e21,0,tmax/5)
    
    
def BEDistrib(x,c):
    """
    return the Bose-Einstein solution from x and c
    """
    return x**2 / (c*np.exp(x)-1)
    

def F(c):
    
    """
    the equation used to find the constant
    """
    eq = quad(BEDistrib, 0, 100, args=(c,)) - ((cl*h)**3*Ndens0) / (8*np.pi*(k*Te)**3)
    return eq[0]


def findCst():
    """
    Function used to find the cst of the Bose-Einstein solution
    """    
    root = fsolve(F,2)
    print("cst = ",root[0])
    return root[0]

####################################################

def findNmax():
    
    Nmax = 8*np.pi*(((k*Te)/(cl*h))**3)*2*zeta(3)
    return Nmax

def defineConstants():
    
    global k,h,Te,cl,sT,Ne,me,qe,mu0,tobs, kTe, Tp, kTp

    #################  CONSTANTS  #################
    
    # Boltzman constant
    k = 1.381e-23
    # Planck constant
    h = 6.626e-34
    # Temperature of the electron field (K)
    # kTe in keV 
    kTe = 1.
    # Te in K
    Te = kTe * 1.602e-16 / k
    cl = 3e8
    # Thompson scattering cross section
    sT = 6.652e-29
    # Electron number density (in m^-3)
    Ne = 1e23
    # Electron mass (in kg)
    me = 9.109e-31
    # Electron charge (in C)
    qe = 1.602e-19
    #permeabilité magnetique du vide
    mu0 = np.pi*4e-7
    
    
    # Photon temperature
    kTp = 2.
    Tp = kTp * 1.602e-16 / k
    

def meshGeneration():
    
    global xb,xa,dxa,dxb,tobs,tmax,dt,dto,N,M, f_photons, e_photons

    #########################  MESH GENERATION  #######################
    
    #### time mesh  ################################
    
    nu = ((k*Te) / (me*cl**2)) * Ne * cl * sT
    
    # observation instants
    # the last instant is taken as the max time in the simulation
    #tobs = 5e-9 * np.arange(12)
    #tobs = np.arange(0,26e-8,5e-8)
    tobs = [0, 1e-7, 3e-7, 4e-7]
    tobs = [0, 2e1] # stabilité atteinte    
    
    tmax= tobs[len(tobs)-1]
    
    dt=1e-2
    dto =  nu * dt
    #print("dt;dto: ",dt,dto)
    N=int(tmax/dt)
    
    #### "energy" mesh  ################################

    M=1000
    
    # the energy carried by the photons (keV)
    emin = 1e-3
    emax = 1e3
    e_photons = np.linspace(emin,emax,M)   
    # frequency vector
    # *1.602e-16 keV ==> J
    f_photons = e_photons * 1.602e-16 / h
    
    # x = h*f / k*Te 
    # with f the frequency 
    # x is a no-dimension quantity
    xmin = emin / kTe
    xmax = emax / kTe
    
    # Valeurs sur les bords des bins (M+1 values)
    #xb = np.logspace(np.log10(xmin),np.log10(xmax),M+1) 
    xb = np.linspace(xmin,xmax,M+1) 
    
    # Largeur des bins (M valeurs)
    dxa = xb[1:] - xb[:-1]
    # valeurs au centre des bins (M valeurs)
    xa = 0.5 * ( xb[1:] + xb[:-1] )
    # distance entre les milieux (M-1 valeurs)
    dxb = xa[1:] - xa[:-1]
    
    
def solveKompaneets():
    
    global A,C,Q,pech,col, u, Ndens0
    
    ################# KOMPANEETS EQUATION  PARAMETERS  #################
    
    # Initialization of the A, C, Q
    # the vectors containing the terms from the Kompaneets equation
    # B is defined in the changCooper function as it evolves with the time
    
    A = xa*xa 
    C =  0.5*(xa[1:]**4+xa[:-1]**4) # sur les bords      (M-1)
    
    # injection term
    Q = np.full_like(xa, 0)            # au milieu          (M)
    #Q=init_cond(xinj)
    
    # sink term
    #n*ro = n/T(x) => T(x) = 1/ro
    ro2=0
    pech = np.full_like(xa, ro2)            # au milieu          (M)
    
    ### initialization of the solution vector : [time]["energy"] (x is a pseudo-energy)
    
    u=np.empty([N,M])
    # photons injection at 0.5 keV, gaussian width of 0.03keV
    einj = 0.5
    ewidth = 0.03
    xinj = einj / kTe
    width = ewidth / kTe
    print(xinj, width)
    #u[0]=init_gaussian(xinj, width)
    u[0] = init_planckian(Tp)
    
    # Cette fois, on impose le nombre de photons
    # STABLE EN DECA DE ~2.4
    #nbPhotons = 1
    
    Ndens0 = 8*np.pi * ((k*Te)/(cl*h))**3 * np.sum(u[0]*dxa*xa**2)
    print("Densité initiale: ",Ndens0)
    
    
    ############  SOLVING KOMPANEETS EQUATION AND PLOTTING FUNCTIONS  ###########
    
    # update the vector u with the solution of the Kompaneets equation for each instant
    changCooper(pech,Q,dto)
    
    # Plotting parameter
    col = [ cm.jet(x) for x in np.linspace(0, 0.3 , len(tobs))]
    
    #plotDensityFreq()
    #plotDensity()
    #plotOccupationRate()
    plotIntensity()
    #plotEnergyDensity()


def main():
    
    defineConstants()
    meshGeneration()
    
    #solveKompaneets()  
    
    synch()
    
    #Nmax=findNmax()
    #print(Nmax)
    

if __name__ == "__main__":
    main()
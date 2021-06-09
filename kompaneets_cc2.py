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
from math import pow
from math import exp
from math import sinh
import tridiagonal_matrix_solver as tms
from plotter import linplot
from plotter import logplot
from scipy.special import iv as Inu
import matplotlib.pyplot as plt
from matplotlib import cm

"""
x: array of the energy x=(h.nu)/(k.T) with nu the frequency
u: solutions of the equation for different instants

"""


def synch(M,x):

    return x*0

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


def setFigureParameters(ylabel, xlabel, ymin, ymax,  xmin, xmax):
    
    # Reglages affichage
    plt.xscale('log')
    plt.yscale('log')
    leg = plt.legend(loc="upper right",prop={'size': 8}, bbox_to_anchor=[1, 1],
                     ncol=1, shadow=True, title="Legend", fancybox=True)
    leg.get_title().set_color("black")
    plt.ylim(ymin,ymax)
    plt.xlim(xmin,xmax)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
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
        
    # plot the occupation rate
    for i,tt in enumerate(tobs):
        print("Tracé ",i+1)
        index = int(tt / dt) - 1
        plt.plot(xa,u[index],color=col[i],label='t={:.2f}s '.format(index*tmax/5))    
    
    setFigureParameters('u(x,t)','x')
    
    
def plotDensity():
                  
    phoDensTh =  (1 / (np.exp((h*f_photons)/(k*Te)) - 1 )) * e_photons**2
    #plt.plot(e_photons, phoDensTh, ":",label='theoritical solution')
    
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
        #phoDens = ((8*np.pi)/cl**3) * (u[index]*f_photons**3) / e_photons
        phoDens = u[index]*e_photons**2
        # multiply by Enorm as we represent the energy in abscisse
        plt.plot(e_photons, phoDens, color = col[i],label='t={:.2f}s '.format(tt))
    
    setFigureParameters('Photon Number Density $ (cm^{-3} . keV^{-1})$','Energy (keV)',1e-6,1e2,1e-4,1e4)

def plotDensityFreq():
                  
    phoDensTh = ((8*np.pi)/cl**3) * ( (1 / (np.exp((h*f_photons)/(k*Te)) - 1 ))   *f_photons**2)
    #plt.plot(f_photons, phoDensTh, ":",label='theoritical solution')

    
    # Calculus of the final density (after 30s)
    # Must be constant (same as the initial one)
    Ndens = np.sum(u[N-1]*xa*xa*dxa)
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
    
    setFigureParameters('Spectral Density $ (cm^{-3} . Hz^{-1})$','Frequency (Hz)',1e2,1e15,1e12,1e21)

def plotIntensity():
 
    intensityTh = ((2*h)/cl**2)  * ( (1 / (np.exp((h*f_photons)/(k*Te)) - 1 ))   *f_photons**3)
    #plt.plot(e_photons, intensityTh, ":",label='theoritical solution')


    for i,tt in enumerate(tobs):
        print("Tracé ",i+1)
        index = int(tt / dt) - 1
        if index==-1:
            index=0
        print(index)

        # plot the photon density rate
        intensity = ((2*h)/cl**2) * (u[index]*f_photons**3)

        # multiply by Enorm as we represent the energy in abscisse
        
        plt.plot(e_photons, intensity, color = col[i],label='t={:.2f}s '.format(tt))
    
    setFigureParameters('Spectral Radiance','Energy (keV)',1e-2,1e4,1e-2,1e3)


def plotEnergyDensity():
    
    t = np.linspace(dt, tmax/5, 25, endpoint=True)
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
        
    setFigureParameters('Photon Energy Density $ (erg . cm^{-3})$','Time (s)',1e12, 1e21,0,tmax/5)
    

####################################################

def main():

    global xb,xa,dxa,dxb,k,h,Te,cl,Ne,tobs,tmax,dt,dto,N,M,u,A,C,Enorm,col, f_photons, e_photons

    #################  CONSTANTS  #################
    
    # Boltzman constant
    k = 1.381e-23
    # Planck constant
    h = 6.626e-34
    # Temperature of the electron field (K)
    # 1keV ~ k*1.1e7 
    Te = 1e8
    cl= 3e8
    # Thompson scattering cross section
    sT = 6.652e-29
    # Electron number density
    Ne=6e22
    # Electron mass
    me=9.109e-31
    
    ##############  MESH GENERATION  ###############
    
    #### time mesh
    
    # observation instants
    # the last instant is takenn as the max time in the simulation
    #tobs = [0.0,0.05,0.1,0.15,0.2,3]
    #tobs = [0,0.1,0.15,0.2,0.25,0.3,0.35,0.4,10]
    #tobs = [0,0.005,0.01,0.02,0.03,0.04,0.05,0.06,0.07]
    #tobs = [0.07,0.1,0.13,0.16,0.19]
    #tobs = [1,2,3,4,5,6,7,8,9]
    #tobs = [0,5e-9,5e-8,1e-7]
    tobs = [0,10]
    tmax= tobs[len(tobs)-1]
    
    dt=5e-3
    dto = ((k*Te) / (me*cl**2)) * Ne * cl * sT * dt
    
    print("dt;dto: ",dt,dto)
    N=int(tmax/dt)
    
    #### energy mesh
    
    M=1000
    
    # the energy carried by the photons (keV)
    emin = 1.e-3
    emax = 1.e3
    #e_photons = np.logspace(np.log10(emin),np.log10(emax),M+1) 
    # car 1ev = 1.602e-19 J ==> 1 kev = 1.602e-16 J
    Enorm = k * Te / 1.602e-16
    #print(Enorm)
    
    
    # x = h*f / k*Te 
    # with f the frequency 
    # x is a no-dimension quantity
    xmin = emin / Enorm
    xmax = emax / Enorm
    #print(xmax)
    
    #xmin = emin
    #xmax = emax
    
    # Valeurs sur les bords des bins (M+1 values)
    xb = np.logspace(np.log10(xmin),np.log10(xmax),M+1) 
    #xb = e_photons/ Enorm
    # Largeur des bins (M valeurs)
    dxa = xb[1:] - xb[:-1]
    # valeurs au centre des bins (M valeurs)
    xa = 0.5 * ( xb[1:] + xb[:-1] )
    # distance entre les milieux (M-1 valeurs)
    dxb = xa[1:] - xa[:-1]
    
    e_photons = xa * Enorm
    # *1.602e-16 keV ==> J
    f_photons = e_photons * 1.602e-16 / h
    

    #################  MAIN & PARAMETERS  #################
    
    # initialization of the solution vector
    
    u=np.empty([N,M])
    # photons injection at 0.5 keV, gaussian width of 0.03keV
    xinj = 0.5 / Enorm
    width = 0.03 / Enorm
    print(xinj, width)
    u[0] = init_gaussian(xinj, width)
    
    # Cette fois, on impose le nombre de photons
    # STABLE EN DECA DE ~2.4
    #nbPhotons = 1
    
    # Calculus of the initial density
    # simple integration
    Ndens = np.sum(u[0]*xa*xa*dxa)
    #print("Densité initiale: ",Ndens)
    
    # Initialization of the A, B, C, Q, w and W vectors
    # the vectors containing the terms from the Kompaneets equation
    
    A = xa*xa 
    C =  0.5*(xa[1:]**4+xa[:-1]**4) # sur les bords      (M-1)
    
    # injection term
    Q = np.full_like(xa, 0)            # au milieu          (M)
    #Q=init_cond(xinj)
    
    # sink term
    #n*ro = n/T(x) => T(x) = 1/ro
    ro2=0
    pech = np.full_like(xa, ro2)            # au milieu          (M)
    
    
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


if __name__ == "__main__":
    main()
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

"""
x: array of the energy x=(h.nu)/(k.T) with nu the frequency
u: solutions of the equation for different instants

"""


# Solution Analytique (Park-Petrosian 1995)
def uth(x,t,xinj):
    q = 3.0
    theta = 1.0
    s = 0.0
    a = 1.0
    alpha = (2-q-s)/2.0
    nu = np.abs((q-1+a)/2/alpha)
    G = 1/np.abs(alpha)/2/t * x**((1-q+a)/2)*xinj**((1-q-a)/2)
    G *= Inu(nu,(x*xinj)**alpha/2/t/alpha**2)
    G *= np.exp(-(x**(2*alpha)+xinj**(2*alpha))/4/alpha**2/t-theta*t)
    return G


def synch(M,x):

    return x*0

def init_gaussian(xinj):

    # on identifie le bin correspondant a xinj
    iinj = np.argwhere(xb<xinj)[-1]
    # on recentre xinj sur le bin en question
    xinj = xa[iinj]

    u0 = norm.pdf(xa , loc = xinj , scale = 0.03 )

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


def setFigureParameters(u,ylabel):
    
    
    # Reglages affichage
    plt.xscale('log')
    plt.yscale('log')
    leg = plt.legend(loc="upper right",prop={'size': 8}, bbox_to_anchor=[1, 1],
                     ncol=1, shadow=True, title="Legend", fancybox=True)
    leg.get_title().set_color("black")
    plt.ylim(1.e-6,1.e2)
    plt.xlim(1.e-4,1.e4)
    plt.xlabel('x')
    plt.ylabel(ylabel)
    plt.show()
    

# solve the equation with the chang-cooper scheme
def changCooper(u0, N, pech, Q, dT):
    
    global B
    
    u=u0
    
    print(dT)

    # simulation for different instants
    for n in range(N):
        
       #we redefine B as it depends on u    
        B = 0.5*(xa[1:]**4+xa[:-1]**4) * (u[1:]+1) # sur les bords      (M-1)
        w = (B/C)*dxb
        lwl = np.abs(w)
        W = lwl * np.exp(-lwl)/(1.0-np.exp(-lwl))
        Wp = W*np.exp(+w/2)
        Wm = W*np.exp(-w/2)
        
        
        # sous-diagonale (taille M-1)
        a = 1/A[1:]/dxa[1:] * C/dxb * Wm
    
        # sur-diagonale (taille M-1)
        c = 1/A[:-1]/dxa[:-1] * C/dxb * Wp
    
        # doagonale (taille M)
        b = np.zeros_like(xa)
        b[0]    = 1/dT + 1/A[   0]/dxa[   0] * (                         0 + C[0 ]/dxb[ 0]*Wm[ 0] ) + pech[0]
        b[1:-1] = 1/dT + 1/A[1:-1]/dxa[1:-1] * ( C[:-1]/dxb[:-1] * Wp[:-1] + C[1:]/dxb[1:]*Wm[1:] ) + pech[1:-1]
        b[-1]   = 1/dT + 1/A[  -1]/dxa[  -1] * ( C[ -1]/dxb[ -1] * Wp[ -1] +                    0 ) + pech[-1]
        
        r = Q + u/dT
        
        u = tms.tridag(-a,b,-c,r)
            
    return u

    
def solveKompaneets_ech_compare():
        
    plt.plot(xa,u0,color='purple',label='t=0s')

    for i,tt in enumerate(tobs):
        n = int(tt/dt)
        u = changCooper(u0,n,np.zeros_like(xa),np.zeros_like(xa))
        uech = changCooper(u0,n,pech,np.full_like(xa, 0))
        
        strlabel1='t={:.2f}s '.format(tt),'avec echappt' 
        strlabel2='t={:.2f}s '.format(tt),'sans echappt' 
        
        plt.plot(xa,norm[i]*uech,color=col[i],label=strlabel1)
        plt.plot(xa,norm[i]*u,":",color=col[i],label=strlabel2 )
    
    setFigureParameters(u)
    
def solveKompaneets_inj_compare():
        
    plt.plot(xa,u0,color='purple',label='t=0s')
    
    for i,tt in enumerate(tobs):
        n = int(tt/dt)
        u = changCooper(u0,n,np.zeros_like(xa),np.zeros_like(xa))
        uinj = changCooper(u0,n,np.zeros_like(xa),Q)

        strlabel1='t={:.2f}s '.format(tt),'avec injection' 
        strlabel2='t={:.2f}s '.format(tt),'sans injection' 
        
        plt.plot(xa,norm[i]*uinj,color=col[i],label=strlabel1)
        plt.plot(xa,norm[i]*u,":",color=col[i],label=strlabel2 )
    
    setFigureParameters(u)


def solveKompaneets():
        
    #plt.plot(xa,u0,color='purple',label='t=0s')
    
    phoDens0 = u0*xa**2
    plt.plot(xa,phoDens0,color='purple',label='t=0s')
    
    
    for i,tt in enumerate(tobs):
        dto = ((k*Te) / (me*cl**2)) * Ne * cl * sT * dt
        n = int(tt/dto)
        u = changCooper(u0,n,pech,Q,dto)
        # Caluculus of the final density (after 30s)
        # Must be constant (same as the initial one)
        Ndens = np.sum(u*xa*xa*dxa)
        print("Densité finale: ",Ndens)

        # plot the occupation rate
        #plt.plot(xa,norm[i]*u,color=col[i],label='t={:.2f}s '.format(tt))
    
        # plot the photon density rate
        phoDens = u*xa**2
        plt.plot(xa,norm[i]*phoDens,color=col[i],label='t={:.2f}s '.format(tt))
    
    #setFigureParameters(u,'u(x,t)')
    setFigureParameters(u,'Photon Number Density')

####################################################
## Mesh generation

# time
#Nmax=30.01
# observation instants
tobs = [0.2,2.,20.]
dt=5e-3


#N=int(Nmax/dt)

# "energy"
M=1000
xmin = 1.e-3
xmax = 1.e3

# Boltzman constant
k = 1.381e-23
# Temperature of the electron field (K)
# 1keV ~ k1.1e7 
Te = 1.1e7
cl= 3e8
# Thompson scattering cross section
sT = 6.652e-29
# Electron number density
Ne=6e22
# Electron mass
me=9.109e-31

#x=np.logspace(-3, 3, M, endpoint=True)


# Valeurs sur les bords des bins (M+1 values)
xb = np.logspace(np.log10(xmin),np.log10(xmax),M+1)
# Largeur des bins (M valeurs)
dxa = xb[1:] - xb[:-1]
# valeurs au centre des bins (M valeurs)
xa = 0.5 * ( xb[1:] + xb[:-1] )
# distance entre les milieux (M-1 valeurs)
dxb = xa[1:] - xa[:-1]


# initialization of the solution vector

#u=np.empty([N,M])
# STABLE EN DECA DE 1.513 (inclus)
xinj=0.5

# Cette fois, on impose le nombre de photons
# STABLE EN DECA DE ~2.4
#nbPhotons = 1

u0 = init_gaussian(xinj)

# Calculus of the initial density
# simple integration
Ndens = np.sum(u0*xa*xa*dxa)
print("Densité initiale: ",Ndens)

# Initialization of the A, B, C, Q, w and W vectors
# the vectors containing the values obtained with the specific function 
# in the equation (34)

A = xa*xa 
C =  0.5*(xa[1:]**4+xa[:-1]**4) # sur les bords      (M-1)

# injection term
Q = np.full_like(xa, 0)            # au milieu          (M)
#Q=init_cond(xinj)

# sink term
#n*ro = n/T(x) => T(x) = 1/ro
ro2=0
pech = np.full_like(xa, ro2)            # au milieu          (M)
symb2=""

col=['r','g','b']
norm = [1,1,1]
fig = plt.figure()

# test comparison sink/no sink
#solveKompaneets_ech_compare()

#test comparison injection/no injection
#solveKompaneets_inj_compare()

# tout parametre compris
solveKompaneets()

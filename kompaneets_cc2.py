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


def init_cond(nbPhotons):
    
    xinj = np.sqrt(nbPhotons)
    print(xinj)
    
    # on identifie le bin correspondant a xinj
    iinj = np.argwhere(xb<xinj)[-1]
    # on recentre xinj sur le bin en question
    xinj = xa[iinj]
    # on initialise avec des zeros partout sauf dans le bin en question
    # ou on met 1/dx de maniere avec avoir une integrale=1
    u0 = np.zeros_like(xa)
    u0[iinj] = 1.0/dxa[iinj]

       
    u0[iinj] = nbPhotons / ((xinj**2) * dxa[iinj])
    print(u0[iinj])
    print("dxa[inj]=",dxa[iinj])
    
    return u0


def setFigureParameters(u):
    
    # Caluculus of the final density (after 30s)
    # Must be constant (same as the initial one)
    Ndens = np.sum(u*xa*xa*dxa)
    print("Densité finale: ",Ndens)
    
    # Reglages affichage
    plt.xscale('log')
    plt.yscale('log')
    leg = plt.legend(loc="upper right",prop={'size': 8}, bbox_to_anchor=[1, 1],
                     ncol=1, shadow=True, title="Legend", fancybox=True)
    leg.get_title().set_color("black")
    plt.ylim(1.e-8,1.e2)
    plt.xlim(1.e-4,1.e4)
    plt.xlabel('x')
    plt.ylabel('u(x,t)')
    plt.show()
    

# solve the equation with the chang-cooper scheme
def changCooper(u0, N,pech,Q):
    
    global B
    
    u=u0

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
        b[0]    = 1/dt + 1/A[   0]/dxa[   0] * (                         0 + C[0 ]/dxb[ 0]*Wm[ 0] ) + pech[0]
        b[1:-1] = 1/dt + 1/A[1:-1]/dxa[1:-1] * ( C[:-1]/dxb[:-1] * Wp[:-1] + C[1:]/dxb[1:]*Wm[1:] ) + pech[1:-1]
        b[-1]   = 1/dt + 1/A[  -1]/dxa[  -1] * ( C[ -1]/dxb[ -1] * Wp[ -1] +                    0 ) + pech[-1]
        
        r = Q + u/dt
        
        u = tms.tridag(-a,b,-c,r)
            
    return u


def solveKompaneets_ech_compare():
        
    plt.plot(xa,u0,color='purple',label='t=0s')
    
    # plotting for 3 instants
    t = [0.3,3.,30.]
    col=['r','g','b']
    norm = [1,1,1]
    for i,tt in enumerate(t):
        n = int(tt/dt)
        u = changCooper(u0,n,np.zeros_like(xa),np.zeros_like(xa))
        uech = changCooper(u0,n,pech,np.full_like(xa, 0))
        
        #print(u)
        #logplot(dt, xa, norm[i]*u, n)
        #logplot(dt,xa, norm[i]*uth(xa,n*dt, 0.1),tt )
        strlabel1='t={:.2f}s '.format(tt),'avec echappt' 
        strlabel2='t={:.2f}s '.format(tt),'sans echappt' 
        
        plt.plot(xa,norm[i]*uech,color=col[i],label=strlabel1)
        plt.plot(xa,norm[i]*u,":",color=col[i],label=strlabel2 )
    
    setFigureParameters(u)
    
def solveKompaneets_inj_compare():
        
    plt.plot(xa,u0,color='purple',label='t=0s')
    
    # plotting for 3 instants
    t = [0.3,3.,30.]
    col=['r','g','b']
    norm = [1,1,1]
    for i,tt in enumerate(t):
        n = int(tt/dt)
        u = changCooper(u0,n,np.zeros_like(xa),np.zeros_like(xa))
        uinj = changCooper(u0,n,np.zeros_like(xa),Q)
        
        #print(u)
        #logplot(dt, xa, norm[i]*u, n)
        #logplot(dt,xa, norm[i]*uth(xa,n*dt, 0.1),tt )
        strlabel1='t={:.2f}s '.format(tt),'avec injection' 
        strlabel2='t={:.2f}s '.format(tt),'sans injection' 
        
        plt.plot(xa,norm[i]*uinj,color=col[i],label=strlabel1)
        plt.plot(xa,norm[i]*u,":",color=col[i],label=strlabel2 )
    
    setFigureParameters(u)


def solveKompaneets():
        
    plt.plot(xa,u0,color='purple',label='t=0s')
    
    # plotting for 3 instants
    t = [0.3,3.,30.]
    col=['r','g','b']
    norm = [1,1,1]
    for i,tt in enumerate(t):
        n = int(tt/dt)
        u = changCooper(u0,n,pech,Q)
        
        #print(u)
        #logplot(dt, xa, norm[i]*u, n)
        #logplot(dt,xa, norm[i]*uth(xa,n*dt, 0.1),tt )
        strlabel1='t={:.2f}s '.format(tt),'avec injection' 
        strlabel2='t={:.2f}s '.format(tt),'sans injection' 
        
        plt.plot(xa,norm[i]*u,color=col[i],label='t={:.2f}s '.format(tt))
    
    setFigureParameters(u)

    

####################################################
## Mesh generation

# time
#Nmax=30.01
t = [0.3,3.,30.]
dt=0.01
#N=int(Nmax/dt)

# energy
M=100
xmin = 1.e-3
xmax = 1.e3

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
#xinj=0.1

# Cette fois, on impose le nombre de photons
# STABLE EN DECA DE ~2.4
nbPhotons = 1
u0 = init_cond(nbPhotons)

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

fig = plt.figure()

# test comparison sink/no sink
#solveKompaneets_ech_compare()

#test comparison injection/no injection
#solveKompaneets_inj_compare()

# tout parametre compris
solveKompaneets()

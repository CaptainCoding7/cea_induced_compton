# -*- coding: utf-8 -*-
"""
Created on Thu Jun  3 13:44:28 2021

@author: Paul
"""

# Script to find a solution for the equation of Kompaneets with the Chang-Cooper scheme

import numpy as np

from constants import k,h,Te,cl,sT,Ne,me,qe,mu0,kTe, Tp, kTp
import changcooper as cc
from plotter import setFigureParameters
from plotter import noLogsetFigureParameters
from plotter import plotIntensity
from thSol import findCst
from thSol import findNmax
from initDistrib import init_gaussian
from initDistrib import init_planckian

"""
x: array of the energy x=(h.nu)/(k.T) with nu the frequency
u: solutions of the equation for different instants

"""



def meshGeneration():
    """
    generation of the time and energy meshes
    """
    
    global N, M, f_photons, e_photons

    #########################  MESH GENERATION  #######################
    
    #### time mesh  ################################
    
    nu = ((k*Te) / (me*cl**2)) * Ne * cl * sT
    
    # observation instants
    # the last instant is taken as the max time in the simulation
    #tobs = 5e-9 * np.arange(12)
    #tobs = np.arange(0,26e-8,5e-8)
    tobs = [1e-7, 3e-7, 4e-7]
    tobs = [2e1] # stabilité atteinte    
    
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
    
    return xa,xb,dxa,dxb,tobs,tmax,dt,dto
    
    
def solveKompaneets(xa,xb,dxa,dxb,tobs,tmax,dt,dto):
    """
    solve the Kompaneets equation with the chang-cooper scheme
    call on of the plotting function
    """
    
    ################# KOMPANEETS EQUATION  PARAMETERS  #################
    
    # Initialization of the A, C, Q
    # the vectors containing the terms from the Kompaneets equation
    # B is defined in the changCooper function as it evolves with the time
    
    A = xa*xa 
    C =  0.5*(xa[1:]**4+xa[:-1]**4) # sur les bords      (M-1)
    
    # injection term
    Q = np.full_like(xa, 0)            # au milieu          (M)
    
    # sink term
    #n*ro = n/T(x) => T(x) = 1/ro
    ro2=0
    pech = np.full_like(xa, ro2)            # au milieu          (M)
    
    ### initialization of the solution vector : [time]["energy"] (x is a pseudo-energy)
    
    #u=np.empty([N,M])
    
    # photons injection at 0.5 keV, gaussian width of 0.03keV
    einj = 0.5
    ewidth = 0.03
    xinj = einj / kTe
    width = ewidth / kTe
    print(xinj, width)
    u0 = init_planckian(e_photons)
    
    # Cette fois, on impose le nombre de photons
    # STABLE EN DECA DE ~2.4
    #nbPhotons = 1
    
    Ndens0 = 8*np.pi * ((k*Te)/(cl*h))**3 * np.sum(u0*dxa*xa**2)
    print("Densité initiale: ",Ndens0)
    
    
    ############  SOLVING KOMPANEETS EQUATION AND PLOTTING FUNCTIONS  ###########
        
    # find the solution of the kompaneets equation with the cc scheme for each instant
    # specified in tobs
    uobs = np.empty([len(tobs)+1,M])
    uobs[0] = u0
    print(len(uobs[0]))
    for i,tt in enumerate(tobs):    
        N=int(tmax/dt)
        uobs[i+1] = cc.changCooper(A, C, xa, dxa, dxb, N, u0, pech, Q, dto)
    
    # Calculus of the final density (after 30s)
    # Must be constant (same as the initial one)
    umax=uobs[len(uobs)-1]
    Ndens = 8*np.pi * ((k*Te)/(cl*h)) * np.sum(umax*xa**2*dxa)
    print("Densité finale: ",Ndens)
    
    C=findCst(Ndens0)
    intensityTh = ((2*h)/cl**2) * f_photons**3 * ( (1 /  (C*np.exp((h*f_photons)/(k*Te)) - 1 ))) 
    
    plotIntensity(Ndens0, uobs, tobs, e_photons, f_photons)
    
    phoDensTh = ((8*np.pi)/cl**3) * f_photons**2 / (findCst(Ndens0)*np.exp((h*f_photons)/(k*Te)) - 1 ) 
    #plotDensityFreq(xa, dxa, tobs, dt)
    
    #phoDensTh = ((8*np.pi)/cl**3) * f_photons**2 * (1 / (findCst(Ndens0)*np.exp((h*f_photons)/(k*Te)) - 1 )) *1.602e-16 / h 
    #plotDensity(xa, dxa, tobs, dt)
    
    #plotOccupationRate(xa, dxa, tobs, dt)


def mainKompaneets():
    
    xa,xb,dxa,dxb,tobs,tmax,dt,dto = meshGeneration()

    solveKompaneets(xa,xb,dxa,dxb,tobs,tmax,dt,dto)  
    
    Nmax=findNmax()
    
    
    #print(Nmax)
    

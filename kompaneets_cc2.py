# -*- coding: utf-8 -*-
"""
Created on Thu Jun  3 13:44:28 2021

@author: Paul
"""

# Script to find a solution for the equation of Kompaneets with the Chang-Cooper scheme

import numpy as np
import matplotlib.pyplot as plt

from constants import k,h,Te,cl,sT,Ne,me,kTe, B, Bcgs, R, Rcgs
import changcooper as cc
from plotter import plotIntensity
from plotter import plotAll
from plotter import setFigureParameters
from thSol import findCst
from thSol import findNmax
from initDistrib import init_gaussian
from initDistrib import init_planckian
from synchrotron import J_nu_theta
from synchrotron import B_nu
from synchrotron import setParameters
from compareBELM import compareBELM
from compareBELM import compareBELMnoCompton

"""
x: array of the energy x=(h.nu)/(k.T) with nu the frequency
u: solutions of the equation for different instants

"""

def meshGeneration():
    """
    generation of the time and energy meshes
    """
    
    #########################  MESH GENERATION  #######################
    
    #### time mesh  ################################
    
    # pseudo-time
    #to = ((k*Te) / (me*cl**2)) * Ne * cl * sT
    to = cl / R
    
    # observation instants
    # the last instant is taken as the max time in the simulation   
    tobs = [1.e-3, 1.e-2, 5e-2,1e-1]
    tobs = [1.e-1]
    dt=1e-4
    dto =  to * dt
    #print("dt;dto: ",dt,dto)
    
    #### "energy" mesh  ################################

    M=100
    
    # the energy carried by the photons (keV)
    emin = 5e-6
    emax = 5e3
    e_pho = np.logspace(np.log10(emin),np.log10(emax),M+1)   
    
    # x = h*f / k*Te 
    # with f the frequency 
    # x is a no-dimension quantity
    xmin = emin / kTe
    xmax = emax / kTe
    
    # Valeurs sur les bords des bins (M+1 values)
    #xb = np.logspace(np.log10(xmin),np.log10(xmax),M+1) 
    xb = np.logspace(np.log10(xmin),np.log10(xmax),M+1) 
    
    # Largeur des bins (M valeurs)
    dxa = xb[1:] - xb[:-1]
    # valeurs au centre des bins (M valeurs)
    xa = 0.5 * ( xb[1:] + xb[:-1] )
    # distance entre les milieux (M-1 valeurs)
    dxb = xa[1:] - xa[:-1]
    
    # frequency vector
    nu = k*Te*xa/h
    print(nu)
    
    return xa,xb,dxa,dxb,tobs,dt,dto, M, e_pho, nu


def initKompParam(xa, nu, e_pho, M):
    """    
    ################# KOMPANEETS EQUATION  PARAMETERS  #################
    
    # Initialization of the A, C, Q parameters
    # the vectors containing the terms from the Kompaneets equation
    # B is defined in the changCooper function as it evolves with the time
    """
    
    # profondeur optique de Thompson
    pT = Ne * R * sT
    A = xa**2 * me * cl**2 / pT / k  / Te
    #C =  0.5*(xa[1:]**4+xa[:-1]**4) # sur les bords      (M-1)    
    C =  (0.5*(xa[1:]+xa[:-1]))**4 # sur les bords      (M-1)
    # ---------------- Synchrotron ----------
        
    # synchrotron emmisivity
    p, nuL, Ub, theta, K2 = setParameters(B)
    # calculus in function of x
    Jnu = J_nu_theta(nu, M, B, p, nuL, Ub, theta, K2)
    
    # verif Jnu
    #plt.plot(e_pho[1:], Jnu)
    #setFigureParameters(r"Tracé de l'émissivité synchrotron $J_\nu(\nu,\theta)$ pour un gaz d'électrons chauffé à {:} keV".format(kTe),
    #                    r'$J_\nu(\nu,\theta)$',r'$\nu$',1e-10,1e1,e_pho[1],e_pho[-1])          
    
    Q = Jnu * cl**2 * R / 2 / h / nu**3
    
    # synchrotron absorption / sink term
    Bnu = B_nu(nu, theta)
    alpha = Jnu / Bnu
    # profondeur optique
    pOpt = alpha*R
    # probabilité d'échappement
    pech = 3/4
    # pOpt*n + pech*n = n/T(x) => T(x) = 1/(pech+pOpt)    
    # synchrotron absorption AND sink term gathered in T
    T = 1 / (pech + pOpt)    
    
    return A, C, Q, T, Jnu, Bnu, pOpt, pech
    

def showDensities(uobs, xa, dxa):
    """
    Compute and print the initial and final densities 
    of the photons distribution.
    """
    Ndens0 = 8*np.pi * ((k*Te)/(cl*h))**3 * np.sum(uobs[0]*dxa*xa**2)
    print("Densité initiale: ",Ndens0)
    
    umax=uobs[len(uobs)-1]
    Ndens = 8*np.pi * ((k*Te)/(cl*h)) * np.sum(umax*xa**2*dxa)
    print("Densité finale: ",Ndens)
    
    return Ndens0, Ndens

def computeLum(u, nu, p_nu):
    
    # intensité specifique
    I_nu = 2*h/cl**2 * u*nu**3 
    #print(I_nu)
    # densité d'énergie en photons
    u_nu = 4 * np.pi / cl * I_nu
    # volume de la source
    V = 4/3 * np.pi * R**3
    # temps moyen d'échappement
    T_nu = R/cl/p_nu
    # énergie total de la source
    E_nu = V*u_nu
    # luminosité de la source
    L_nu = E_nu / T_nu
    
    return I_nu, L_nu


def solveKompaneets(xa,xb,dxa,dxb,tobs,dt,dto, M, e_pho, nu):
    """
    solve the Kompaneets equation with the chang-cooper scheme
    call of the plotting function
    """
    ####### initialization of the parameters of the Kompaneets equation  ######
    A, C, Q, T, Jnu, Bnu, pOpt, pech = initKompParam(xa, nu, e_pho, M)
    e_pho = e_pho[1:]

    ####### initialization of the photons distribution ######
    # photons injection at 0.5 keV, gaussian width of 0.03keV
    einj = 0.5
    ewidth = 0.03
    xinj = einj / kTe
    width = ewidth / kTe
    u0 = init_planckian(e_pho)
    #u0 = init_gaussian(xa,xb, xinj, width)
        
    ############  SOLVING KOMPANEETS EQUATION AND PLOTTING FUNCTIONS  ###########
        
    # find the solution of the kompaneets equation with the cc scheme for each instant
    # specified in tobs
    #------------------------------------------------------------------------#
    u, a, b, c, r = cc.changCooper(A, C, T, Q, xa, dxa, dxb, u0, M, tobs, dt, dto) #
    ueq = u[-1]
    #------------------------------------------------------------------------#

    # the initial and final densities of the photons distribution
    Ndens0, Ndens = showDensities(u, xa, dxa)    
    # theoritical solution of the kompaneets equation WITHOUT the synchrotron emission/absorption
    #C=findCst(Ndens0)
    #intensityTh = ((2*h)/cl**2) * f_photons[1:]**3 * ( (1 /  (C*np.exp((h*f_photons[1:])/(k*Te)) - 1 ))) 
     
    # synchrotron contribution ==> au regime stable du/dto = 0
    uS = Q / (1 + pOpt)

    #compton contribution
    m = np.diag(-a, -1) + np.diag(b, 0) + np.diag(-c, 1)
    I = np.identity(len(m))+np.identity(len(m))+np.identity(len(m))
    m = (I - m)/dto + np.diag(pOpt) + pech*I
    #m = (m-I)/dto - pOpt - pech*I
    #uC = ueq * m[50]
    uC = np.matmul(m, ueq)
    print(uC)
    
    toErg = 1e7
    
    L = np.empty([len(tobs)+1,M])
    I = np.empty([len(tobs)+1,M])
    for i in range(len(u)):
        I[i],L[i] = computeLum(u[i]*toErg,nu, pech)
    Is, Ls = computeLum(uS*toErg, nu, pech)
    Ic, Lc = computeLum(uC*toErg, nu, pech)
    
    """
    plotAll(Bnu, nu*Ls, nu*Lc, nu*L, tobs, e_pho, nu,
                  'Spectre de la couronne (boule ; B = {:.0E} G, kT = {:} keV, R = {:.0E} cm)'.format(Bcgs,kTe,Rcgs),
                  r'$\nu.L_\nu$ $(erg.s^{-1})$','Energy (keV)',1e28,5e32,1e-5,1e4)
    """
    L_Bnu = 4 * np.pi**2 * R**2  * Bnu * toErg
    compareBELM(nu*L_Bnu,nu*Ls, nu*Lc, nu*L, tobs, e_pho, nu,
                 'Spectre de la couronne (boule ; B = {:.0E} G, kT = {:} keV, R = {:.0E} cm)'.format(Bcgs,kTe,Rcgs),
                  r'$\nu.L_\nu$ $(erg.s^{-1})$','Energy (keV)',1e28,5e38,1e-10,1e4)


def mainKompaneets():
    
    xa,xb,dxa,dxb,tobs,dt,dto, M, e_pho, nu = meshGeneration()

    solveKompaneets(xa, xb, dxa, dxb, tobs, dt, dto, M, e_pho, nu)
    

mainKompaneets()
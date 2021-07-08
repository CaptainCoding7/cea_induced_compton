# -*- coding: utf-8 -*-
"""
Created on Thu Jun  3 13:44:28 2021

@author: Paul
"""

# Script to find a solution for the equation of Kompaneets with the Chang-Cooper scheme

import numpy as np
import matplotlib.pyplot as plt

from constants import k,h,cl,me
import changcooper as cc
from plotter import plotIntensity
from plotter import plotAll
from plotter import setFigureParameters
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

    M=300
    
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
    
    return xa,xb,dxa,dxb,tobs,dt,dto, M, e_pho, nu


def initKompParam(xa, nu, e_pho, M):
    """    
    ################# KOMPANEETS EQUATION  PARAMETERS  #################
    
    # Initialization of the A, C, Q parameters
    # the vectors containing the terms from the Kompaneets equation
    # B is defined in the changCooper function as it evolves with the time
    """

    A = xa**2 * me * cl**2 / pT / k  / Te
    #C =  0.5*(xa[1:]**4+xa[:-1]**4) # sur les bords      (M-1)    
    C =  (0.5*(xa[1:]+xa[:-1]))**4 # sur les bords      (M-1)
   
    # ---------------- Synchrotron ----------
        
    # synchrotron emmisivity
    p, nuL, Ub, theta, K2 = setParameters(B, Te)
    
    # calculus in function of x
    Jnu = J_nu_theta(nu, M, B, p, nuL, Ub, theta, K2, Ne)
    
    Q = Jnu * cl**2 * R / 2 / h / nu**3
    
    # synchrotron absorption / sink term
    Bnu = B_nu(nu, Te)
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

def computeLum(u, nu):
    
    # intensité specifique
    I_nu = 2*h/cl**2 * u*nu**3 
    #print(I_nu)
    # densité d'énergie en photons
    u_nu = 4 * np.pi / cl * I_nu
    # volume de la source
    V = 4/3 * np.pi * R**3
    # temps d'échappement T_nu = R / c
    T_nu = R/cl
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
    dtoinv=1/dto
    #dtoinv=0
    # with the induced compton effect:
    u_ind, a, b, c, r = cc.changCooper(A, C, T, Q, xa, dxa, dxb, u0, M, tobs, dt, dto, True) #
    # without the induced compton effect
    u, a, b, c, r = cc.changCooper(A, C, T, Q, xa, dxa, dxb, u0, M, tobs, dt, dto, False) #
    u_ind_eq = u_ind[-1]
    ueq = u[-1]
    #------------------------------------------------------------------------#

    # the initial and final densities of the photons distribution
    Ndens0, Ndens = showDensities(u, xa, dxa)    
    
    ### The compton and synchrotron contribution are calculated for the 
    #   solution found with the induced compton effect

    # synchrotron contribution ==> au regime stable du/dto = 0
    uS = Q / (pOpt - pech)
    uS = (Q - pOpt*u_ind_eq)
    
    #compton contribution
    m = np.diag(-a, -1) + np.diag(b, 0) + np.diag(-c, 1)
    I = np.identity(len(m))
    m = I*dtoinv - m + np.diag(pOpt+pech)
    uC = np.matmul(m, u_ind_eq)
    
    toErg = 1e7
    
    ### Calculus of the intensities/luminosities  ###
    I_ind,L_ind = computeLum(pech*u_ind_eq*toErg, nu)
    I, L = computeLum(pech*ueq*toErg, nu)
    Is, Ls = computeLum(uS*toErg, nu)
    Ic, Lc = computeLum(uC*toErg, nu)
    L_Bnu = 4 * np.pi**2 * R**2  * Bnu * toErg
    
    # plot the luminosities*nu
    
    plotAll(nu*L_Bnu, nu*Ls, nu*Lc, nu*L, nu*L_ind, tobs, e_pho, nu,
                  r'Spectre de la couronne (boule ; B = {:.0E} G, kT = {:} keV, R = {:.0E} cm, $\tau_p$ = {:.2E})'.format(Bcgs,kTe,Rcgs,pT),
                  r'$\nu.L_\nu$ $(erg.s^{-1})$','Energy (keV)',ymin,ymax,xmin,xmax)
    
    """
    # same plots with curves from BELM Model
    compareBELM(nu*L_Bnu,nu*Ls, nu*Lc, nu*L, tobs, e_pho, nu,
                 r'Spectre de la couronne (boule ; B = {:.0E} G, kT = {:} keV, R = {:.0E} cm, $\tau_p$ = {:.2E})'.format(Bcgs,kTe,Rcgs,pT),
                  r'$\nu.L_\nu$ $(erg.s^{-1})$','Energy (keV)',ymin,ymax,xmin,xmax)
    """
    
def mainKompaneetsFromGui(r,b,kT,pt,xMin,xMax,yMin,yMax):
    
    import constants as cst
    import plotter as pr

    cst.Rcgs = r
    cst.R = cst.Rcgs*1e-2
    cst.Bcgs = b
    cst.B = cst.Bcgs * 1e-4
    cst.kTe = kT
    cst.Te = cst.kTe * 1.602e-16 / cst.k
    cst.pT = pt
    cst.Ne = cst.pT / cst.R / cst.sT

    pr.xmin = xMin
    pr.xmax = xMax
    pr.ymin = yMin
    pr.ymax = yMax

    from constants import Te, kTe, B, Bcgs, R, Rcgs, pT, Ne
    from plotter import xmin,xmax,ymin,ymax
    global Te, kTe, B, Bcgs, R, Rcgs, pT, Ne,xmin,xmax,ymin,ymax
    
    
    xa,xb,dxa,dxb,tobs,dt,dto, M, e_pho, nu = meshGeneration()

    solveKompaneets(xa, xb, dxa, dxb, tobs, dt, dto, M, e_pho, nu)
    
def mainKompaneets():
    
    global xmin,xmax,ymin,ymax
    
    xmin = 1e-5
    xmax=1e4
    ymin=1e28
    ymax=5e33    
    from constants import Te, kTe, B, Bcgs, R, Rcgs, pT, Ne
    global Te, kTe, B, Bcgs, R, Rcgs, pT, Ne
    
    xa,xb,dxa,dxb,tobs,dt,dto, M, e_pho, nu = meshGeneration()

    solveKompaneets(xa, xb, dxa, dxb, tobs, dt, dto, M, e_pho, nu)
    
#mainKompaneets()
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
from plotter import plotIntensityAll
from plotter import setFigureParameters
from thSol import findCst
from thSol import findNmax
from initDistrib import init_gaussian
from initDistrib import init_planckian
from synchrotron import J_nu_theta
from synchrotron import B_nu
from synchrotron import setParameters

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
    dt=1e-5
    dto =  to * dt
    #print("dt;dto: ",dt,dto)
    
    #### "energy" mesh  ################################

    M=100
    
    # the energy carried by the photons (keV)
    emin = 5e-9
    emax = 5e3
    e_photons = np.logspace(np.log10(emin),np.log10(emax),M+1)   
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
    xb = np.logspace(np.log10(xmin),np.log10(xmax),M+1) 
    
    # Largeur des bins (M valeurs)
    dxa = xb[1:] - xb[:-1]
    # valeurs au centre des bins (M valeurs)
    xa = 0.5 * ( xb[1:] + xb[:-1] )
    # distance entre les milieux (M-1 valeurs)
    dxb = xa[1:] - xa[:-1]
    
    return xa,xb,dxa,dxb,tobs,dt,dto, M, e_photons, f_photons


def initKompParam(xa, f_photons, e_photons, M):
    
    ################# KOMPANEETS EQUATION  PARAMETERS  #################
    
    # Initialization of the A, C, Q
    # the vectors containing the terms from the Kompaneets equation
    # B is defined in the changCooper function as it evolves with the time
    
    # profondeur optique de Thompson
    pT = Ne * R * sT
    A = xa**2 * me * cl**2 / pT / k /Te
    C =  0.5*(xa[1:]**4+xa[:-1]**4) # sur les bords      (M-1)
    
    # ---------------- Synchrotron ----------
        
    # synchrotron emmisivity
    p, nuL, Ub, theta, K2 = setParameters(B)
    # calculus in function of x
    Jnu = J_nu_theta(f_photons, M, B, p, nuL, Ub, theta, K2)
    
    # verif Jnu
    #plt.plot(e_photons[1:], Jnu)
    #setFigureParameters(r"Tracé de l'émissivité synchrotron $J_\nu(\nu,\theta)$ pour un gaz d'électrons chauffé à {:} keV".format(kTe),
    #                    r'$J_\nu(\nu,\theta)$',r'$\nu$',1e-10,1e1,e_photons[1],e_photons[-1])          
    
    Q = Jnu * cl**2 * R / 2 / h / f_photons[1:]**3
    
    # synchrotron absorption / sink term
    Bnu = B_nu(f_photons, theta)
    alpha = Jnu / Bnu
    # profondeur optique
    pOpt = alpha*R
    # pOpt*n + n = n/T(x) => T(x) = 1/(1+pOpt)    
    # synchrotron absorption AND sink term gathered in T
    T = 1 / (1 + pOpt)    
    
    return A, C, Q, T, Jnu, Bnu, pOpt
    

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


def solveKompaneets(xa,xb,dxa,dxb,tobs,dt,dto, M, e_photons, f_photons):
    """
    solve the Kompaneets equation with the chang-cooper scheme
    call of the plotting function
    """
    ####### initialization of the parameters of the Kompaneets equation  ######
    A, C, Q, T, Jnu, Bnu, pOpt = initKompParam(xa, f_photons, e_photons, M)

    ####### initialization of the photons distribution ######
    # photons injection at 0.5 keV, gaussian width of 0.03keV
    einj = 0.5
    ewidth = 0.03
    xinj = einj / kTe
    width = ewidth / kTe
    #u0 = init_planckian(e_photons[1:])
    u0 = init_gaussian(xb, xinj, width)
        
    ############  SOLVING KOMPANEETS EQUATION AND PLOTTING FUNCTIONS  ###########
        
    # find the solution of the kompaneets equation with the cc scheme for each instant
    # specified in tobs
    #------------------------------------------------------------------------#
    uobs = cc.changCooper(A, C, T, Q, xa, dxa, dxb, u0, M, tobs, dt, dto) #
    #------------------------------------------------------------------------#

    # the initial and final densities of the photons distribution
    Ndens0, Ndens = showDensities(uobs, xa, dxa)    
    # theoritical solution of the kompaneets equation WITHOUT the synchrotron emission/absorption
    #C=findCst(Ndens0)
    #intensityTh = ((2*h)/cl**2) * f_photons[1:]**3 * ( (1 /  (C*np.exp((h*f_photons[1:])/(k*Te)) - 1 ))) 
    
    # synchrotron only
    N=int(tobs[-1]/dt)
    uS = cc.changCooper2(xb, N, u0, T, Q, dto)
    #compton only
    T= np.zeros_like(xa)
    Q = np.zeros_like(xa)
    uobsComp = cc.changCooper(A, C, T, Q, xa, dxa, dxb, u0, M, tobs, dt, dto)
    
    plotIntensityAll(Jnu, Bnu, uS, uobsComp, uobs, tobs, e_photons[1:], f_photons[1:],
                  'Intensité de la couronne (B = {:.0E} G, kT = {:} keV, R = {:.0E} cm)'.format(Bcgs,kTe,Rcgs),'Spectral Radiance $(J.m^{-2}.s^{-1}.Hz^{-1}.str^{-1})$',
                  'Energy (keV)',1e-15,5e10,1e-10,1e10)
    
    #plotIntensity(0, uobs, tobs, e_photons[1:], f_photons[1:],
    #             'Intensité de la couronne (B = {:.0E} G, kT = {:} keV, R = {:.0E} cm)'.format(Bcgs,kTe,Rcgs),'Spectral Radiance $(keV.m^{-2}.s^{-1}.Hz^{-1}.str^{-1})$',
    #             'Energy (keV)',1e-15,5e10,1e-10,1e10)
    
    #phoDensFreqTh = ((8*np.pi)/cl**3) * f_photons[1:]**2 / (findCst(Ndens0)*np.exp((h*f_photons[1:])/(k*Te)) - 1 ) 
    #plotDensityFreq(xa, dxa, tobs, dt)
    #phoDensTh = ((8*np.pi)/cl**3) * f_photons**2 * (1 / (findCst(Ndens0)*np.exp((h*f_photons)/(k*Te)) - 1 )) *1.602e-16 / h 
    #plotDensity(xa, dxa, tobs, dt)
    #plotOccupationRate(xa, dxa, tobs, dt)


def mainKompaneets():
    
    xa,xb,dxa,dxb,tobs,dt,dto, M, e_photons, f_photons = meshGeneration()

    solveKompaneets(xa, xb, dxa, dxb, tobs, dt, dto, M, e_photons, f_photons)
    

mainKompaneets()
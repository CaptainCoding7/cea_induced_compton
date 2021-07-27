# -*- coding: utf-8 -*-
"""
Created on Thu Jul  8 15:21:30 2021

@author: pa267340
"""
from numpy import pi, exp, logspace, log10
from math import sqrt
from constants import k_cgs,k,h_cgs,h,cl,me,qe_cgs,Rcgs,R,sT
from matplotlib.pyplot import plot
from plotter import setFigureParameters
from changcooper import changCooper2
from initDistrib import init_gaussian

def j_Brem(T, pT, Rcgs, nu):
    """
    Applied the formula with 
    => n_i x n_e = Ne (neutral middle)
    => the ion charge Z = 1 (one proton on each ion)
    => the Gaunt factor = 1
    => the division by 4pi to have the emission per str
    The emission returned is in erg/s/Hz/cm^3/str (cgs)
    """
    
    qe = qe_cgs
    c = cl*1e2
    m = me * 1e3
    k = k_cgs
    h=h_cgs
    Ne = pT / Rcgs / sT / 1e4# (cm^-3)
    
    a = 2**5*pi*qe**6 / 3/m/c**3
    b = sqrt(2*pi/3/k/m)
    print(a*b)
     
    return a*b/sqrt(T)*Ne**2*exp(-h*nu/k/T)/4/pi



def meshGeneration(Te, Ne):
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
    tobs = [3e-3] # stabilité atteinte    
    
    tmax = tobs[-1]
    
    dt=1e-6
    dto = dt
    dto =  to * dt
    #print("dt;dto: ",dt,dto)
    N=int(tmax/dt)
    
    #### frequency mesh  ################################

    M=200
    numin=1e12
    numax=1e20
    #valeurs au bords des bins (M+1 valeurs)
    nu = logspace(log10(numin),log10(numax),M+1)
    
    return nu, tobs,tmax,dt,dto, M


def solveTransfertEq(nu,tobs,tmax,dt,dto,M, Te, kTe, Ne,pT):
    """
    solve the Transfert equation with the chang-cooper scheme
    call on of the plotting function
    """
    nua = nu[1:]
    j = j_Brem(Te, pT, Rcgs, nua)
    j_si = j * 1e-7 * 1e6
    print(j)
    print(j_si)
    
    a = h*nua/k/Te
    Bnu = ((2*h)/cl**2) * nua**3 *  (exp(-a) / (1 - exp(-a)))
    alpha = j / Bnu

    ################# TRANSFERT EQUATION  PARAMETERS  #################
    # we use the chang-cooper scheme without taking account of the A,B,C vectors
    
    # emission term
    Q = cl*j
    
    # absorption term
    T = 1/cl/alpha
    
    # photons injection at 0.5 keV, gaussian width of 0.03keV
    xinj = 1e14
    width = 1e2
    I0=init_gaussian(nu[1:], nu, xinj, width)
    epho = nu * h / 1.602e-16
    #I0 = init_planckian(epho[1:])
    
    # find the solution of the transfert equation with the cc scheme for each instant
    # specified in tobs
    N=int(tobs[-1]/dt)
    I = changCooper2(nu, N, I0,T, Q, dt)

    #t=np.logspace(np.log10(1e-14),np.log10(1e-10),100)
    t = tobs[-1]
    Ith = I0 * exp(-cl*alpha*t) + Bnu*(1-exp(-cl*alpha*t))
    #print(alpha)
    
    intensity = ((2*h)/cl**2) * (I*nu[1:]**3) 

    plot(nu[1:],Ith, '+', label = 'Ith')
    plot(nu[1:],I0, label = 'I0')
    plot(nu[1:],I, label = 'I_nu')
    plot(nu[1:],Bnu, label = 'B_nu')
    setFigureParameters(r"Intensité spécifique du Bremsstralhung pour un plasma chauffé à {:} keV".format(kTe)
                        ,r'$I_\nu (erg/s/Hz/cm^3/str)$',r'$\nu$',1e-10,1e5, nu[0],nu[-1])      
    
    return I


def testBrem():

    kTe = 10 #(keV)
    Te = kTe * 1.602e-16 / k #(K)
    pT = 2.36   
    Ne = pT / Rcgs / sT / 1e4# (cm^-3)
    #Ne = pT / R / sT # (m^-3)
    emin = 5e-6
    emax = 5e3
    e_pho = logspace(log10(emin),log10(emax),200)   
    nu = e_pho / h * 1.602e-16
    
    """
    j = j_Brem(Te, Ne, nu)    
    a = h*nu/k/Te
    Bnu = ((2*h)/cl**2) * nu**3 *  (exp(-a) / (1 - exp(-a)))
    alpha = j / Bnu
    print(j)
    plot(e_pho,j)
    setFigureParameters("Bremsstrahlung emission", r"$j_\nu$ (erg/s/Hz/m^3/str)", "Energy (keV)", 1e-20, 1e-8, 1e-5, 1e4)
    """
    
    nu,tobs,tmax,dt,dto, M = meshGeneration(Te, Ne)
    #alpha_nu_theta(B, p, nuL, Ub, theta, K2)
    
    solveTransfertEq(nu, tobs, tmax, dt, dto, M, Te, kTe, Ne,pT)    

    
#testBrem()    
    
    
    
    
    
    
    
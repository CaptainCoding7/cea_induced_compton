# -*- coding: utf-8 -*-
"""
Created on Tue Jun 22 11:47:30 2021

@author: pa267340
"""

import numpy as np
from scipy.stats import norm
from scipy.signal import unit_impulse
from constants import Tp, kTp

def init_gaussian(xb,xinj, width):
    """
    initialization of the photons distribution with a gaussian profile
    """

    # on identifie le bin correspondant a xinj
    iinj = np.argwhere(xb<xinj)[-1]
    # on recentre xinj sur le bin en question
    xinj = xb[iinj]
    u0 = norm.pdf(xb[1:] , loc = xinj , scale = width )
    
    return u0


def init_dirac(xa,xb,dxa,xinj):
    """
    initialization of the photons distribution with a dirac profile
    """
        
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

def init_planckian(e_photons):
    """
    initialization of the photons distribution with a planckian profile
    """
    print(kTp)
    
    #u0 = (1 /  (np.exp((h*f_photons)/(k*Tp)) - 1 ))
    # reshaping of the expression to avoid overflows
    u0 = np.exp(-e_photons/kTp) / (1 - (np.exp(-e_photons/kTp) ))
    return u0


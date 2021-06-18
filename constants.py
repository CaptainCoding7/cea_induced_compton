# -*- coding: utf-8 -*-
"""
Created on Thu Jun 17 14:35:07 2021

@author: pa267340
"""

import numpy as np
    
#################  CONSTANTS  #################

# Boltzman constant
k = 1.381e-23
# Planck constant
h = 6.626e-34
# Temperature of the electron field (K)
# kTe in keV 
kTe = 100.
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
    
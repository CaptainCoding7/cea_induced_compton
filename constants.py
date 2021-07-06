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
cl = 3e8
# Thompson scattering cross section
sT = 6.652e-29
# Electron mass (in kg)
me = 9.109e-31
# Electron charge (in C)
qe = 1.602e-19
#permeabilité magnetique du vide
mu0 = np.pi*4e-7


#################  PARAMETERS  #################

# Schwarzschild radius (Rs = 5e7 cm )
Rcgs = 5e7
R = Rcgs * 1e-2
# Thompson optical depth
pT = 2.36
# Electron number density (in m^-3)
Ne = pT / R / sT
# Temperature of the electron field
# kTe in keV 
kTe = 80
# Te in K
Te = kTe * 1.602e-16 / k
# Magnetic field (in T)
Bcgs=1e5
B=Bcgs*1e-4

# Photon temperature (used to initialize the planckian photon distribution)
kTp = 2
Tp = kTp * 1.602e-16 / k



    
    


    
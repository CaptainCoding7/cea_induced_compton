# -*- coding: utf-8 -*-
"""
Created on Wed Jun 23 16:32:34 2021

@author: Paul
"""

import numpy as np
from scipy.integrate import quad
from scipy.optimize import fsolve
from scipy.special import zeta

from constants import k,h,cl

    
def BEDistrib(x,c):
    """
    return the Bose-Einstein solution from x and c
    """
    return x**2 / (c*np.exp(x)-1)
    

def F(c, Ndens0, Te):
    """
    the equation used to find the constant
    """
    eq = quad(BEDistrib, 0, 100, args=(c,)) - ((cl*h)**3*Ndens0) / (8*np.pi*(k*Te)**3)
    return eq[0]


def findCst(Ndens0):
    """
    Function used to find the cst of the theoritical Bose-Einstein solution
    """    
    root = fsolve(F,2, args=(Ndens0))
    print("cst = ",root[0])
    return root[0]

####################################################

def findNmax(Te):
    """
    the theoritical maximal density below which there is no BE condensate
    """
    Nmax = 8*np.pi*(((k*Te)/(cl*h))**3)*2*zeta(3)
    return Nmax
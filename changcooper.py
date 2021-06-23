# -*- coding: utf-8 -*-
"""
Created on Tue Jun 22 10:02:49 2021

@author: pa267340

Functions to perform the chang-cooper scheme
"""

import numpy as np


def tridag(a, b, c, r):
    """
    Solve linear system with tridiagonal coefficient matrix.
    a is the lower band, b is the diagonal, c is the upper band, and r
    is the right hand side. The solution is returned in u.
    [b1 c1 0 ...      ] [u1]  [r1] 
    [a1 b2 c2 0 ...    ] [ :]  [ :] 
    [ 0 a2 b3 c3 0 ...   ] [ ] = [ ] 
    [           ] [ ]  [ ] 
    [ ... 0 an-2 bn-1 cn-1] [ :]  [ :] 
    [    ... 0 an-1 bn ] [un]  [rn] 
    """ 
    
    n = len(b)
    tmp = np.zeros(n-1) # necessary temporary array 
    if b[0] == 0: 
      raise RuntimeError('System is effectively order N-1') 
    beta = b[0]
    u=np.zeros_like(r)
    u[0] = r[0] / beta 
      
    for i in range(1, n): # Decompose and forward substitution 
      tmp[i-1] = c[i-1] / beta 
      beta = b[i] - a[i-1] * tmp[i-1] 
      if beta == 0: 
        raise RuntimeError('Method failure')
      u[i] = (r[i] - a[i-1] * u[i-1]) / beta 
      
    for i in range(n-1, 0, -1): # Backward substitution 
      u[i-1] -= tmp[i-1] * u[i] 
      
    return u


def changCooper(A, C, xa, dxa, dxb, N, u, pech, Q, dT):
    """
    solve the equation with the chang-cooper scheme
    """
    
    #print(dT)
    
    #dT=dt

    for n in range(N-1):
        
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
    
        # diagonale (taille M)
        b = np.zeros_like(xa)
        b[0]    = 1/dT + 1/A[   0]/dxa[   0] * (                         0 + C[0 ]/dxb[ 0]*Wm[ 0] ) + pech[0]
        b[1:-1] = 1/dT + 1/A[1:-1]/dxa[1:-1] * ( C[:-1]/dxb[:-1] * Wp[:-1] + C[1:]/dxb[1:]*Wm[1:] ) + pech[1:-1]
        b[-1]   = 1/dT + 1/A[  -1]/dxa[  -1] * ( C[ -1]/dxb[ -1] * Wp[ -1] +                    0 ) + pech[-1]
        
        r = Q + u/dT
                
        u = tridag(-a,b,-c,r)
        
        #print("---- u = ",u)
            
    return u


def changCooper2(xa, N, u, pech, Q, dT):
    """
    solve the equation with a simplified chang-cooper scheme
    (A, B, C parameters are null in the original equation)
    """

    for n in range(N-1):
              
        # sous-diagonale (taille M-1)
        a = np.zeros(len(xa)-1)
    
        # sur-diagonale (taille M-1)
        c = np.zeros(len(xa)-1)
    
        # diagonale (taille M)
        b = np.zeros_like(xa)
        b[0]    = 1/dT + pech[0]
        b[1:-1] = 1/dT + pech[1:-1]
        b[-1]   = 1/dT + pech[-1]
        
        r = Q + u/dT
                
        u = tridag(-a,b,-c,r)
        
        #print("---- u = ",u)
            
    return u
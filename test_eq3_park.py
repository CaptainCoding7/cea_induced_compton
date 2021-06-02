# Script to find a solution for the equation (24) of the Park paper with the Chang-Cooper scheme

import numpy as np
from scipy.stats import norm
from scipy.signal import unit_impulse
from scipy.stats import lognorm
from math import pow
from math import exp
from math import sinh
import tridiagonal_matrix_solver as tms
from plotter import linplot
from plotter import logplot
from scipy.special import iv as Inu
import matplotlib.pyplot as plt

"""
x: array of the energy x=(h.nu)/(k.T) with nu the frequency
u: solutions of the equation for different instants

"""


# Solution Analytique (Park-Petrosian 1995)
def uth(x,t,xinj):
    q = 3.0
    theta = 1.0
    s = 0.0
    a = 1.0
    alpha = (2-q-s)/2.0
    nu = np.abs((q-1+a)/2/alpha)
    G = 1/np.abs(alpha)/2/t * x**((1-q+a)/2)*xinj**((1-q-a)/2)
    G *= Inu(nu,(x*xinj)**alpha/2/t/alpha**2)
    G *= np.exp(-(x**(2*alpha)+xinj**(2*alpha))/4/alpha**2/t-theta*t)
    return G


def compute_am(dt, Am, dxa, dxa_inf1_2, Cm_inf1_2, Wm_inf1_2):

    return ((dt/(Am*dxa)) * (Cm_inf1_2/dxa_inf1_2) * Wm_inf1_2)

def compute_bm(dt, Am, dxa, dxa_inf1_2, dxa_sup1_2, Cm_inf1_2, Cm_sup1_2, Wp_inf1_2, Wm_sup1_2, Tm):

    return ( 1 + dt/(Am*dxa) * ((Cm_inf1_2/dxa_inf1_2) * Wp_inf1_2 + (Cm_sup1_2/dxa_sup1_2) * Wm_sup1_2) + dt/Tm)

def compute_cm(dt, Am, dxa, dxa_sup1_2, Cm_sup1_2, Wp_sup1_2):
    
    return ((dt/(Am*dxa)) * (Cm_sup1_2/dxa_sup1_2) * Wp_sup1_2)
    
    
def compute_rm(dt, Qm, u_n_m):
    
    return dt*Qm + u_n_m


def synch(M,x):

    return x*0


def init_sol(x,N):
    
    temperature=70
    xinj=0.1

    # initialization of the solution vector
    u=np.empty([N,M])
    #print(x)
    
    # we search the index of the x array with the closest value from xinj
    # calculate the difference array
    difference_array = np.absolute(x-xinj)
    # find the index of minimum element from the array
    index = difference_array.argmin()
    print(x[index], index)
    
    u[0] = temperature*unit_impulse(M, index)
    #u[0] = lognorm.pdf(x , loc = temperature , s = 0.01 )
    
    return u

def init_cond(xinj):
    # on identifie le bin correspondant a xinj
    iinj = np.argwhere(xb<xinj)[-1]
    # on recentre xinj sur le bin en question
    xinj = xa[iinj]
    # on initialise avec des zeros partout sauf dans le bin en question
    # ou on met 1/dx de maniere avec avoir une integrale=1
    u0 = np.zeros_like(xa)
    u0[iinj] = 1.0/dxa[iinj]
    return u0

# solve the equation with the chang-cooper scheme
def changCooper(u0, N):

    u=u0

    """
  
    # other vectors for the next calculus
    #w = (B/C)*dxa
    W = np.empty(M)
    Wm = np.empty(M)
    Wp = np.empty(M)
    # as we have sinh and exp operations, we can't do operations between sequances,
    # a for loop is necessary
    for m in range(len(w)):
        
        # W[m] takes two different values according to the value of w[m]
        abs_w = abs(w[m])
        if abs_w < 0.1:
            W[m] = 1 / ( 1 + pow(w[m],2)/24 + pow(w[m],4)/1920 )
        else:
            W[m] = ( abs_w * exp(-abs_w/2) ) / ( 1 - exp(-abs_w))
        
        #W[m] = (w[m]/2)/sinh(w[m]/2)
        Wm[m] = W[m] * exp(-w[m]/2)
        Wp[m] = W[m] * exp(w[m]/2)
    
    # the vectors used to solve the linear system
    a=np.empty(M)
    b=np.empty(M)
    c=np.empty(M)
    r=np.empty(M)
            
# computation of the a, b, c and r vectors

    
    for m in range(0, M-1):
    
        #print("x:",x[m])
        
        # now we evaluate the +1/2 and -1/2 terms
        if m!=M-2:
            
            #+1/2 terms
            dxa_sup1_2=(xa[m]+xa[m+1])/2
            Cm_sup1_2=(C[m]+C[m+1])/2
            Wp_sup1_2 = (Wp[m]+Wp[m+1])/2
            Wm_sup1_2 = (Wm[m]+Wm[m+1])/2
        else: 
            # if m=M-1, all +1/2 terms are = 0
            # the value of dxa_sup_1_2 is not important as it will only
            # appear as denominator in fractions with numerator = 0
            dxa_sup1_2 = 1
            Cm_sup1_2 = Wp_sup1_2 = Wm_sup1_2 = 0

        if m!=0:
            
            # -1/2 terms
            dxa_inf1_2=(xa[m]+xa[m-1])/2
            Cm_inf1_2=(C[m]+C[m-1])/2
            Wm_inf1_2 = (Wm[m]+Wm[m-1])/2
            Wp_inf1_2 = (Wp[m]+Wp[m-1])/2
        
        else:
            # if m==0, all -1/2 terms are = 0
            dxa_inf1_2 = 1
            Cm_inf1_2 = Wm_inf1_2 = Wp_inf1_2 = 0

        
        #then we evaluate the values of the vectors for each m point
            
        a[m] = compute_am(dt, A[m], dxa[m], dxa_inf1_2, Cm_inf1_2, Wm_inf1_2)
        
        b[m] = compute_bm(dt, A[m], dxa[m], dxa_inf1_2, dxa_sup1_2, Cm_inf1_2, Cm_sup1_2, Wp_inf1_2, Wm_sup1_2, T[m])
    
        c[m] = compute_am(dt, A[m], dxa[m], dxa_sup1_2, Cm_sup1_2, Wp_sup1_2)
    
    """
        
    w = (B/C)*dxb
    lwl = np.abs(w)
    W = lwl * np.exp(-lwl)/(1.0-np.exp(-lwl))
    Wp = W*np.exp(+w/2)
    Wm = W*np.exp(-w/2)
    
    
    # sous-diagonale (taille M-1)
    a = 1/A[1:]/dxa[1:] * C/dxb * Wm

    # sur-diagonale (taille M-1)
    c = 1/A[:-1]/dxa[:-1] * C/dxb * Wp

    # doagonale (taille M)
    b = np.zeros_like(xa)
    b[0]    = 1/dt + 1/A[   0]/dxa[   0] * (                         0 + C[0 ]/dxb[ 0]*Wm[ 0] ) + 1/T[0]
    b[1:-1] = 1/dt + 1/A[1:-1]/dxa[1:-1] * ( C[:-1]/dxb[:-1] * Wp[:-1] + C[1:]/dxb[1:]*Wm[1:] ) + 1/T[1:-1]
    b[-1]   = 1/dt + 1/A[  -1]/dxa[  -1] * ( C[ -1]/dxb[ -1] * Wp[ -1] +                    0 ) + 1/T[-1]
    


    # simulation for different instants
    for n in range(N):
        r = Q + u/dt
        u = tms.tridag(-a,b,-c,r)
            
    return u



####################################################
## Mesh generation

# time
#Nmax=30.01
t = [0.3,3.,30.]
dt=0.01
#N=int(Nmax/dt)

# energy
M=100
xmin = 1.e-3
xmax = 1.e3

#x=np.logspace(-3, 3, M, endpoint=True)


# Valeurs sur les bords des bins (M+1 values)
xb = np.logspace(np.log10(xmin),np.log10(xmax),M+1)
# Largeur des bins (M valeurs)
dxa = xb[1:] - xb[:-1]
# valeurs au centre des bins (M valeurs)
xa = 0.5 * ( xb[1:] + xb[:-1] )
# distance entre les milieux (M-1 valeurs)
dxb = xa[1:] - xa[:-1]

"""
# we evaluate Δx
dxa=np.empty(M)

for m in range(M):
    if m!=0 and m!=M-1:
        dxa[m] = (x[m+1]-x[m-1])/2
        #dxa[m] = x[m+1]-x[m]

    elif m==0:
        dxa[m] = x[m+1]-x[m]
        
    elif m==M-1:
        dxa[m] = x[m]-x[m-1] 
"""        

# initialization of the solution vector
#u=np.empty([N,M])
u0=init_cond(0.1)


# Initialization of the A, B, C, Q, w and W vectors
# the vectors containing the values obtained with the specific function 
# in the equation (34)

A = np.ones_like(xa)            # au milieu des bins (M)
B = -0.5*(xa[1:]**2+xa[:-1]**2) # sur les bords      (M-1)
C =  0.5*(xa[1:]**3+xa[:-1]**3) # sur les bords      (M-1)
Q = np.zeros_like(xa)           # au milieu          (M)
T = np.ones_like(xa)            # au milieu          (M)

"""
# sink term
#n*ro = n/T(x) => T(x) = 1/ro
ro=0.2
T=np.full(M,1)

A = np.ones(M)
B = -x*x
C = x*x*x
Q = synch(M,x)
"""
    
logplot(dt, xa, u0, 0)

# plotting for 3 instants
t = [0.3,3.,30.]
col=['r','g','b']
norm = [1,1,1.e11]
for i,tt in enumerate(t):
    n = int(tt/dt)
    u = changCooper(u0,n)
    print(u)
    #logplot(dt, xa, norm[i]*u, n)
    #logplot(dt,xa, norm[i]*uth(xa,n*dt, 0.1),tt )
    plt.plot(xa,norm[i]*u,':',color=col[i],label='t={:.2f}'.format(tt))
    plt.plot(xa,norm[i]*uth(xa,n*dt,0.1),color=col[i])

# Reglages affichage
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.ylim(1.e-8,1.e2)
plt.xlim(1.e-4,1.e4)
plt.xlabel('x')
plt.ylabel('u(x,t)')
plt.show()


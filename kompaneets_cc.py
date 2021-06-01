# Script to find a solution for the Kompaneets equation with the Chang-Cooper scheme

import numpy as np
from scipy.stats import norm
from scipy.signal import unit_impulse
from scipy.stats import lognorm
import matplotlib.pyplot as plt
import seaborn as sb
from math import pow
from math import exp
from math import sinh
import tridiagonal_matrix_solver as tms

"""
x: array of the energy x=(h.nu)/(k.T) with nu the frequency
u: solutions of the Kompaneets equation for different instants

"""

def compute_am(dt, Am, dxm, dxm_inf1_2, Cm_inf1_2, Wless_inf1_2):

    return ((dt/(Am*dxm)) * (Cm_inf1_2/dxm_inf1_2) * Wless_inf1_2)



def compute_bm(dt, Am, dxm, dxm_inf1_2, dxm_sup1_2, Cm_inf1_2, Cm_sup1_2, Wplus_inf1_2, Wless_sup1_2, Tm):

    return ( 1 + dt/(Am*dxm) * ((Cm_inf1_2/dxm_sup1_2) * Wplus_inf1_2 + (Cm_sup1_2/dxm_sup1_2) * Wless_sup1_2))# + dt/Tm)

def compute_cm(dt, Am, dxm, dxm_sup1_2, Cm_sup1_2, Wplus_sup1_2):
    
    return ((dt/(Am*dxm)) * (Cm_sup1_2/dxm_sup1_2) * Wplus_sup1_2)
    
    
def compute_rm(dt, Qm, u_n_m):
    
    return dt*Qm + u_n_m


def synch(M,x):

    return x*0
    

def linplot(x, u, N):

    for n in range(0, N):
       plt.plot(x, u[n], label="t=0.{0}s".format(n))

    plt.plot(x, u[0], label="t=0.1s")

    leg = plt.legend(loc="upper right", bbox_to_anchor=[1, 1],
                     ncol=2, shadow=True, title="Legend", fancybox=True)
    leg.get_title().set_color("black")

    plt.xlim(1e-3, 3)
    plt.xlabel('Energy (kEv)')
    plt.ylabel('Particle Number Density')

def logplot(x, u, N):

    for n in range(0, N):
       plt.plot(x, u[n], label="t=0.{0}s".format(n))

    plt.plot(x, u[0], label="t=0.1s")

    leg = plt.legend(loc="upper right", bbox_to_anchor=[1, 1],
                     ncol=2, shadow=True, title="Legend", fancybox=True)
    leg.get_title().set_color("black")

    #plt.plot(np.log10(data), pdf)
    plt.xscale('log')
    plt.yscale("log")
    plt.ylim(1e-10,10)
    #plt.xlim(1e-1, 10)
    plt.xlabel('Energy (kEv)')
    plt.ylabel('Particle Number Density')

# solve the kompaneets equation with the chang-cooper scheme
def kompaneets():

    ####################################################
    ## Mesh generation and initialization of the parameters

    M=100
    N=10

    x=np.logspace(-3, 3, M, endpoint=True)
    #x=np.linspace(0.001, 1000, M, endpoint=True)
            
    dt=0.1
    temperature=1

    # initialization of the solution vector
    u=np.empty([N,M])
    u[0] = norm.pdf(x , loc = temperature , scale = 0.5 )
    #u[0] = 10*unit_impulse(M, 50)
    

    # sink term
    #n*ro = n/T(x) => T(x) = 1/ro
    ro=0.2
    T=np.full(M,1/ro)
    
    # Initialization of the A, B, C, Q, w and W vectors
    # the vectors containing the values obtained with the specific function 
    # in the Kompaneets equation
    A = x*x
    B = x*x*x*x
    C = x*x*x*x
    
    #synchrotron term
    Q = synch(M,x)
    
    # we evaluate Δx
    dxm=np.empty(M)
    
    for m in range(M):
        if m!=0 and m!=M-1:
            dxm[m] = (x[m+1]-x[m-1])/2
            #dxm[m] = x[m+1]-x[m]

        elif m==0:
            dxm[m] = x[m+1]-x[m]
            
        elif m==M-1:
            dxm[m] = x[m]-x[m-1] 
            
    w = (B/C)*dxm

    # other vectors for the next calculus
    W = np.empty(100)
    Wless = np.empty(100)
    Wplus = np.empty(100)
    # as we have sinh and exp operations, we can't do operations between sequances,
    # a for loop is necessary
    for m in range(M):
        
        # W[m] takes two different values according to the value of w[m]
        abs_w = abs(w[m])
        if abs_w < 0.1:
            W[m] = 1 / ( 1 + pow(w[m],2)/24 + pow(w[m],4)/1920 )
        else:
            W[m] = ( abs_w * exp(-abs_w/2) ) / ( 1 - exp(-abs_w))
        
        #W[m] = (w[m]/2)/sinh(w[m]/2)
        Wless[m] = W[m] * exp(-w[m]/2)
        Wplus[m] = W[m] * exp(w[m]/2)
        

    # simulation for different instants
    for n in range(N):
        
        # the vectors used to solve the linear system
        a=np.empty(100)
        b=np.empty(100)
        c=np.empty(100)
        r=np.empty(100)
                
    # computation of the a, b, c and r vectors
        for m in range(0, M):
        
            #print("x:",x[m])
            
            # now we evaluate the +1/2 and -1/2 terms
            if m!=M-1:
                
                #+1/2 terms
                dxm_sup1_2=(x[m]+x[m+1])/2
                Cm_sup1_2=(C[m]+C[m+1])/2
                Wplus_sup1_2 = (Wplus[m]+Wplus[m+1])/2
                Wless_sup1_2 = (Wless[m]+Wless[m+1])/2
            else: 
                # if m=M-1, all +1/2 terms are = 0
                # the value of dxm_sup_1_2 is not important as it will only
                # appear as denominator in fractions with numerator = 0
                dxm_sup1_2 = 1
                Cm_sup1_2 = Wplus_sup1_2 = Wless_sup1_2 = 0

            if m!=0:
                
                # -1/2 terms
                dxm_inf1_2=(x[m]+x[m-1])/2
                Cm_inf1_2=(C[m]+C[m-1])/2
                Wless_inf1_2 = (Wless[m]+Wless[m-1])/2
                Wplus_inf1_2 = (Wplus[m]+Wplus[m-1])/2
            
            else:
                # if m==0, all -1/2 terms are = 0
                dxm_inf1_2 = 1
                Cm_inf1_2 = Wless_inf1_2 = Wplus_inf1_2 = 0

                
            #then we evaluate the values of the vectors for each m point
                
            a[m] = compute_am(dt, A[m], dxm[m], dxm_inf1_2, Cm_inf1_2, Wless_inf1_2)
            
            b[m] = compute_bm(dt, A[m], dxm[m], dxm_inf1_2, dxm_sup1_2, Cm_inf1_2, Cm_sup1_2, Wplus_inf1_2, Wless_sup1_2, T[m])

            c[m] = compute_am(dt, A[m], dxm[m], dxm_sup1_2, Cm_sup1_2, Wplus_sup1_2)

            # we evaluate r with the values of u found at the previous isntant
            # now = n+1 // past = n
            r[m] = compute_rm(dt, Q[m], u[n][m])

        if n<N-1:
            u[n+1] = tms.tridag(a,b,c,r)
            #print(len(u[n+1]))
            #print(len(x))
    logplot(x,u,N)
    

kompaneets()
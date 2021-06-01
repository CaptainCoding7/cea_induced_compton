# Script to find a solution for the Kompaneets equation with the Chang-Cooper scheme

import numpy as np
from scipy.stats import norm
from scipy.stats import lognorm
import matplotlib.pyplot as plt
import seaborn as sb
from math import pow
from math import exp
from math import sinh
import tridiagonal_matrix_solver as tms

"""
nx: size of the frequencies array
x: array of the energy x=(h.nu)/(k.T) with nu the frequency


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

    return 0
    

def plot(x, u, N):

    for n in range(0, N):
        plt.plot(x, u[n], label="t=0.{0}s".format(n))
    leg = plt.legend(loc="upper right", bbox_to_anchor=[1, 1],
                     ncol=2, shadow=True, title="Legend", fancybox=True)
    leg.get_title().set_color("black")

    #plt.plot(np.log10(data), pdf)
    plt.xscale('log')
    #plt.yscale("log")
    plt.xlabel('Energy (kEv)')
    plt.ylabel('Particle Number Density')

# solve the kompaneets equation with the chang-cooper scheme
def kompaneets():

    # Mesh generation

    M=100
    N=10

    x=np.logspace(-3, 3, M, endpoint=True)
    #x=np.linspace(0.001, 1000, M, endpoint=True)

            
    dt=0.01
    temperature=1

    # initialization of the solution vector
    u=[ [ 0 for m in range(M) ] for n in range(N) ]
    u[0] = norm.pdf(x , loc = temperature , scale = 0.5 )

    
    # the vectors containing the values obtained with the specific function 
    # in the Kompaneets equation
    A = [ 0 for m in range(M)]
    B = [ 0 for m in range(M)]
    C = [ 0 for m in range(M)]
    
    # synchrotron
    Q = [ 0 for m in range(M)]

    # diffusive term
    #n*ro = n/T(x) => T(x) = 1/ro
    ro=0.2
    T=[ 1/ro for m in range(M)]
    #T=[ 1 for m in range(M)]

    # other vectors for the next calculus
    w = [ 0 for m in range(M)]
    W = [ 0 for m in range(M)]
    Wless = [ 0 for m in range(M)]
    Wplus = [ 0 for m in range(M)]


    
    # Initialization of the A, B, C, Q, w and W vectors
    for m in range(M-1):
        A[m] = x[m]*x[m]
        B[m] = pow(x[m],4)
        C[m] = pow(x[m],4)
        
        Q[m] = synch(M,x)
        
        dxm = x[m+1]-x[m]        
        w[m] = (B[m]/C[m])*dxm
        W[m] = (w[m]/2)/sinh(w[m]/2)
        Wless[m] = W[m] * exp(-w[m]/2)
        Wplus[m] = W[m] * exp(w[m]/2)

        

    # simulation for different instants
    for n in range(N):
        
        # the vectors used to solve the linear system
        a=[]
        b=[]
        c=[]
        r=[]
        
    # computation of the a, b, c and r vectors
        for m in range(1, M-1):
            
            
            print("x:",x[m])
            # First we evaluate all the middle points
            dxm=x[m]-x[m-1]
            dxm_inf1_2=(x[m]+x[m-1])/2
            dxm_sup1_2=(x[m]+x[m+1])/2

            Cm_inf1_2=(C[m]+C[m-1])/2
            Cm_sup1_2=(C[m]+C[m+1])/2

            Wless_inf1_2 = (Wless[m]+Wless[m-1])/2
            Wplus_sup1_2 = (Wplus[m]+Wplus[m+1])/2
            Wless_sup1_2 = (Wless[m]+Wless[m+1])/2
            Wplus_inf1_2 = (Wplus[m]+Wplus[m-1])/2

            #then we evaluate the values of the vectors for each m point
                
            a.append(compute_am(dt, A[m], dxm, dxm_inf1_2, Cm_inf1_2, Wless_inf1_2))
            
            b.append(compute_bm(dt, A[m], dxm, dxm_inf1_2, dxm_sup1_2, Cm_inf1_2, Cm_sup1_2, Wplus_inf1_2, Wless_sup1_2, T[m]))

            c.append(compute_am(dt, A[m], dxm, dxm_sup1_2, Cm_sup1_2, Wplus_sup1_2))

            # we evaluate r with the values of u found at the previous isntant
            # now = n+1 // past = n
            r.append(compute_rm(dt, Q[m], u[n][m] ))

    # as we can't do the same calculus for the values at the boundaries (for example a[-1] doesn't exists)
    # we take initialize the boundary values with the close ones
        a.insert(0,-a[0])
        b.insert(0,-b[0])
        c.insert(0,-c[0])        
        r.insert(0,-r[0])
                
        a.append(a[M-2]+ (a[M-2]-a[M-3]))
        b.append(b[M-2]+ (b[M-2]-b[M-3]))
        c.append(c[M-2]+ (c[M-2]-c[M-3]))
        r.append(r[M-2]+ (r[M-2]-r[M-3]))

        if n<N-1:
            u[n+1] = tms.tridag(a,b,c,r)
            #print(len(u[n+1]))
            #print(len(x))
    plot(x,u,N)
    
    
    
def testDistrib():

    M=100
    
    N=10
    x = np.linspace(0, 10, M, endpoint=True)
    #x = np.logspace(0, 10, M, endpoint=True)
    #np.arange(1,10,0.01)
    #t=np.logspace(0, 10, nbT, endpoint=True)
    
    #dx=[ ((x[m+1] + x[m-1])/2) for m in range(1,M-1)]
    dt=0.01
    
    # initialization of the solution vector
    u=[ [ 0 for m in range(M) ] for n in range(N) ]
    #u[0] = lognorm.pdf(x , s=1, loc = 5 , scale = 1 )
    u[0] = norm.pdf(x, loc = 5 , scale = 1 )
    
    u[0][M-1]=0

    plot(x,u[0],'black')
    

kompaneets()
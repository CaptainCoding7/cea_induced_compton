# -*- coding: utf-8 -*-
"""
Created on Tue Jun  1 16:26:31 2021

@author: Paul
"""
from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np
from constants import cl, h, k , Te


def setFigureParameters(title, ylabel, xlabel, ymin, ymax,  xmin, xmax):
    
    # Reglages affichage
    plt.xscale('log')
    plt.yscale('log')
    leg = plt.legend(loc="upper left",prop={'size': 7}, bbox_to_anchor=[0, 1],
                     ncol=1, shadow=True, title="Legend", fancybox=True)
    leg.get_title().set_color("black")
    plt.title(title)
    plt.ylim(ymin,ymax)
    plt.xlim(xmin,xmax)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.show()
    
def noLogsetFigureParameters(title, ylabel, xlabel, ymin, ymax,  xmin, xmax):
    
    # Reglages affichage
    plt.yscale('log')
    leg = plt.legend(loc="upper right",prop={'size': 7}, bbox_to_anchor=[1, 1],
                     ncol=1, shadow=True, title="Legend", fancybox=True)
    leg.get_title().set_color("black")
    plt.ylim(ymin,ymax)
    plt.xlim(xmin,xmax)
     #   plt.xticks([1, 2, 3, 4, 5])
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.show()


def plotOccupationRate(uobs, tobs, xa):
    col = [ cm.jet(x) for x in np.linspace(0, 0.3 , len(tobs))]

    print("Plotting occupation rate...")
    tmax= tobs[len(tobs)-1]

    # plot the occupation rate
    for i,tt in enumerate(tobs):
        print("Tracé ",i+1)
        plt.plot(xa,uobs[i],color=col[i],label='t={:.2f}s '.format(tt))    
    
    setFigureParameters('','u(x,t)','x')
    
    
def plotDensityFreq(phoDensTh, uobs, tobs, f_photons, e_photons):
    """
    plot the evolving photons density as a function of the frequency
    """
    col = [ cm.jet(x) for x in np.linspace(0, 0.3 , len(tobs))]
    print("Plotting spectral density...")
                  
    #plt.plot(f_photons, phoDensTh, "+",color='red',label='theoritical solution')

    for i,tt in enumerate(tobs):        
        # plot the photon density rate
        phoDens = ((8*np.pi)/cl**3) * (uobs[i]*f_photons**2)

        # multiply by Enorm as we represent the energy in abscisse
        plt.plot(f_photons, phoDens, color = col[i],label='t={:.2f}s '.format(tt))
    
    setFigureParameters('','Spectral Density $ (nbPhotons . m^{-3} . Hz^{-1})$','Frequency (Hz)',1e2,1e13,1e13,1e21)
    
    
def plotDensity(phoDensTh, uobs, tobs, f_photons, e_photons):
    """
    plot the evolving photons density as a function of the energy (keV)
    """
    col = [ cm.jet(x) for x in np.linspace(0, 0.3 , len(tobs))]
    print("Plotting density...")
                 
    plt.plot(e_photons, phoDensTh, "+",color='red',label='theoritical solution')

    for i,tt in enumerate(tobs):      
        # plot the photon density rate
        """
        two ways to convert the density from nbPho/m3/Hz into nbPho/m3/kev
        - 1kev = 2.4e17 Hz
        - E=h.nu ==> divide by h (and  multiply by 1.602e-16)
        """
        #phoDens = ((8*np.pi)/cl**3) * (uobs[i]*(f_photons**2)) * 2.4e17
        phoDens = ((8*np.pi)/cl**3) * (uobs[i]*(f_photons**2)*1.602e-16) / h 
        plt.plot(e_photons, phoDens, color = col[i],label='t={:.1E}s'.format(tt))
    
    setFigureParameters('','Photon Number Density $ (nbPhotons . m^{-3} . keV^{-1})$','Energy (keV)',1e20,1e30,1e-4,1e3)


def plotIntensity(iTh, uobs, tobs, e_photons, f_photons,
                  titre, ylab,xlab,ymin,sup,xmin,xsup):
    """
    plot the sepctral radiance (intensité specifique) as a function of the energy
    """
    col = [ cm.jet(x) for x in np.linspace(0, 0.3 , len(tobs))]
    
    # ANALYTICAL SOLUTION
    #plt.plot(e_photons, iTh, "+", color='red',label='th solution (no Synchrotron)')
    
    # plot the intensity of the initial photon field
    intensity = ((2*h)/cl**2) * (uobs[0]*f_photons**3)
    plt.plot(e_photons, intensity, color = 'black',label='t=0s')

    for i,tt in enumerate(tobs):
        # plot the intensity of the photon field
        intensity = ((2*h)/cl**2) * (uobs[i+1]*f_photons**3) 
        plt.plot(e_photons, intensity, color = col[i],label='t={:.1E}s'.format(tt))    

    setFigureParameters(titre, ylab,xlab,ymin,sup,xmin,xsup)


def plotAll(Bnu, nuLs, nuLc, nuL, tobs, e_pho, nu,
                  titre, ylab,xlab,ymin,sup,xmin,xsup):
    """
    plot the sepctral radiance (intensité specifique) as a function of the energy
    """
    col = [ cm.jet(x) for x in np.linspace(0, 0.3 , len(tobs))]
    
    # Synchrotron only
    plt.plot(e_pho, nuLs, color='red',label='Synchrotron radiation')

    # compton only
    #plt.plot(e_pho, nuLc, color = 'orange',label='Compton scattering')   

    # Blackbody emission
    #plt.plot(e_pho,Bnu, color='black',label='Blackbody emission')
    
    # plot the intensity of the initial photon field
    plt.plot(e_pho, nuL[0], color = 'purple',label='t=0s')

    for i,tt in enumerate(tobs):
        # plot the intensity of the photon field
        plt.plot(e_pho, nuL[i+1], color = col[i],label='t={:.1E}s'.format(tt))    
        
    setFigureParameters(titre, ylab,xlab,ymin,sup,xmin,xsup)







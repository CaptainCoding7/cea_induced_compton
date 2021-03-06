# -*- coding: utf-8 -*-
"""
Created on Thu Jul  1 15:02:31 2021

@author: pa267340
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from plotter import setFigureParameters
from plotter import noLogsetFigureParameters


def plotBELM(nuLB, nuLs, synchB, nuLc, compB, nuLeq, spectrumB, tobs, e_pho, eB,
                  titre, ylab,xlab,ymin,sup,xmin,xsup):
    """
    plot the sepctral radiance (intensité specifique) as a function of the energy
    """
    col = [ cm.jet(x) for x in np.linspace(0, 0.3 , len(tobs))]
    
    # Synchrotron only
    plt.plot(e_pho, nuLs, color='blue',label='Synchrotron')

    """
    nuLcNeg = np.empty_like(nuLc)
    for i,u in enumerate(nuLc):
        if u<0:
            nuLcNeg[i]=-u
            
    # Compton only
    plt.plot(e_pho, nuLcNeg, '--', color='orange',label='Negative Compton')
"""
    plt.plot(e_pho, nuLc, color='orange',label='Compton')

    # BELM Synchrotron 
    plt.plot(eB, synchB, color='purple',label='BELM Synchrotron ')

    # BELM compton
    plt.plot(eB, compB, color = 'red',label='BELM Compton')   

    # BELM Spectra
    plt.plot(eB, spectrumB, color = 'darkblue',label='BELM Spectrum')   

    # Blackbody emission
    plt.plot(e_pho,nuLB, color='black',label='Blackbody emission')
    
    # plot the intensity of the initial photon field
    #plt.plot(e_pho, nuL[0], color = 'black',label='t=0s')

    plt.plot(e_pho, nuLeq, color = 'green',label='Equilibrium spectrum')    
        
    setFigureParameters(titre, ylab,xlab,ymin,sup,xmin,xsup)


def plotBELMnoCompton(nuLB, nuLs, synchB, nuLeq, spectrumB, tobs, e_pho, eB,
                  titre, ylab,xlab,ymin,sup,xmin,xsup):
    """
    plot the sepctral radiance (intensité specifique) as a function of the energy
    """
    col = [ cm.jet(x) for x in np.linspace(0, 0.3 , len(tobs))]
    
    # Synchrotron only
    plt.plot(e_pho, nuLs, color='blue',label='Synchrotron')

    # BELM Synchrotron 
    plt.plot(eB, synchB, color='purple',label='BELM Synchrotron ')

    # BELM Spectra
    plt.plot(eB, spectrumB, color = 'darkblue',label='BELM Spectrum')   

    # Blackbody emission
    plt.plot(e_pho,nuLB, color='black',label='Blackbody emission')
    
    # plot the intensity of the initial photon field
    #plt.plot(e_pho, nuL[0], color = 'black',label='t=0s')

    # equilibrium spectrum
    plt.plot(e_pho, nuLeq, color = 'green',label='Equilibrium spectrum')    

    setFigureParameters(titre, ylab,xlab,ymin,sup,xmin,xsup)
    
    
def compareBELM(nuLB,nuLs,nuLc,nuL,tobs,e_pho,nu,
                titre, ylab,xlab,ymin,sup,xmin,xsup):

    """
    compare the results entered in arguments with the results get from the
    BELM model
    """
    
    firstline=41
    lastline=552
    nblines=512
    datafile = open("spectra.txt")
    lines_to_read = np.linspace(firstline,lastline,nblines)
    lines = np.empty((nblines,9))
    eB = np.empty(nblines)
    spectrumB = np.empty(nblines)
    synchB = np.empty(nblines)
    compB = np.empty(nblines)


    for i, line in enumerate(datafile):
    
        if i in lines_to_read:
            lines[i-firstline] = line.split()
    

#### We store the BELM spectra's data in arrays named as xB

    for i in range(len(lines)):
        eB[i] = lines[i][2]
        spectrumB[i] = lines[i][5]
        compB[i] = lines[i][6]
        synchB[i] = lines[i][7]
        
    
    plotBELM(nuLB,nuLs, synchB, nuLc, compB, nuL, spectrumB, tobs, e_pho, eB, 
             titre, ylab, xlab, ymin, sup, xmin, xsup)
    
def compareBELMnoCompton(nuLB,nuLs,nuLc,nuL,tobs,e_pho,nu,
                titre, ylab,xlab,ymin,sup,xmin,xsup):

    """
    compare the results entered in arguments with the results get from the
    BELM model
    the Compton scattering is not used in the simulation
    """
    firstline=41
    lastline=552
    nblines=512
    datafile = open("synch_escape.txt")
    lines_to_read = np.linspace(firstline,lastline,nblines)
    lines = np.empty((nblines,7))
    eB = np.empty(nblines)
    spectrumB = np.empty(nblines)
    synchB = np.empty(nblines)


    for i, line in enumerate(datafile):
    
        if i in lines_to_read:
            lines[i-firstline] = line.split()
    

#### We store the BELM spectra's data in arrays named as xB

    for i in range(len(lines)):
        eB[i] = lines[i][2]
        spectrumB[i] = lines[i][4]
        synchB[i] = lines[i][5]
        
    
    plotBELMnoCompton(nuLB,nuLs, synchB, nuL, spectrumB, tobs, e_pho, eB, 
             titre, ylab, xlab, ymin, sup, xmin, xsup)
    
    
        
def compareBELMIC(nuLIC,nuLnoIC,e_pho, titre, ylab,xlab,ymin,ymax,xmin,xmax):

    """
    compare the results entered in arguments with the results get from the
    BELM model
    plot both the cases where the induced compton is take into account and where it's not
    """
    
    firstline=41
    lastline=4136
    nblines=4096
    fileIC = open("Spectres/SSC_cygX-1_with_induced.txt")
    fileNoIC = open("Spectres/SSC_cygX-1_without_induced.txt")
    """
    fileIC = open("Spectres/BSC_GX339-4_avecinduit.txt")
    fileNoIC = open("Spectres/BSC_GX339-4_sansinduit.txt")
    """
    
    lines_to_read = np.linspace(firstline,lastline,nblines)
    linesIC = np.empty((nblines,8))
    linesNoIC = np.empty((nblines,8))
    eB = np.empty(nblines)
    spectB_IC = np.empty(nblines)
    spectBnoIC = np.empty(nblines)

    for i, line in enumerate(fileIC):
        if i in lines_to_read:
            linesIC[i-firstline] = line.split()
            
    for i, line in enumerate(fileNoIC):
        if i in lines_to_read:
            linesNoIC[i-firstline] = line.split()
        

#### We store the BELM spectra's data in arrays named as xB

    for i in range(nblines):
        eB[i] = linesIC[i][2]
        spectB_IC[i] = linesIC[i][4]
        spectBnoIC[i] = linesNoIC[i][4]
        
        

    # BELM Spectra
    plt.plot(eB, spectB_IC, color = 'darkblue',label='BELM Spectrum')   

    # BELM Spectra without IC
    plt.plot(eB, spectBnoIC, color = 'lightblue',label='BELM Spectrum without IC')   
    
    # equilibrium spectrum 
    plt.plot(e_pho, nuLIC, color = 'green',label='Equilibrium spectrum')    

    # equilibrium spectrum without IC
    plt.plot(e_pho, nuLnoIC, color = 'lightgreen',label='Equilibrium spectrum without IC')    

    setFigureParameters(titre, ylab,xlab,ymin,ymax,xmin,xmax)        
        
    plt.plot(eB, spectBnoIC/spectB_IC)
    plt.xscale('log')
    plt.grid()
    noLogsetFigureParameters("Rapport des spectres sans Compton induit / avec Compton induit", 
                             r'$\nu.L_\nu / \nu.Lind_\nu $','Energie (keV)', 0, 2, eB[0], eB[-1])
      
    
    #plotBELM(nuLB,nuLs, np.zeros(nblines), nuLc, np.zeros(nblines), nuL, spectB_IC, spectBnoIC, tobs, e_pho, eB, 
    #         titre, ylab, xlab, ymin, sup, xmin, xsup)

"""

xmin = 1e-5
xmax = 1e-6
ymin = 1e25
ymax = 1e35

z= np.zeros(100)

compareBELMIC(z,z,z,'BELM Spectra comparison', r'$\nu.L_\nu$ $(erg.s^{-1})$','Energy (keV)',
              ymin,ymax,xmin,xmax)

"""
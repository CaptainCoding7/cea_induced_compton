# -*- coding: utf-8 -*-
"""
Created on Thu Jul  1 15:02:31 2021

@author: pa267340
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from plotter import setFigureParameters


def plotBELM(nuLB, nuLs, synchB, nuLc, compB, nuL, spectrumB, tobs, e_pho, eB,
                  titre, ylab,xlab,ymin,sup,xmin,xsup):
    """
    plot the sepctral radiance (intensité specifique) as a function of the energy
    """
    col = [ cm.jet(x) for x in np.linspace(0, 0.3 , len(tobs))]
    
    # Synchrotron only
    plt.plot(e_pho, nuLs, color='blue',label='Synchrotron')

    # Compton only
    plt.plot(e_pho, nuLc, color='orange',label='Compton')

    # BELM Synchrotron 
    plt.plot(eB, synchB, color='purple',label='BELM Synchrotron ')

    # BELM compton
    plt.plot(eB, compB, color = 'red',label=' BELM Compton')   

    # BELM Spectra
    plt.plot(eB, spectrumB, color = 'darkblue',label=' BELM Spectrum')   

    # Blackbody emission
    plt.plot(e_pho,nuLB, color='black',label='Blackbody emission')
    
    # plot the intensity of the initial photon field
    #plt.plot(e_pho, nuL[0], color = 'black',label='t=0s')

    for i,tt in enumerate(tobs):
        # plot the intensity of the photon field
        plt.plot(e_pho, nuL[i+1], color = 'green',label='t={:.1E}s'.format(tt))    
        
    setFigureParameters(titre, ylab,xlab,ymin,sup,xmin,xsup)


def plotBELMnoCompton(nuLB, nuLs, synchB, nuL, spectrumB, tobs, e_pho, eB,
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
    plt.plot(eB, spectrumB, color = 'darkblue',label=' BELM Spectrum')   

    # Blackbody emission
    plt.plot(e_pho,nuLB, color='black',label='Blackbody emission')
    
    # plot the intensity of the initial photon field
    #plt.plot(e_pho, nuL[0], color = 'black',label='t=0s')

    for i,tt in enumerate(tobs):
        # plot the intensity of the photon field
        plt.plot(e_pho, nuL[i+1], color = 'green',label='t={:.1E}s'.format(tt))    
        print()
    setFigureParameters(titre, ylab,xlab,ymin,sup,xmin,xsup)
    
    
def compareBELM(nuLB,nuLs,nuLc,nuL,tobs,e_pho,nu,
                titre, ylab,xlab,ymin,sup,xmin,xsup):

    """
    compare the results entered in arguments with the results get from the
    BELM model
    """
    
    datafile = open("spectra.txt")
    lines_to_read = np.linspace(41,296,256)
    lines = np.empty((256,9))
    eB = np.empty(256)
    spectrumB = np.empty(256)
    compB = np.empty(256)
    synchB = np.empty(256)


    for i, line in enumerate(datafile):
    
        if i in lines_to_read:
            lines[i-41] = line.split()
    

#### We store the BELM spectra's data in arrays named as xB

    for i in range(len(lines)):
        eB[i] = lines[i][2]
        spectrumB[i] = lines[i][5]
        compB[i] = lines[i][6]
        synchB[i] = lines[i][7]
        
    
    plotBELM(nuLB,nuLs, synchB, nuLc, compB, nuL, spectrumB, tobs, e_pho, eB, 
             titre, ylab, xlab, ymin, sup, xmin, xsup)
    
def compareBELMnoCompton(nuLB,nuLs,nuL,tobs,e_pho,nu,
                titre, ylab,xlab,ymin,sup,xmin,xsup):

    """
    compare the results entered in arguments with the results get from the
    BELM model
    """
    firstline=41
    lastline=224
    nblines=184
    datafile = open("synch_escape.txt")
    lines_to_read = np.linspace(firstline,lastline,nblines)
    lines = np.empty((nblines,7))
    eB = np.empty(nblines)
    spectrumB = np.empty(nblines)
    synchB = np.empty(nblines)


    for i, line in enumerate(datafile):
    
        if i in lines_to_read:
            lines[i-41] = line.split()
    

#### We store the BELM spectra's data in arrays named as xB

    for i in range(len(lines)):
        eB[i] = lines[i][2]
        spectrumB[i] = lines[i][4]
        synchB[i] = lines[i][5]
        
    
    plotBELMnoCompton(nuLB,nuLs, synchB, nuL, spectrumB, tobs, e_pho, eB, 
             titre, ylab, xlab, ymin, sup, xmin, xsup)
    
    
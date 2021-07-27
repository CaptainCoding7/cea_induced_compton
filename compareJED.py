# -*- coding: utf-8 -*-
"""
Created on Tue Jul 20 11:42:48 2021

@author: pa267340
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from plotter import setFigureParameters
from plotter import noLogsetFigureParameters


def plotJED(nuLB, nuLs, nuLb, nuLc, nuLeq, spectrumB, tobs, e_pho, eJ,
                  titre, ylab,xlab,ymin,sup,xmin,xsup):
    """
    plot the sepctral radiance (intensité specifique) as a function of the energy
    """
    col = [ cm.jet(x) for x in np.linspace(0, 0.3 , len(tobs))]
    
    # Synchrotron only
    plt.plot(e_pho, nuLs, color='red',label='Synchrotron contribution')
    
    # Bremsstralhung only
    plt.plot(e_pho, nuLb, color='purple',label='Bremstralhung contribution')

    # compton only
    plt.plot(e_pho, nuLc, color='orange',label='Compton contribution')

    # JED Spectra
    plt.plot(eJ, spectrumB, color = 'blue',label='JED Spectrum')   

    # Blackbody emission
    plt.plot(e_pho,nuLB, color='black',label='Blackbody emission')
    
    # plot the intensity of the initial photon field
    #plt.plot(e_pho, nuL[0], color = 'black',label='t=0s')

    plt.plot(e_pho, nuLeq, color = 'green',label='Equilibrium spectrum')    
        
    setFigureParameters(titre, ylab,xlab,ymin,sup,xmin,xsup)
    


def compareJED(nuLB,nuLs,nuLb,nuLc,nuL,tobs,e_pho,nu,
                titre, ylab,xlab,ymin,sup,xmin,xsup):

    """
    compare the results entered in arguments with the results get from the
    BELM model
    """
    
    firstline=3
    lastline=514
    nblines=512
    datafile = open("SSC_Spec__nuFnu_Flux_99.txt")
    lines_to_read = np.linspace(firstline,lastline,nblines)
    lines = np.empty(nblines)
    
    keVfile = open("keV.txt")
    eJ = np.empty(nblines)
    for i, line in enumerate(keVfile):
        eJ[i] = line
    spectJED = np.empty(nblines)

    # distance au trou noir en cm
    D = 3.086e+22

    for i, line in enumerate(datafile):
    
        if i in lines_to_read:
            lines[i-firstline] = 4*np.pi*D**2*float(line)/10
    

#### We store the BELM spectra's data in arrays named as xB

    for i in range(len(lines)):
        spectJED[i] = lines[i]

    
    plotJED(nuLB,nuLs, nuLb, nuLc, nuL, spectJED, tobs, e_pho, eJ, 
             titre, ylab, xlab, ymin, sup, xmin, xsup)
    

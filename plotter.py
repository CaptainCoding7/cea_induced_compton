# -*- coding: utf-8 -*-
"""
Created on Tue Jun  1 16:26:31 2021

@author: Paul
"""
import matplotlib.pyplot as plt


def setFigureParameters(title, ylabel, xlabel, ymin, ymax,  xmin, xmax):
    
    # Reglages affichage
    plt.xscale('log')
    plt.yscale('log')
    leg = plt.legend(loc="upper right",prop={'size': 7}, bbox_to_anchor=[1, 1],
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
    #plt.xscale('log')
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



# plotting functions

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

def logplot(dt, x, u, n):
    

    plt.plot(x, u, label="t={:.2f}s".format(dt*float(n)))
       #print(u[n])
       #print(x)

    leg = plt.legend(loc="upper right", prop={'size': 8}, bbox_to_anchor=[1, 1],
                     ncol=2, shadow=True, title="Legend", fancybox=True)
    leg.get_title().set_color("black")

    #plt.plot(np.log10(data), pdf)
    plt.xscale('log')
    plt.yscale("log")
    plt.ylim(1e-8,1e3)
    #plt.xlim(1e-1, 0.5)
    plt.xlabel('Energy (kEv)')
    plt.ylabel('Particle Number Density')
    
    

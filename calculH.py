# -*- coding: utf-8 -*-
"""
Created on Tue Jul 20 13:25:07 2021

@author: pa267340
"""

from constants import k,me,cl,sT
import numpy as np

"""
model = "BSC spectre 91"
Te = 7.111e8
kTe = 61.2778
pT = 6.384
r = 5.3
e = 1.349e-1
"""
model = "BSC spectre 99"
Te = 1.255e9
kTe = 108.1475
pT = 6.824
r = 2.2
e = 1.934e-1
"""

model = "SSC spectre 99"
Te = 7.049e9
kTe = Te * k / 1.603E-16
pT = 3.266e-1
r = 2.2
e = 2.836e-1
"""

# calcul de H
G = 6.67e-11
M = 1.998e31
c = 3e8
Rg = G*M/c**2 *1e2
Rc = r*Rg # distance de la couronne au trou noir
H = e*Rc
R = H/2


# calcul de B

mu = 0.5
k_cgs = k * 1e7
me_cgs = me * 1e3
c_cgs = cl * 1e2
sT_cgs = sT * 1e4

theta = k_cgs*Te/me_cgs/c_cgs**2
print(theta)
lb = mu * pT * theta
print(lb)
Blb = np.sqrt(lb*8*np.pi*c_cgs**2*me_cgs / sT_cgs/R)

print("Les paramètres à rentrer pour le modèle", model, "sont:")
print('kT = {:.3E} keV'.format(kTe))
print('pT = {:.3E}'.format(pT))
print('R = {:.3E} cm'.format(R))
print("B = {:.3E}".format(Blb))

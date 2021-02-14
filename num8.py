# -*- coding: utf-8 -*-
"""
Created on Sun Feb 14 00:01:38 2021

@author: franc
"""
from scipy.stats import moyal
import matplotlib.pyplot as plt 
import numpy as np

#Densité de probabilité moyal avec maximum à 150
moyal = moyal(loc=150, scale=4)

#Énergie du proton incident
T = np.linspace(moyal.ppf(0.01), moyal.ppf(0.99), 10000)

#Densité de probabilité de l'énergie
p = moyal.pdf(T)

#Nombres aléatoires suivant la densité de probabilité moyal
r = moyal.rvs(size=10000)


plt.plot(T, p, "k--", label="Densité de probabilité moyal")
plt.hist(r, density=True, bins = 40, label="Variables aléatoires générées")
plt.xlabel("Énergie du proton incident [MeV]")
plt.ylabel("Densité de probabilité [-]")
plt.legend()
plt.title("Variables aléatoires générées suivant \n la distribution moyal")
plt.show()


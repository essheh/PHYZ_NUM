# -*- coding: utf-8 -*-
"""
Created on Sun Feb 14 00:01:38 2021

@author: franc
"""
import timeit
import numpy as np
import scipy
import sympy
from scipy.misc import derivative as der
import scipy.constants as cst
from numpy import sqrt, pi, log


################################ Num 5

eau_liq = {"Densité" : 1.00000E+00,
                       # Mean Excitation Energy [MeV]
           "MEE" : 75.000000E-06,
                               # numéro atomique
           "Composition" : np.array([[1, 8],
                                     # proportion dans le milieu
                                     [0.111894, 0.888106],
                                     # masse molaire [g/mol] (source ptable.com)
                                     [1.008, 15.999]])}
os_compact = {"Densité" : 1.85000E+00,
                      # Mean Excitation Energy [MeV]
              "MEE" : 91.900000E-06,
                                              # numéro atomique
              "Composition" : np.array([[1, 6, 7, 8, 12, 15, 16, 20],
                                                # proportion dans le milieu
                                        [0.063984, 0.278000, 0.027000, 0.410016, 0.002000, 0.070000, 0.002000, 0.147000],
                                                        # masse molaire [g/mol] (source ptable.com)
                                        [1.008, 12.011, 14.007, 15.999, 24.305, 30.738, 32.065, 40.078]])}
def n_e(num_ato, proportion, mass_mol, densité):
    return cst.Avogadro * densité * np.sum(proportion * num_ato / mass_mol)
eau_liq["n_e"] = n_e(num_ato = eau_liq["Composition"][0], proportion = eau_liq["Composition"][1],
                     mass_mol = eau_liq["Composition"][2], densité = eau_liq["Densité"])
os_compact["n_e"] = n_e(num_ato = os_compact["Composition"][0], proportion = os_compact["Composition"][1],
                        mass_mol = os_compact["Composition"][2], densité = os_compact["Densité"])
r_e = 2.8179403227E-13
m_e_c2 = 0.51099895000
m_p_c2 = 938.27208816
def gamma(T):
    return (T/(m_p_c2)) + 1
def beta(T):
    beta = sqrt((gamma(T)**2 - 1) / gamma(T)**2)
    return beta
a = 2 * m_e_c2
b = 1 + (m_e_c2 / m_p_c2)**2
d = 2*(m_e_c2 / m_p_c2)
def Te_max(T):
    return a*(gamma(T)**2 - 1) / (b + d * gamma(T))
def t1(T):
    return 2 * pi * r_e**2 *m_e_c2 / (beta(T)**2)
def t2(T):
    return 2 * m_e_c2 * beta(T)**2 * gamma(T)**2

def Scol(n_e, T, I):
    return t1(T) * n_e * (log(t2(T)*Te_max(T) / I**2) - 2*beta(T)**2)



#############################

#Num 8 
from scipy.stats import moyal
import matplotlib.pyplot as plt 
import numpy as np

#Densité de probabilité moyal avec maximum à 150
moyal = moyal(loc=150, scale=4)

#Énergie du proton incident
E = np.linspace(moyal.ppf(0.01), moyal.ppf(0.99), 10000)

#Densité de probabilité de l'énergie
p = moyal.pdf(E)

#10 000 valeurs d'énergie suivant la densité de probabilité moyal
T = moyal.rvs(size=10000)


plt.plot(E, p, "k--", label="Densité de probabilité moyal")
plt.hist(T, density=True, bins = 40, label="Variables aléatoires générées")
plt.xlabel("Énergie du proton incident [MeV]")
plt.ylabel("Densité de probabilité [-]")
plt.legend()
plt.title("Variables aléatoires générées suivant \n la distribution moyal")
plt.show()


#Num 9
#Fonction à intégrer tout au long du numéro 
def to_integrate(T):
    return 1/(Scol(eau_liq["n_e"], T, eau_liq["MEE"]) / eau_liq["Densité"])

#Temps d'intégration par scipy.integrate.quad


def test_1(T):
    result = np.zeros(10000)
    for i, value in enumerate(T):
         result[i] = scipy.integrate.quad(to_integrate, a = 3, b = value, limit = 20, epsabs = 10**-9)[0]
    return result
    
  
    
x = test_1(T)
plt.hist(x, bins = 100)
plt.ylabel("Nombre de protons")
plt.xlabel("Portée dans l'eau [cm]")
plt.show()

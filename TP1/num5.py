# Méthode des trapèzes.
# Méthode de Simpson.

## Num1
import numpy as np
import sympy
from sympy import *
import math
from scipy.misc import derivative as der
import scipy.constants as cst
from numpy import sqrt, pi, log


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

## Num5
JOULES_TO_EV = cst.value("joule-electron volt relationship")

def trapèze(a,b,func,N):
    h = (b - a) / N
    s = 0.5 * func(a) + 0.5 * func(b)
    for k in range(1,N):
        s += func(a+k*h)
    return s*h


def error(a,b,func,N0,precision):
    iter = 0
    i0 = trapèze(a,b,func,N0)
    #Deuxième itération
    N = 2 * N0
    iN = trapèze(a,b,func,N)
    err = abs((iN-i0)/3)
    while err >= precision:
        i = iN
        iter += 1
        N = 2 * N
        iN = trapèze(a,b,func,N)
        errN = abs((iN - i) / 3)
    return iN,err,N

# Valeurs initiales et bornes
a= 0.1
precision = 10**-16
N0=100000000
b = 1000000 * 150 / JOULES_TO_EV

#Calcul eau
def func_eau(T):
    Scol_eau = t1(T) * eau_liq["n_e"] * (log(t2(T)*Te_max(T) / eau_liq["MEE"]**2) - 2*beta(T)**2)
    return eau_liq["Densité"]/Scol_eau

i0 = trapèze(a,b,func_eau,N0)
    #Deuxième itération
N = 2 * N0
iN = trapèze(a,b,func_eau,N)
err = abs((iN-i0)/3)
# int_eau,err_eau,Neau = error(a,b,func_eau,N0,precision)
print(err)
import matplotlib.pylab as plt
import numpy as np
import timeit
from scipy import constants

hbar = constants.hbar #Planck's constant divided by 2*pi [m^2 kg / s]
M = constants.electron_mass #Mass of the electron [kg]
L = 1e-08 #Length of the box [m]
x_0 = L/2 # [m]
sigma = 1e-10 #[m^-1]
kappa = 5e-10 #[m^-1]
######Number 1#######
N = 1000 #Number of spatial slice 
a = N/L  #Distance between each slice    
h = 10**(-18)  #Step [s]
a1 = 1 + 1j * h * hbar / (2 * M * a ** 2)
a2 = - h * hbar * 1j / (4 * M * a ** 2)
b1 = 1 - 1j * h * hbar / (2 * M * a ** 2)
b2 = - a2

def psi_0(x):
    return np.exp(-(x-x_0)**2/2*sigma**2)*np.exp(1j*kappa*x) 


def initial_condition(N, L, h): 
    a = N/L
    x = np.linspace(0, L, N+1) 
    psi = np.array(psi_0(x)) #Initial condition  
    psi[0] = 0 
    psi[-1] = 0
    return psi 
    
    
def v(psi):
    v = np.zeros(len(psi), dtype = 'complex_')
    for i in range(len(psi)):
        if i == len(psi)-1:
            v[i] = b1 * psi[i] + b2 * (psi[i-1])
        else:
            v[i] = b1 * psi[i] + b2 * (psi[i+1] + psi[i-1]) 
    return v

psi = initial_condition(N, L, h)
v = v(psi)
print(v)
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
    v = np.zeros(len(psi), dtype = 'complex_') #Initialize the vector v 
    n = n = len(psi) - 1 #Last element of vector v 
    v[0] = b1 * psi[0] + b2 * (psi[1])
    v[1:n-1] = b1 * psi[1:n-1] + b2 * (psi[2:n] + psi[0:n-2])
    v[n] = b1 * psi[n] + b2 * (psi[n-1])
    
    return v

psi = initial_condition(N, L, h)



######Number 2#######
def A(N):
    A = np.zeros([3, N + 1], dtype = 'complex') #We write the tridiagonal matrix in the good form 
    A[0, 0] = 0
    A[0, 1:] = a2
    A[1, :] = a1
    A[2, 0:N-1] = a2
    A[2, N] = 0
    return A 
    
def Thomas(A, v, N): 
    x = np.zeros(N+1, dtype = 'complex')
    a_ = A[0, 0:] 
    b_ = A[1, 0:]
    c_ = A[2, 0:]
    d_ = v[:]
    c_prime = np.zeros(N + 1, dtype = 'complex' )
    d_prime = np.zeros(N + 1, dtype = 'complex' )
    
    c_prime[0] = c_[0]/b_[0]
    c_prime[1: N - 1] = c_[1: N - 1]/(b_[1: N - 1] - a_[1: N - 1] * c_prime[0: N - 2])
    d_prime[0] = d_[0]/b_[0] 
    d_prime[1: N] = (d_[1: N] - a_[1: N] * d_prime[0: N - 1 ])/(b_[1: N] - a_[1: N] * c_prime[0: N - 1])
    x[N] = d_[N]/b_[N]
    x[N-1: :-1] = (d_[N-1: :-1] - c_[N-1: :-1]*x[N:0 :-1])/b_[N-1: :-1] 
    
    
    
    return x 

A = A(N)
v = v(psi)
print(Thomas(A, v, N))
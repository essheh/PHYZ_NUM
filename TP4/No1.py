import matplotlib.pylab as plt
import numpy as np
import timeit
from scipy import constants

hbar = constants.hbar #Planck's constant divided by 2*pi [m^2 kg / s]
M = constants.electron_mass #Mass of the electron [kg]
L = 1e-08 #Length of the box [m]
x_0 = L/2 # [m]
sigma = 1e-10 #[m]
kappa = 5e10 #[m^-1]
######Number 1#######
N = 1000 #Number of spatial slice 
a = N/L  #Distance between each slice    
h = 10**(-18)  #S       tep [s]
a1 = 1 + 1j * h * hbar / (2 * M * a ** 2)
a2 = - h * hbar * 1j / (4 * M * a ** 2)
b1 = 1 - 1j * h * hbar / (2 * M * a ** 2)
b2 = - a2

def psi_0(x):
    return np.exp(-(x-x_0)**2/(2*sigma)**2)*np.exp(1j *kappa *x)    


def initial_condition(N, L, h): 
    a = N/L
    x = np.linspace(0, L, N+1) 
    psi = np.array(psi_0(x)) #Initial condition  
    psi[0] = 0 
    psi[-1] = 0 
    return psi 
    
    
def v(psi):
    v = np.zeros(len(psi), dtype = 'complex_') #Initialize the vector v 
    n = len(psi) - 1 #Last element of vector v 
    v[0] = b1 * psi[0] + b2 * (psi[1])
    v[1:n] = b1 * psi[1:n] + b2 * (psi[2:n+1] + psi[0:n-1])
    v[n] = b1 * psi[n] + b2 * (psi[n-1])
    
    return v





######Number 2#######
def A_tri(N):
    A = np.zeros([3, N + 1], dtype = 'complex') #We write the tridiagonal matrix in the good form 
    A[0, 0] = 0
    A[0, 1:] = a2
    A[1, :] = a1
    A[2, 0:N] = a2
    A[2, N] = 0
    return A 
    
def A(N):
    A = np.zeros([N+1, N+1], dtype = 'complex')
    for i in range(N+1):
        A[i, i] = a1
    for i in range(N):     
        A[i, i+1] = a2 
    for i in range(N):
        A[i+1, i] = a2  
    return A
    
def Thomas(A_tri, v, N): 
    c_prime = np.zeros(N, dtype = 'complex')
    d_prime = np.zeros(N + 1, dtype = 'complex')
    x = np.zeros(N + 1, dtype = 'complex')
    
    a_ = A_tri[0, :] 
    b_ = A_tri[1, :]
    c_ = A_tri[2, :]
    d_ = v[:]
    
    c_prime[0] = c_[0]/b_[0] #First element of the new c' vector
    d_prime[0] = d_[0]/b_[0] #First element of the new d' vector
    
    for i in range(1, N):   
        c_prime[i] = c_[i]/(b_[i] - a_[i] * c_prime[i - 1])
    
    for i in range(1, N+1): 
        d_prime[i] = (d_[i] - a_[i] * d_prime[i-1])/(b_[i] - a_[i] * c_prime[i-1])
    
    x[N] = d_prime[N]
    
    for i in range(N-1, -1, -1):
        x[i] = d_prime[i] - c_prime[i]*x[i+1] 
        
    return x

#psi = initial_condition(N, L, h)
#v = v(psi)  
#A_tri = A_tri(N)    
#x = Thomas(A_tri, v, N)
#A = A(N)
#Ax = np.matmul(A, x)
#print(x) 
#print(Ax)
#print(v)

#comparison = Ax == v
#equal_arrays = comparison.all()
#print(equal_arrays)


#c_prime[0] = c_[0]/b_[0]
#c_prime[1: N] = c_[1: N]/(b_[1: N] - a_[1: N] * c_prime[0: N - 1])
#d_prime[0] = d_[0]/b_[0] 
#d_prime[1: N+1] = (d_[1: N+1] - a_[1: N+1] * d_prime[0: N])/(b_[1: N+1] - a_[1: N+1] * c_prime[0: N])
#x[N] = d_prime[N]
#x[N-1: :-1] = (d_prime[N-1: :-1] - c_prime[N-1: :-1]*x[N:0 :-1])

#Numéro 3 
def psi(N, L, h, t):
    psi = initial_condition(N, L, h) #Initial wave function at time t = 0
    x = np.linspace(0, L, N+1) 
    total_steps = int(t//h)  #Number of steps to reach time t
    v_vector = np.zeros(len(psi))
    #solution =  [psi]
    for i in range(total_steps):
        v_vector = v(psi)
        psi = Thomas(A_tri(N), v_vector, N) #New wave function at time t + h  
        #solution.append(psi) 
    return x, psi  

t = 100*h     
x, psi = psi(N, L, h, t)
plt.plot(x, abs(psi)**2) 
plt.show()
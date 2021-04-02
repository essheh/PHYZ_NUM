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
a = L/N  #Distance between each slice    
h = 1e-18  #Step [s]
a1 = 1 + 1j * h * hbar / (2 * M * a ** 2)
a2 = - h * hbar * 1j / (4 * M * a ** 2)
b1 = 1 - 1j * h * hbar / (2 * M * a ** 2)
b2 = - a2

#Parameters of the simulation

N = 1000 #Number of spatial slice 




def psi_0(x): #Return the initial wave function psi (x,0) at time t = 0 
    return np.exp(-(x-x_0)**2/(2*sigma**2))*np.exp(1j *kappa *x)  


def initial_condition(N, L, h): #Set a wave function at time t = 0, with each point in x separated by a
    x = np.linspace(0, L, N+1) #Each point separated by a
    psi = np.array(psi_0(x)) #Initial condition  
    psi[0] = 0 #The wave function equal 0 at discontinuity point i.e. infinite potential value 
    psi[-1] = 0 #Same as last 
    return psi 
    
    
x = np.linspace(0, L, N+1) #Each point separated by a
psi_t0 = initial_condition(N, L, h)
plt.plot(x, abs(psi_t0)**2)
#plt.show()
    
    
def v(psi): #return the vector v
    b1 = 1 - 1j * h * hbar / (2 * M * a ** 2)
    b2 = h * hbar * 1j / (4 * M * a ** 2)
    
    v = np.zeros(len(psi), dtype = complex) #Initialize the vector v 
    n = len(psi) - 1 #Last element of vector v 
    
    v[0] = b1 * psi[0] + b2 * (psi[1])
    v[1:n] = b1 * psi[1:n] + b2 * (psi[2:n+1] + psi[0:n-1])
    v[n] = b1 * psi[n] + b2 * (psi[n-1])
    
    return v


v_0 = v(psi_t0)
#print(v_0)



######Number 2#######
def A_tri(N):
    a1 = 1 + 1j * h * hbar / (2 * M * a ** 2)
    a2 = - h * hbar * 1j / (4 * M * a ** 2)
    A = np.zeros([3, N + 1], dtype = complex) #We write the tridiagonal matrix in the good form 
    
    A[0, 0] = 0
    A[0, 1:] = a2
    A[1, :] = a1
    A[2, 0:N] = a2
    A[2, N] = 0
    return A 

M1 = A_tri(N)  
    
def A_normal(N):
    A = np.zeros([N+1, N+1], dtype = complex)
    for i in range(N+1):
        A[i, i] = a1
    for i in range(N):     
        A[i, i+1] = a2 
    for i in range(N):
        A[i+1, i] = a2  
    return A
    
def Thomas(A_tri, v): 
    """
    The function returns the vector x of the system of equation Ax = v 
    param 1 A: is a triadiagonal matrix written in the form a_{i} * x_{i-1} + b_{i} * x_{i}+ c_{i} * x_{i+1}
    param 2 v: v is the vector of the RHS    
    return: x, the solution of the system  
    """
    N = len(v) - 1 
    c_prime = np.zeros(N, dtype = complex) #Initialization of the vector c' 
    d_prime = np.zeros(N + 1, dtype = complex) #Initialization of the vector d' 
    x = np.zeros(N + 1, dtype = complex) #Initialization of the vector x 
    
    a_ = A_tri[0, :] #Lower diagonal
    b_ = A_tri[1, :] #Middle diagonal
    c_ = A_tri[2, :] #Upper diagonal
    d_ = v[:] #RHS 
    
    c_prime[0] = c_[0]/b_[0] #First element of the new c' vector
    d_prime[0] = d_[0]/b_[0] #First element of the new d' vector
    
    for i in range(1, N):   
        c_prime[i] = c_[i]/(b_[i] - a_[i] * c_prime[i - 1]) #New values of c' 
    
    for i in range(1, N+1): 
        d_prime[i] = (d_[i] - a_[i] * d_prime[i-1])/(b_[i] - a_[i] * c_prime[i-1]) #New values of d'
    
    x[N] = d_prime[N] #Last element of the vector x 
    
    for i in range(N-1, -1, -1):
        x[i] = d_prime[i] - c_prime[i]*x[i+1] #New values of x by substitution of c' and d', from the before last to the first element

        
    return x

#First we must write the matrix A in its matrix form 

def A(N):
    
    A = np.zeros([N+1, N+1], dtype = complex) #We initialize the size of the square matrix
    a1 = 1 + 1j * h * hbar / (2 * M * a ** 2)
    a2 = - h * hbar * 1j / (4 * M * a ** 2)
    
    for i in range(N+1): #Middle diagonal
        A[i, i] = a1
    
    for i in range(N):   #Upper diagonal  
        A[i, i+1] = a2 
    
    for i in range(N):   #Lower diagonal
        A[i+1, i] = a2  
    
    return A

M2 = A(N)

#We check if the matrix product Ax = v 

x = Thomas(M1, v_0) #The vector x found by Thomas algorithm 
Ax = np.matmul(M2, x) #The matrix product

comparison = Ax == v_0
equal_arrays = comparison.all() #If the value of equal_arrays is true then all the elements of Ax and v are the same 


if equal_arrays:
    print("Le vecteur de la multiplication matricielle Ax est égale au vecteur v.")





#M2 = A_normal(N)
#Ax = np.matmul(M2, x)
#comparison = Ax == v_vector
#equal_arrays = comparison.all()
#print(equal_arrays)


#c_prime[0] = c_[0]/b_[0]
#c_prime[1: N] = c_[1: N]/(b_[1: N] - a_[1: N] * c_prime[0: N - 1])
#d_prime[0] = d_[0]/b_[0] 
#d_prime[1: N+1] = (d_[1: N+1] - a_[1: N+1] * d_prime[0: N])/(b_[1: N+1] - a_[1: N+1] * c_prime[0: N])
#x[N] = d_prime[N]
#x[N-1: :-1] = (d_prime[N-1: :-1] - c_prime[N-1: :-1]*x[N:0 :-1])

def psi(N, L, h, t):
    A = A_tri(N) 
    psi = initial_condition(N, L, h) #Initial wave function at time t = 0
    x = np.linspace(0, L, N+1) 
    total_steps = int(t//h)  #Number of steps to reach time t
    solution =  [psi]
    for i in range(total_steps):
        v_vector = v(psi)
        psi = Thomas(A, v_vector) #New wave function at time t + h  
    
        solution.append(psi) 
    return x, solution 

t = 1200*h     
x, solution = psi(N, L, h, t)
def psi(N, L, h, t):
    A = A_tri(N) 
    psi = initial_condition(N, L, h) #Initial wave function at time t = 0
    x = np.linspace(0, L, N+1) 
    total_steps = int(t//h)  #Number of steps to reach time t
    solution =  [psi]
    for i in range(total_steps):
        v_vector = v(psi)
        psi = Thomas(A, v_vector) #New wave function at time t + h  
    
        solution.append(psi) 
    return x, solution 

t = 1600*h     
x, solution = psi(N, L, h, t)
for i in range(16):
    plt.plot(x, abs(solution[i*100])**2) 

plt.show()

def psi_linlag(N, L, h, t):
    A = A(N) 
    psi = initial_condition(N, L, h) #Initial wave function at time t = 0
    x = np.linspace(0, L, N+1) 
    total_steps = int(t//h)  #Number of steps to reach time t
    solution =  [psi]
    for i in range(total_steps):
        v_vector = v(psi)
        psi = np.linalg.solve(A, v_vector) #New wave function at time t + h  
        solution.append(psi) 
    return x, solution 

t = 1600*h     
x, solution = psi(N, L, h, t)
for i in range(16):
    plt.plot(x, abs(solution[i*100])**2) 
plt.show()

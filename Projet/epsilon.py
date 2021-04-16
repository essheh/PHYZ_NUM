import numpy as np
import mpmath
from mpmath import *
import matplotlib.pyplot as plt

def leibniz(n):
    terms_leibniz = []
    part_sum_leibniz = []
    t_sum = 0
    for i in range(n):
        term = (-1) ** i * 4/(2*i+1)
        terms_leibniz.append(term)
        t_sum = t_sum + term
        part_sum_leibniz.append(t_sum)
  
# Erreur relative en pourcentage (pct) selon la valeur de pi de numpy
    err_pi_pct = 100*abs(np.pi - part_sum_leibniz[n-1]) / np.pi    
    
    return terms_leibniz, part_sum_leibniz, err_pi_pct


def ln3(n):
    terms_ln3 = []
    part_sum_ln3 = []
    t_sum = 0
    for i in range(n):
        term = ((-1)**(i)) * 2**(i+1) / (i+1)
        terms_ln3.append(term)
        t_sum = t_sum + term
        part_sum_ln3.append(t_sum)
        
# Erreur relative en pourcentage (pct) selon la valeur de pi de numpy
    err_ln3_pct = 100*abs(np.log(3) - part_sum_ln3[n-1]) / np.log(3)    
    
    return terms_ln3, part_sum_ln3, err_ln3_pct

def epsilon(partial):
    """
    The function returns the Sum of n terms of a serie with the epsilon algorithm.
    param 1 S: array of terms of a serie with the size n. 
    return: sum, Sum of the serie.   
    """
    e = np.zeros((n + 1, n + 1))

    for i in range(1, n + 1):
        e[i, 1] = partial[i - 1]

    for i in range(3, n + 2):
        for j in range(3, i + 1):
            e[i - 1, j - 1] = e[i - 2, j - 3] + 1 / (e[i - 1, j - 2] - e[i - 2, j - 2])

    sumation = e[:, 1:n + 1:2]
    return sumation[-1,-1]

k = 100
sum_leibniz = []
sum_ln3 = []
epsi_leibniz = []
epsi_ln3 = []

for N in range(1,k):
    terms_leibniz, part_sum_leibniz, err_pi = leibniz(N)
    sum_leibniz.append(part_sum_leibniz[N-1])
    epsi = epsilon(part_sum_leibniz,N)
    epsi_leibniz.append(epsi)
    
for N in range(1,k):
    terms_ln3, part_sum_ln3, err_ln3 = ln3(N)
    sum_ln3.append(part_sum_ln3[N-1])
    epsi = epsilon(part_sum_ln3,N)
    epsi_ln3.append(epsi)
    
f = list(range(1, k))

plt.plot(f,epsi_leibniz, label = 'epsilon')
plt.plot(f,sum_leibniz,label = 'somme')
plt.ylim(2.5,5)
plt.legend()
plt.show()

plt.plot(f,epsi_ln3, label = 'epsilon')
plt.plot(f,sum_ln3,label = 'somme')
plt.ylim(-100,100)
plt.xlim(0,10)
plt.legend()
plt.show()




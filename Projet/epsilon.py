import numpy as np
import mpmath
from mpmath import *
import matplotlib.pyplot as plt
# s=0
# ps=[]
# n=20
# for k in range(0, n):
#     term = ((-1)**(k)) * 2**(k+1) / (k+1)
#     s += term
#     ps.append(s)
    
def epsilon(S,n):
    """
    The function returns the Sum of n terms of a serie with the epsilon algorithm.
    param 1 S: array of terms of a serie with the size n. 
    return: sum, Sum of the serie.   
    """
    partials= []
    s=0
    for k in range(0, n):
        term = S(k)
        s += term
        partials.append(s)
    e = np.zeros((n + 1, n + 1))

    for i in range(1, n + 1):
        e[i, 1] = partials[i - 1]

    for i in range(3, n + 2):
        for j in range(3, i + 1):
            e[i - 1, j - 1] = e[i - 2, j - 3] + 1 / (e[i - 1, j - 2] - e[i - 2, j - 2])

    sumation = e[:, 1:n + 1:2]
    return sumation[-1,-1]


leibnizS = lambda i: (-1) ** i * 4/(2*i+1)
ln3S = lambda i: ((-1)**(i)) * 2**(i+1) / (i+1)

leibniz_epsi = np.zeros(49)

for n in  range(2,50):
    print(n)
    leibniz_epsi[n-1] = epsilon(leibnizS,n)
    
N = range(1,50)
print(N)



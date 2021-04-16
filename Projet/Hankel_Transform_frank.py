# -*- coding: utf-8 -*-
"""
Created on Wed Apr 14 19:49:02 2021

@author: franc
"""
import scipy.integrate
import numpy as np
import scipy.special as sp
import matplotlib.pyplot as plt

def epsilon(terms):
    """
    The function returns the Sum of n terms of a serie with the epsilon algorithm.
    param 1 S: array of terms of a serie with the size n. 
    return: sum, Sum of the serie.   
    """
    n = len(terms)
    e = np.zeros((n + 1, n + 1))

    for i in range(1, n + 1):
        e[i, 1] = terms[i - 1]

    for i in range(3, n + 2):
        for j in range(3, i + 1):
            e[i - 1, j - 1] = e[i - 2, j - 3] + 1 / (e[i - 1, j - 2] - e[i - 2, j - 2])

    sumation = e[:, 1:n + 1:2]
    return sumation[-1,-1]



def simpson(func, a, b, N):
    """
    Integration by Simpson method
    param 1 func: func is the function to integrate 
    param 2 a: a is the lower boundary of the integral 
    param 3 b: b is the upper boundary of the integral  
    param 4 N: N is the number of slic
    return: I, the integration of the function   
    """
    h = abs(b - a) / N  # Size of the slices 
    x = np.linspace(a, b, N+1) #N+1 points in the interval [a,b]
    y = func(x)
    
    y_odd = y[1:-1:2] #Sum of odd 
    s1 = np.sum(y_odd)
    
    y_even = y[2:-2:2] #Sum of even
    s2 = np.sum(y_even)
    
    I = h/3*(y[0] + y[-1] + s1*4+s2*2) #Area under the curve 
    
    return I 
    
def Hankel_transform(g, p, n, L):
    """
    Parameters
    ----------
    g : The function to transform
    p : The p value in the Fourrier-Bessel space 
    n : The order of the Bessel function of the first kind (Integer)
    L : The number of term for the partial sum (Integer)
    Returns
    -------
    The partial sum L of the Hankel transform of the function g
    """
    zeros = sp.jn_zeros(n, L+2) # array of the L+1 first zeros of bessel order n  
    value = 0.0 #Initialize the value 
    partial_sums = []
    #function to integrate
    def f(x):
        return x*g(x/p)/p**2 * sp.jv(n, x) 
    
    for i in range(L+1):
       #I = simpson(f, zeros[i], zeros[i+1], 1000)
       I = scipy.integrate.quad(f, zeros[i], zeros[i+1])[0]
       value += I
       partial_sums.append(value)
    epsilon_value = epsilon(partial_sums)
    return value, partial_sums, epsilon_value

def f1(x): 
    return x
p = 50
n = 0
L = 10

p_value, partial_sums, epsilon_value = Hankel_transform(f1, p, n, L)



F = []
for i in range(1, p+1):
    F.append(Hankel_transform(f1, i, n, L)[2])

x = np.arange(1, 51, 1)
plt.plot(x, F)
plt.plot(x, -1/x**3)
plt.show()

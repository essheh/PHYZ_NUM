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
    
def Hankel_transform(func, p, n, L):
    """
    Parameters
    ----------
    g : function that receives one argument and returns one 
        g the function to transform
    p :The p value in the Fourrier-Bessel space 
    n : integer
        The order of the Bessel function of the first kind.
    L : The number of term for the partial sum 
    Returns
    -------
    The partial sum L of the hankel transform of the function func 
    """
    zeros = sp.jn_zeros(n, L+1) # array of the L+1 first zeros of bessel order n  
    value = 0 #Initialize the value 
    terms = []
    #function to integrate
    def f(x):
        return x*func(x/p)/p**2 * sp.jv(n, x) 
    
    for i in range(L):
       #I = simpson(f, zeros[i], zeros[i+1], 10000)
       I = scipy.integrate.quad(f, zeros[i], zeros[i+1])[0]
       value += I
       terms.append(I)
    epsilon_value = epsilon(terms)
    return value, terms, epsilon_value

def f1(x): 
    return -2*x
p = 1
n = 0
L = 50

p_value, terms, epsilon_value = Hankel_transform(f1, p, n, L)

print(epsilon_value)
n = np.arange(len(terms))
plt.plot(n, terms)

def ln3(n):
    terms_ln3 = []
    t_sum = 0
    for i in range(n):
        term = ((-1)**(i)) * 2**(i+1) / (i+1)
        terms_ln3.append(term)
        t_sum = t_sum + term
        part_sum_ln3 = t_sum
        
# Erreur relative en pourcentage (pct) selon la valeur de pi de numpy
    err_ln3_pct = 100*abs(np.log(3) - part_sum_ln3) / np.log(3)    
    
    return terms_ln3, part_sum_ln3, err_ln3_pct

n = 100

terms_ln3, part_sum_ln3, err_ln3 = ln3(n)

print(epsilon(terms_ln3))